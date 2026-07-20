// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "SimulationEngine.h"
#include "core/utils/safety/numericalSafetyLogger.h"
#include "core/utils/Profiler.h"
#include "core/utils/RngUtils.h"
#include "core/log/Logger.h"
#include "core/types/IonEnsemble.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include <cmath>
#include <algorithm>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <chrono>
#include <cstdlib>  // for setenv (NUMA thread placement)
#include <limits>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ICARION {
namespace integrator {

SimulationEngine::SimulationEngine(
    const config::FullConfig& config,
    std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries,
    std::shared_ptr<IIntegrationStrategy> integrator,
    std::shared_ptr<physics::ICollisionHandler> collision_handler,
    std::shared_ptr<physics::IReactionHandler> reaction_handler
) : config_(config),
    force_registries_(std::move(force_registries)),
    integrator_(integrator),
    collision_handler_(collision_handler),
    reaction_handler_(reaction_handler),
    deep_collision_diagnostics_(config_.output)
{
    // Validate force registries (must have one per domain)
    if (force_registries_.size() != config_.domains.size()) {
        throw std::invalid_argument(
            "SimulationEngine: force_registries size (" + std::to_string(force_registries_.size()) + 
            ") must match domains size (" + std::to_string(config_.domains.size()) + ")"
        );
    }
    
    for (size_t i = 0; i < force_registries_.size(); ++i) {
        if (!force_registries_[i]) {
            throw std::invalid_argument(
                "SimulationEngine: ForceRegistry[" + std::to_string(i) + "] cannot be null"
            );
        }
        // Rebind registries to the engine-local config copy (SSOT at runtime).
        force_registries_[i]->set_domain(&config_.domains[i]);
    }
    
    if (!integrator_) {
        throw std::invalid_argument("SimulationEngine: IntegrationStrategy cannot be null");
    }

    parallel_enabled_ = config_.simulation.enable_openmp;

    // Capture RK45 runtime settings for metadata
    if (auto* rk45 = dynamic_cast<RK45Strategy*>(integrator_.get())) {
        auto cfg = rk45->get_config();
        config_.simulation.rk45_runtime_settings = config::SimulationConfig::RK45RuntimeSettings{
            cfg.atol,
            cfg.rtol,
            cfg.safety_factor,
            cfg.min_step_factor,
            cfg.max_step_factor,
            cfg.max_step_increase,
            cfg.max_step_decrease,
            cfg.absolute_min_step_s
        };
        rk45->enable_stats(!parallel_enabled_);
    }

    integrator_->set_parallel_enabled(parallel_enabled_);
    
    // Create domain manager
    domain_manager_ = std::make_unique<DomainManager>(
        config_.domains,
        config_.simulation.rng_seed
    );
    // Wire field models into force registries (same ordering as domains) if not set
    for (size_t i = 0; i < force_registries_.size(); ++i) {
        if (force_registries_[i] && force_registries_[i]->field_model() == nullptr) {
            force_registries_[i]->set_field_model(domain_manager_->field_model(static_cast<int>(i)));
        }
    }
    
    // Initialize OpenMP thread settings (NUMA-aware)
    initialize_openmp_settings();
    
    // Create output manager
    std::string hdf5_path = config_.output.folder + "/" + config_.output.trajectory_file;
    std::string log_path = "";  // TODO(v1.1): Add log file path to OutputConfig
    if (config_.output.print_progress) {
        log_path = config_.output.folder + "/simulation.log";
    }
    
    output_manager_ = std::make_unique<OutputManager>(
        hdf5_path,
        log_path,
        config_.simulation.dt_s * config_.simulation.write_interval,  // Write interval in seconds
        50  // Buffer max
    );
    if (config_.output.buffer_byte_cap > 0) {
        output_manager_->set_buffer_byte_cap(config_.output.buffer_byte_cap);
        output_manager_->log_progress("Output buffer byte cap set to " + std::to_string(config_.output.buffer_byte_cap) + " bytes");
    }

    if (!parallel_enabled_ && integrator_ && integrator_->is_adaptive()) {
        output_manager_->log_progress("OpenMP disabled for adaptive integrator to avoid shared state races (RK45 not thread-safe)");
    }
    
    // Initialize numerical safety logger if enabled
    if (config_.simulation.enable_safety_logging) {
        std::string safety_log = config_.output.folder + "/numerical_safety.log";
        safety::NumericalSafetyLogger::getInstance(safety_log, config_.simulation.verbose_safety);
        
        // Configure logger (enable_logging, verbose_mode, buffer_size, max_history)
        safety::NumericalSafetyLogger::getInstance().configure(
            true,  // enable_logging
            config_.simulation.verbose_safety,  // verbose_mode
            1000,  // buffer_size (default)
            10000  // max_history (default)
        );
        
        output_manager_->log_progress("Numerical safety logging enabled: " + safety_log);
    }
}



void SimulationEngine::initialize(const core::IonEnsemble& ensemble) {
    // Initialize output system without per-step AoS conversions
    output_manager_->initialize(config_, ensemble);
    
    output_manager_->log_progress("Simulation engine initialized (SoA)");
    
    std::ostringstream msg;
    msg << "Configuration: " << ensemble.size() << " ions, "
        << config_.domains.size() << " domains, "
        << "dt = " << config_.simulation.dt_s * 1e9 << " ns, "
        << "t_max = " << config_.simulation.total_time_s * 1e6 << " µs";
    output_manager_->log_progress(msg.str());
    
#ifdef ICARION_USE_GPU
    initialize_gpu(config_.simulation.enable_gpu);
#endif

    // Warn about space charge update cadence
    for (const auto& registry : force_registries_) {
        if (registry && registry->space_charge_model()) {
            space_charge_stale_warned_ = true;
            if (integrator_ && integrator_->is_adaptive()) {
                if (!adaptive_space_charge_enabled()) {
                    output_manager_->log_progress("Warning: Space charge + adaptive RK45 disabled via ICARION_ADAPTIVE_SC=0; fields updated once per macro-step.");
                } else {
                    output_manager_->log_progress("Info: Space charge + adaptive RK45 will rebuild fields at each RK stage (performance heavy).");
                }
            } else {
                output_manager_->log_progress("Warning: Space-charge fields are updated once per macro-step; fast-changing clouds may be inaccurate unless using RK4/RK45 stage-refresh paths.");
            }
            break;
        }
    }

    // Warn about environment cache approximation for collisions/reactions
    if (config_.domains.size() > 1 &&
        config_.physics.collision_model != config::CollisionModel::NoCollisions) {
        output_manager_->log_progress(
            "Warning: Collisions/reactions use pre-step environment; mid-step boundary crossings do not refresh gas properties (approximate SSOT). Use small dt across sharp domain changes.");
    }

    if (config_.physics.collision_multi_event_mode) {
        output_manager_->log_progress(
            "Info: physics.collision_multi_event_mode enabled. Collisions use approximate micro-subcycling; for highest accuracy keep unresolved multiple collisions rare.");
    }
}

void SimulationEngine::initialize_openmp_settings() {
#ifdef _OPENMP
    if (!parallel_enabled_) {
        return;  // OpenMP disabled
    }
    
    // Use all available threads (respects OMP_NUM_THREADS environment variable)
    int num_threads = omp_get_max_threads();
    omp_set_num_threads(num_threads);
    
    // NUMA-aware thread placement (Linux only)
    #ifdef __linux__
    // Check if user already set OMP_PLACES
    const char* omp_places = std::getenv("OMP_PLACES");
    if (!omp_places) {
        // Set default: bind threads to physical cores (prevents migration)
        // OMP_PLACES=cores: One thread per physical core
        // OMP_PROC_BIND=close: Bind threads to nearby cores (NUMA-aware)
        setenv("OMP_PLACES", "cores", 0);
        setenv("OMP_PROC_BIND", "close", 0);
        
        // Note: Performance impact depends on system topology:
        // - NUMA systems (AMD Threadripper, EPYC): measured +20-30% speedup in internal runs
        // - Uniform memory (single-socket): measured +5-10% speedup
        // - Small systems (laptop): minimal impact
    }
    #endif
    
#endif  // _OPENMP
}

void SimulationEngine::log_progress(double t) {
    const int total_steps = static_cast<int>(config_.simulation.total_time_s / config_.simulation.dt_s);
    const int log_interval = std::max(1, total_steps / 10);
    
    if (current_step_ % log_interval == 0 || current_step_ == total_steps - 1) {
        double percent = 100.0 * t / config_.simulation.total_time_s;
        
        std::ostringstream msg;
        msg << std::fixed << std::setprecision(0)
            << percent << "% completed (t = " << t * 1e3 << " ms, step " 
            << current_step_ << "/" << total_steps << ")";
        
        output_manager_->log_progress(msg.str());
    }
}

core::IonEnsemble SimulationEngine::run(core::IonEnsemble& ensemble) {
    const auto wall_start = std::chrono::steady_clock::now();
    
    // 1. Initialize subsystems
    initialize(ensemble);
    
    // 2. Main time loop using SoA
    current_time_ = 0.0;
    current_step_ = 0;
    dt_per_ion_.assign(ensemble.size(), config_.simulation.dt_s);
    deep_collision_diagnostics_.reset(ensemble.size());
    collision_runtime_stats_.reset();
    
    output_manager_->log_progress("Starting main simulation loop (SoA)");
    
    // Helper to check if we should continue
    auto should_continue = [&]() -> bool {
        if (current_time_ >= config_.simulation.total_time_s) {
            return false;
        }
        const auto* born_flags = ensemble.born_data();
        for (size_t i = 0; i < ensemble.size(); ++i) {
            if (ensemble.is_active(i)) {
                return true;
            }
            // Continue looping if a birth is scheduled in the future
            if (!born_flags[i] && ensemble.birth_time(i) >= current_time_ &&
                ensemble.birth_time(i) <= config_.simulation.total_time_s) {
                return true;
            }
        }
        return false;
    };
    
    while (should_continue()) {
        // Process one timestep using SoA
        double new_time = process_timestep(ensemble);
        
        // OutputManager owns the time-based sampling schedule. Combining this
        // with a separate step-modulo trigger produces adjacent duplicate
        // snapshots at every write interval.
        if (current_step_ == 0 || output_manager_->should_write(new_time)) {
            PROFILE_SCOPE_IF_ENABLED("Output Writing");
            output_manager_->log_step(new_time, ensemble);
        }
        
        // Update time and step counter
        current_time_ = new_time;
        current_step_++;
        
        // Progress logging (every 10%)
        {
            PROFILE_SCOPE_IF_ENABLED("Progress Update");
            log_progress(current_time_);
        }
    }
    
    // 3. Finalization
    if (integrator_->is_adaptive()) {
        if (auto* rk45 = dynamic_cast<RK45Strategy*>(integrator_.get())) {
            const auto& stats = rk45->get_stats();
            double avg_dt = (stats.accepted_steps > 0)
                ? stats.sum_step_used / static_cast<double>(stats.accepted_steps)
                : 0.0;
            double min_dt = (stats.min_step_used == std::numeric_limits<double>::max())
                ? 0.0
                : stats.min_step_used;
            std::ostringstream rk_msg;
            rk_msg << "RK45 stats: accepted=" << stats.accepted_steps
                   << ", rejected=" << stats.rejected_steps
                   << ", dt_avg=" << avg_dt
                   << ", dt_min=" << min_dt
                   << ", dt_max=" << stats.max_step_used
                   << ", avg_error=" << stats.avg_error;
            output_manager_->log_progress(rk_msg.str());
        }
    }
    output_manager_->log_progress("Simulation completed (SoA)");
    
    // Count active ions
    size_t active_count = 0;
    for (size_t i = 0; i < ensemble.size(); ++i) {
        if (ensemble.is_active(i)) {
            active_count++;
        }
    }
    
    std::ostringstream msg;
    msg << "Final state: " << active_count << "/" << ensemble.size() << " ions active";
    output_manager_->log_progress(msg.str());

    log_collision_runtime_stats();

    const auto wall_end = std::chrono::steady_clock::now();
    last_wall_runtime_s_ = std::chrono::duration<double>(wall_end - wall_start).count();
    
    // Direct SoA finalization (no conversion overhead)
    output_manager_->finalize(current_time_, ensemble);
    if (deep_collision_diagnostics_.enabled()) {
        const std::string hdf5_path = config_.output.folder + "/" + config_.output.trajectory_file;
        deep_collision_diagnostics_.write_hdf5(hdf5_path);
    }
    
    // Safety report
    if (config_.simulation.enable_safety_logging) {
        std::string report_file = config_.output.folder + "/numerical_safety_report.txt";
        safety::NumericalSafetyLogger::getInstance().generateSafetyReport(report_file);
        
        auto stats = safety::NumericalSafetyLogger::getInstance().getStatistics();
        std::ostringstream safety_msg;
        safety_msg << "Numerical safety: " << stats.total_violations << " violations detected";
        if (stats.recovery_attempts > 0) {
            safety_msg << ", " << stats.successful_recoveries << "/" 
                      << stats.recovery_attempts << " recoveries successful";
        }
        output_manager_->log_progress(safety_msg.str());
        output_manager_->log_progress("Safety report written: " + report_file);
    }
    
    // Return final SoA ensemble
    return ensemble;
}

void SimulationEngine::log_collision_runtime_stats() {
    if (!collision_runtime_stats_.has_activity()) {
        return;
    }

    const std::string coll_msg = collision_runtime_stats_.summary_message(current_step_);
    output_manager_->log_progress(coll_msg);
    log::Logger::main()->info("{}", coll_msg);

    if (!config_.physics.collision_multi_event_mode &&
        collision_runtime_stats_.should_warn_single_collision_timestep_load()) {
        const std::string warn_msg = collision_runtime_stats_.single_collision_timestep_warning();
        output_manager_->log_progress(warn_msg);
        log::Logger::main()->warn("{}", warn_msg);
    }
}

double SimulationEngine::process_timestep(core::IonEnsemble& ensemble) {
    const size_t n_ions = ensemble.size();
    if (dt_per_ion_.size() != n_ions) {
        dt_per_ion_.assign(n_ions, config_.simulation.dt_s);
    }

    update_dynamic_environments(current_time_);
    refresh_ensemble_environment_cache(ensemble);

    if (has_space_charge_model() && integrator_->is_adaptive() && !adaptive_space_charge_enabled()) {
        throw std::runtime_error("Space charge with adaptive integrator (RK45) disabled via ICARION_ADAPTIVE_SC=0.");
    }

    std::vector<double> dt_used_per_ion = dt_per_ion_;
    std::vector<double> dt_next_per_ion = dt_per_ion_;

    if (rng_by_ion_.size() != n_ions) {
        PROFILE_SCOPE_IF_ENABLED("RNG Initialization");
        utils::sync_rng_pool_for_ensemble(
            rng_by_ion_, rng_fingerprints_, ensemble, config_.simulation.rng_seed);
    }

    update_space_charge_models(ensemble);

    const std::vector<int> integration_domains = prepare_ions_for_integration(ensemble);

    double max_dt_used = 0.0;
    {
        PROFILE_SCOPE_IF_ENABLED("Integrator (Total)");
        max_dt_used = perform_integration(
            ensemble,
            current_time_,
            dt_per_ion_,
            integration_domains,
            dt_used_per_ion,
            dt_next_per_ion);
    }

    // Apply stochastic processes using the actual dt used by the integrator
    perform_collisions(ensemble, dt_used_per_ion, integration_domains);
    perform_reactions(ensemble, dt_used_per_ion, integration_domains);

    finalize_ions_after_integration(ensemble, dt_used_per_ion, dt_next_per_ion);

    dt_per_ion_.swap(dt_next_per_ion);

    const double new_time = next_engine_time_after_step(
        ensemble,
        current_time_,
        config_.simulation.total_time_s,
        config_.simulation.dt_s,
        max_dt_used);
    rng_fingerprints_ = utils::build_ion_rng_fingerprints(ensemble);
    return new_time;
}

}  // namespace integrator
}  // namespace ICARION

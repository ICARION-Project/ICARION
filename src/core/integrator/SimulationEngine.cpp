// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "SimulationEngine.h"
#include "core/utils/safety/numericalSafetyGuards.h"
#include "core/utils/safety/numericalSafetyLogger.h"
#include "core/utils/Profiler.h"
#include "core/utils/RngUtils.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/log/Logger.h"
#include "core/types/IonEnsemble.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/integrator/strategies/GPUIntegrationStrategy.h"
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
#include <atomic>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef ICARION_VERSION
#define ICARION_VERSION "unknown"
#endif

namespace ICARION {
namespace integrator {

namespace {
constexpr int kOmpChunk = 128;
}

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

#ifdef ICARION_USE_GPU
void SimulationEngine::initialize_gpu(bool enable_gpu) {
    PROFILE_SCOPE_IF_ENABLED("GPU Initialization");
    
    // Check if GPU is enabled in config
    if (!enable_gpu) {
        output_manager_->log_progress("GPU: Disabled in configuration, using CPU-only");
        return;
    }

    output_manager_->log_progress(
        std::string("GPU: Disabled for v") + ICARION_VERSION +
        " (experimental path remains built but not used). Running CPU-only.");
    return;

    output_manager_->log_progress("GPU: Experimental path (E/B-only). Space-charge, damping, magnetic forces, and multi-domain batches are NOT supported; falling back to CPU in those cases.");
    bool friction_damping_present = false;

    // Validate force setup: only a single ElectricFieldForce per domain, no space charge
    auto gpu_forces_supported = [&]() -> bool {
        if (config_.physics.enable_space_charge) return false;
        for (const auto& reg : force_registries_) {
            if (!reg) continue;
            if (reg->space_charge_model()) return false;
            const auto& forces = reg->forces();
            for (const auto& f : forces) {
                if (dynamic_cast<physics::DampingForce*>(f.get())) {
                    const auto name = f->name();
                    if (name.find("Friction") != std::string::npos) {
                        friction_damping_present = true;
                    }
                }
            }
            if (forces.size() != 1) return false;
            if (forces.front()->name().find("ElectricField") != 0) return false;
        }
        return true;
    };

    if (!gpu_forces_supported()) {
        output_manager_->log_progress("GPU: Disabled — unsupported force combination (requires exactly one ElectricFieldForce per domain, no space charge/damping/magnetic).");
        return;
    }
    
    try {
        // Check if CUDA is available
        if (!icarion::gpu::GPUContext::is_cuda_available()) {
            output_manager_->log_progress("GPU: CUDA not available, using CPU-only");
            return;
        }
        
        // Create GPU context (device 0)
        gpu_context_ = icarion::gpu::GPUContext::create(0);
        if (!gpu_context_) {
            output_manager_->log_progress("GPU: Failed to create context, using CPU-only");
            return;
        }
        
        // Log GPU properties
        const auto& props = gpu_context_->get_properties();
        std::ostringstream gpu_msg;
        gpu_msg << "GPU: " << props.name 
                << " (Compute " << props.compute_capability_major << "." << props.compute_capability_minor
                << ", " << props.total_memory / (1024*1024) << " MB)";
        output_manager_->log_progress(gpu_msg.str());
        
        // Create GPU collision helper (if collision handler exists)
        if (collision_handler_) {
            // Determine collision model from config
            auto collision_model_enum = config_.physics.collision_model;
            std::string collision_model_str;
            
            // Only create GPU helper for supported models (HSS, EHSS)
            bool gpu_supported = false;
            if (collision_model_enum == config::CollisionModel::HSS) {
                collision_model_str = "HSS";
                gpu_supported = true;
            } else if (collision_model_enum == config::CollisionModel::EHSS) {
                collision_model_str = "EHSS";
                gpu_supported = true;
            }
            
            if (gpu_supported) {
                gpu_collision_helper_ = icarion::gpu::GPUCollisionHelper::create(
                    *gpu_context_, 
                    gpu_collision_threshold_,
                    collision_model_str,
                    config_.simulation.rng_seed
                );
                
                if (gpu_collision_helper_) {
                    std::ostringstream collision_msg;
                    collision_msg << "GPU: Collision processing enabled for N >= " 
                                 << gpu_collision_threshold_ << " ions (" << collision_model_str << ")";
                    output_manager_->log_progress(collision_msg.str());
                    
                    // For EHSS, upload geometry data
                    if (collision_model_enum == config::CollisionModel::EHSS) {
                        // TODO(v1.1): Extract geometry from species database and upload to GPU
                        // This requires access to species_db molecular geometries
                        output_manager_->log_progress("GPU: EHSS geometry upload deferred (using HSS fallback for now)");
                    }
                }
            }
        }
        
        // NOTE: GPU Space Charge (P³M) initialization deferred to first use
        // via lazy initialization in try_gpu_space_charge().
        // This avoids coupling to SpaceChargeConfig which may not exist yet.
        // Full integration with DomainConfig.space_charge pending.

        (void)friction_damping_present; // GPU runtime disabled
    }
    catch (const std::exception& e) {
        output_manager_->log_progress(std::string("GPU: Initialization failed: ") + e.what());
        gpu_space_charge_.reset();
        gpu_collision_helper_.reset();
        gpu_context_.reset();
    }
}

#endif

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
        
        // Log trajectory snapshot (write every write_interval steps)
        if (current_step_ == 0 || current_step_ % config_.simulation.write_interval == 0 ||
            output_manager_->should_write(new_time)) {
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

void SimulationEngine::update_space_charge_models(core::IonEnsemble& ensemble) {
    PROFILE_SCOPE_IF_ENABLED("Space Charge Update");
    std::vector<const physics::ISpaceChargeModel*> updated;
    updated.reserve(force_registries_.size());

    for (const auto& registry : force_registries_) {
        if (!registry) continue;
        auto* model = registry->space_charge_model();
        if (!model) continue;

        if (std::find(updated.begin(), updated.end(), model) != updated.end()) {
            continue;  // already updated (shared model)
        }
        model->update_fields(ensemble, current_time_);
        updated.push_back(model);
    }
}

void SimulationEngine::update_dynamic_environments(double t) {
    for (auto& domain : config_.domains) {
        if (!domain.environment.has_dynamic_pressure()) {
            continue;
        }

        try {
            domain.environment.update_time_dependent(t, config_.waveforms);
        } catch (const std::exception& e) {
            throw std::runtime_error(
                "Failed to evaluate pressure waveform for domain '" + domain.name + "': " + e.what());
        }
    }
}

double SimulationEngine::process_timestep(core::IonEnsemble& ensemble) {
    const size_t n_ions = ensemble.size();
    if (dt_per_ion_.size() != n_ions) {
        dt_per_ion_.assign(n_ions, config_.simulation.dt_s);
    }

    update_dynamic_environments(current_time_);

    // Refresh environment cache from SSOT domains (prevents drift)
    {
        const auto* temp = ensemble.temperature_data();
        (void)temp;  // suppress unused in release
        for (size_t i = 0; i < n_ions; ++i) {
            int dom = ensemble.domain_index(i);
            if (dom < 0 || static_cast<size_t>(dom) >= config_.domains.size()) continue;
            const auto& env = config_.domains[static_cast<size_t>(dom)].environment;
            ensemble.temperature_data()[i] = env.temperature_K;
            ensemble.gas_density_data()[i] = env.particle_density_m_3;
            ensemble.neutral_mass_data()[i] = env.gas_mass_kg;
        }
    }

    bool has_space_charge = false;
    for (const auto& reg : force_registries_) {
        if (reg && reg->space_charge_model()) {
            has_space_charge = true;
            break;
        }
    }
    if (has_space_charge && integrator_->is_adaptive() && !adaptive_space_charge_enabled()) {
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

    // Capture pre-integration positions for boundary intersection
    std::vector<Vec3> pos_before;
    pos_before.reserve(n_ions);
    const auto* pos_x_before = ensemble.pos_x_data();
    const auto* pos_y_before = ensemble.pos_y_data();
    const auto* pos_z_before = ensemble.pos_z_data();
    for (size_t i = 0; i < n_ions; ++i) {
        pos_before.emplace_back(pos_x_before[i], pos_y_before[i], pos_z_before[i]);
    }

    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* active = ensemble.active_data();
    auto* born = ensemble.born_data();

    std::vector<int> integration_domains(n_ions, -1);

    const bool use_omp = parallel_enabled_;

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int i = 0; i < static_cast<int>(n_ions); ++i) {
            if (!born[i] && current_time_ >= ensemble.birth_time(i)) {
                born[i] = 1;
                active[i] = 1;
            }

            if (!active[i] || !born[i]) {
                continue;
            }

            Vec3 pos(pos_x[i], pos_y[i], pos_z[i]);
            int domain_idx = -1;
            {
                PROFILE_SCOPE_IF_ENABLED("Domain Finding");
                domain_idx = domain_manager_->find_domain_index(pos);
            }
            if (domain_idx < 0) {
                active[i] = 0;
                ensemble.set_death_time(i, current_time_);
                integration_domains[i] = -1;
                continue;
            }

            integration_domains[i] = domain_idx;
            const auto& domain_config = config_.domains[domain_idx];

            if (ensemble.domain_index(i) != domain_idx) {
                ensemble.update_domain_cache(i, domain_idx,
                    domain_config.environment.temperature_K,
                    domain_config.environment.particle_density_m_3,
                    domain_config.environment.gas_mass_kg);
            }
        }
    }

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

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int i = 0; i < static_cast<int>(n_ions); ++i) {
            if (!active[i] || !born[i]) {
                continue;
            }

            Vec3 pos_after(pos_x[i], pos_y[i], pos_z[i]);
            int new_domain_idx;
            {
                PROFILE_SCOPE_IF_ENABLED("Boundary Checks");
                new_domain_idx = domain_manager_->find_domain_index(pos_after);
            }
            if (new_domain_idx < 0) {
                const int prev_dom = ensemble.domain_index(i);
                new_domain_idx = domain_manager_->forward_axial_bridge_domain(prev_dom, pos_after);
                if (new_domain_idx < 0) {
                    active[i] = 0;
                    ensemble.set_death_time(i, current_time_);
                    continue;
                }
            } else {
                // Still inside current domain; allow seamless hand-off at a shared boundary
                const int cur_dom_idx = ensemble.domain_index(i);
                const int handoff_domain = domain_manager_->shared_boundary_handoff_domain(
                    cur_dom_idx, new_domain_idx, pos_after);
                if (handoff_domain >= 0) {
                    // Hand off to next domain without deactivation
                    new_domain_idx = handoff_domain;
                }
            }
            if (new_domain_idx != ensemble.domain_index(i)) {
                const auto& new_dom = config_.domains[new_domain_idx];
                ensemble.update_domain_cache(i, new_domain_idx,
                    new_dom.environment.temperature_K,
                    new_dom.environment.particle_density_m_3,
                    new_dom.environment.gas_mass_kg);
            }

            ensemble.set_time(i, ensemble.time(i) + dt_used_per_ion[i]);
            if (dt_next_per_ion[i] <= 0.0) {
                // Safety: avoid zeros propagating into next step
                dt_next_per_ion[i] = dt_used_per_ion[i] > 0.0
                    ? dt_used_per_ion[i]
                    : config_.simulation.dt_s;
            }

            Vec3 vel_check = ensemble.get_vel(i);
            bool position_valid = ICARION::safety::is_finite(pos_after);
            bool velocity_valid = ICARION::safety::is_finite(vel_check);

            if (!position_valid || !velocity_valid) {
                active[i] = 0;
                ensemble.set_death_time(i, current_time_);

                if (!config_.simulation.enable_safety_logging) {
                    std::cerr << "Warning: Ion " << i << " has invalid state at t = "
                              << ensemble.time(i) << " s, deactivating" << std::endl;
                }
            }
        }
    }

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

double SimulationEngine::perform_integration(core::IonEnsemble& ensemble,
                                           double t,
                                           const std::vector<double>& dt_per_ion,
                                           const std::vector<int>& domain_indices,
                                           std::vector<double>& dt_used_per_ion,
                                           std::vector<double>& dt_next_per_ion) {
    double max_dt_used = 0.0;

    auto has_space_charge = [&]() {
        for (const auto& reg : force_registries_) {
            if (reg && reg->space_charge_model()) return true;
        }
        return false;
    }();

    // Only run batch if dt is uniform across active ions
    bool uniform_dt = false;
    double dt_batch = 0.0;
    if (integrator_->supports_batch()) {
        uniform_dt = true;
        for (size_t i = 0; i < domain_indices.size(); ++i) {
            if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
            if (dt_batch == 0.0) {
                dt_batch = dt_per_ion[i];
            } else if (dt_per_ion[i] != dt_batch) {
                uniform_dt = false;
                break;
            }
        }
        if (uniform_dt && dt_batch > 0.0) {
            PROFILE_SCOPE_IF_ENABLED("Integrator Batch Step");
            if (integrator_->step_batch(ensemble, t, dt_batch, force_registries_, domain_indices)) {
                std::fill(dt_used_per_ion.begin(), dt_used_per_ion.end(), dt_batch);
                std::fill(dt_next_per_ion.begin(), dt_next_per_ion.end(), dt_batch);
                return dt_batch;
            }
        }
    }

    if (has_space_charge) {
        auto refresh_space_charge_stage = [&](double stage_time_s) {
            const double saved_time = current_time_;
            current_time_ = stage_time_s;
            update_space_charge_models(ensemble);
            current_time_ = saved_time;
        };

        if (auto* rk4 = dynamic_cast<RK4Strategy*>(integrator_.get())) {
            if (uniform_dt && dt_batch > 0.0) {
                PROFILE_SCOPE_IF_ENABLED("Integrator RK4 SC Batch");
                if (rk4->step_batch_with_stage_refresh(
                        ensemble, t, dt_batch, force_registries_, domain_indices,
                        refresh_space_charge_stage)) {
                    std::fill(dt_used_per_ion.begin(), dt_used_per_ion.end(), dt_batch);
                    std::fill(dt_next_per_ion.begin(), dt_next_per_ion.end(), dt_batch);
                    dt_per_ion_.assign(ensemble.size(), dt_batch);
                    return dt_batch;
                }
            }
            output_manager_->log_progress(
                "Warning: Space-charge present without uniform RK4 dt; falling back to per-ion integration with stale fields.");
        } else if (auto* rk45 = dynamic_cast<RK45Strategy*>(integrator_.get())) {
            double dt_seed = dt_batch;
            if (dt_seed <= 0.0) {
                for (size_t i = 0; i < ensemble.size(); ++i) {
                    if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                    dt_seed = dt_per_ion[i];
                    break;
                }
            }
            if (dt_seed <= 0.0) {
                dt_seed = config_.simulation.dt_s;
            }

            if (integrator_->is_adaptive()) {
                PROFILE_SCOPE_IF_ENABLED("Integrator RK45 SC Adaptive");
                const auto result = rk45->step_batch_adaptive(
                    ensemble, t, dt_seed, force_registries_, domain_indices,
                    refresh_space_charge_stage);
                if (result.accepted) {
                    for (size_t i = 0; i < ensemble.size(); ++i) {
                        if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                        dt_used_per_ion[i] = result.dt_used;
                        dt_next_per_ion[i] = result.dt_next;
                    }
                    dt_per_ion_.assign(ensemble.size(), result.dt_next);
                    return result.dt_used;
                }
            } else if (uniform_dt && dt_batch > 0.0) {
                PROFILE_SCOPE_IF_ENABLED("Integrator RK45 SC Batch");
                if (rk45->step_batch_fixed(
                        ensemble, t, dt_batch, force_registries_, domain_indices,
                        refresh_space_charge_stage)) {
                    std::fill(dt_used_per_ion.begin(), dt_used_per_ion.end(), dt_batch);
                    std::fill(dt_next_per_ion.begin(), dt_next_per_ion.end(), dt_batch);
                    dt_per_ion_.assign(ensemble.size(), dt_batch);
                    return dt_batch;
                }
            }
        }
    }

    auto* rk45 = integrator_->is_adaptive()
        ? dynamic_cast<RK45Strategy*>(integrator_.get())
        : nullptr;

    const size_t n = ensemble.size();
    for (size_t i = 0; i < n; ++i) {
        if (domain_indices[i] < 0 || !ensemble.is_active(i)) {
            continue;
        }
        const int dom = domain_indices[i];
        if (dom < 0 || dom >= static_cast<int>(force_registries_.size())) {
            continue;
        }
        const auto& registry = force_registries_[dom];
        if (!registry) {
            continue;
        }
        const double dt = dt_per_ion[i];
        PROFILE_SCOPE_IF_ENABLED("Integrator Step");
        integrator_->step(ensemble, i, t, dt, *registry);

        if (rk45) {
            double used = rk45->last_dt_used();
            double next = rk45->last_dt_suggested();
            if (used > 0.0) {
                dt_used_per_ion[i] = used;
                max_dt_used = std::max(max_dt_used, used);
            }
            if (next > 0.0) {
                dt_next_per_ion[i] = next;
            }
        } else {
            dt_used_per_ion[i] = dt;
            dt_next_per_ion[i] = dt;
            max_dt_used = std::max(max_dt_used, dt);
        }
    }
    return max_dt_used;
}

void SimulationEngine::perform_collisions(core::IonEnsemble& ensemble,
                                          const std::vector<double>& dt_used_per_ion,
                                          const std::vector<int>& domain_indices) {
    if (!collision_handler_ || domain_indices.empty()) {
        return;
    }

    PROFILE_SCOPE_IF_ENABLED("Collision Handling");

    const bool has_batch = collision_handler_->supports_batch();
    const bool multi_event_mode = config_.physics.collision_multi_event_mode;
    const bool split_centered = config_.physics.collision_time_centered;
    const bool split_randomized = config_.physics.collision_time_randomized;
    const int configured_subcycles = std::max(1, config_.physics.collision_subcycles_per_step);
    const int max_events = std::max(1, config_.physics.collision_max_events_per_step);
    const int subcycles = multi_event_mode
        ? std::max(configured_subcycles, max_events)
        : configured_subcycles;
    const bool needs_cpu_split = split_centered || split_randomized || subcycles > 1;
    const bool deep_collision_enabled = deep_collision_diagnostics_.enabled();
    const size_t domain_count = config_.domains.size();
    auto per_domain = group_active_indices_by_domain(ensemble, domain_indices, domain_count);

    for (size_t dom = 0; dom < domain_count; ++dom) {
        auto& indices = per_domain[dom];
        if (indices.empty()) {
            continue;
        }
        const auto& env = config_.domains[dom].environment;
        bool handled = false;
        if (has_batch && !multi_event_mode && !needs_cpu_split && !deep_collision_enabled) {
            // Use batch path only if all dt are equal (avoids bias)
            double dt_batch = dt_used_per_ion[indices.front()];
            bool uniform_dt = true;
            for (size_t idx : indices) {
                if (dt_used_per_ion[idx] != dt_batch) {
                    uniform_dt = false;
                    break;
                }
            }
            if (uniform_dt) {
                handled = collision_handler_->handle_batch(
                    ensemble, indices, dt_batch, env, rng_by_ion_);
                if (handled) {
                    collision_runtime_stats_.add_incomplete_batch_attempts(indices.size());
                }
            }
        }
        if (!handled) {
            handle_collisions_cpu(ensemble, dt_used_per_ion, indices, env, static_cast<int>(dom));
        }
    }
}

void SimulationEngine::handle_collisions_cpu(core::IonEnsemble& ensemble,
                                             const std::vector<double>& dt_used_per_ion,
                                             const std::vector<size_t>& indices,
                                             const config::EnvironmentConfig& env,
                                             int domain_index) {
    if (!collision_handler_ || indices.empty()) {
        return;
    }

    const auto* active = ensemble.active_data();
    const auto* born = ensemble.born_data();
    const bool use_omp = parallel_enabled_;
    const bool multi_event_mode = config_.physics.collision_multi_event_mode;
    const bool split_centered = config_.physics.collision_time_centered;
    const bool split_randomized = config_.physics.collision_time_randomized;
    const int configured_subcycles = std::max(1, config_.physics.collision_subcycles_per_step);
    const int max_events = std::max(1, config_.physics.collision_max_events_per_step);
    const int collision_substeps = multi_event_mode
        ? std::max(configured_subcycles, max_events)
        : configured_subcycles;
    uint64_t macro_attempts_local = 0;
    uint64_t substep_attempts_local = 0;
    uint64_t events_local = 0;
    std::atomic<bool> collision_exception{false};
    std::string collision_exception_msg;

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk) reduction(+:macro_attempts_local,substep_attempts_local,events_local)
        for (int k = 0; k < static_cast<int>(indices.size()); ++k) {
            if (collision_exception.load(std::memory_order_relaxed)) {
                continue;
            }
            size_t ion_idx = indices[static_cast<size_t>(k)];
            if (!active[ion_idx] || !born[ion_idx]) {
                continue;
            }
            macro_attempts_local += 1;
            const double dt_total = dt_used_per_ion[ion_idx];
            std::array<double, 2> primary_segments{dt_total, 0.0};
            int primary_count = 1;
            if (split_randomized) {
                const double alpha = rng_by_ion_[ion_idx].uniform01();
                primary_segments[0] = dt_total * alpha;
                primary_segments[1] = dt_total * (1.0 - alpha);
                primary_count = 2;
            } else if (split_centered) {
                primary_segments[0] = 0.5 * dt_total;
                primary_segments[1] = 0.5 * dt_total;
                primary_count = 2;
            }

            double elapsed = 0.0;
            for (int seg = 0; seg < primary_count; ++seg) {
                const double seg_dt = primary_segments[seg];
                if (seg_dt <= 0.0) {
                    continue;
                }
                const double collision_dt = seg_dt / static_cast<double>(collision_substeps);
                for (int substep = 0; substep < collision_substeps; ++substep) {
                    if (collision_exception.load(std::memory_order_relaxed)) {
                        break;
                    }
                    substep_attempts_local += 1;
                    auto view = ensemble.collision_data(ion_idx);
                    physics::CollisionEventDiagnostics event_diag{};
                    const bool collect_deep = deep_collision_diagnostics_.enabled();
                    const Vec3 pos_before = collect_deep ? view.kin.pos() : Vec3{};
                    const Vec3 vel_before = collect_deep ? view.kin.vel() : Vec3{};
                    bool collided = false;
                    try {
                        collided = collision_handler_->handle_collision(
                            view,
                            collision_dt,
                            rng_by_ion_[ion_idx],
                            env,
                            collect_deep ? &event_diag : nullptr);
                    } catch (const std::exception& e) {
                        if (!collision_exception.exchange(true)) {
                            std::ostringstream oss;
                            oss << "Collision handling failed in domain '"
                                << config_.domains[static_cast<size_t>(domain_index)].name
                                << "' for ion index " << ion_idx
                                << " (species='" << ensemble.species_id(ion_idx) << "'): "
                                << e.what();
                            #pragma omp critical(simengine_collision_exception)
                            {
                                collision_exception_msg = oss.str();
                            }
                        }
                        break;
                    } catch (...) {
                        if (!collision_exception.exchange(true)) {
                            std::ostringstream oss;
                            oss << "Collision handling failed in domain '"
                                << config_.domains[static_cast<size_t>(domain_index)].name
                                << "' for ion index " << ion_idx
                                << " (species='" << ensemble.species_id(ion_idx) << "'): unknown exception";
                            #pragma omp critical(simengine_collision_exception)
                            {
                                collision_exception_msg = oss.str();
                            }
                        }
                        break;
                    }
                    if (collided) {
                        events_local += 1;
                    }
                    if (collect_deep) {
                        const double event_time = current_time_ + elapsed + 0.5 * collision_dt;
                        deep_collision_diagnostics_.note_collision(
                            ion_idx,
                            domain_index,
                            collided,
                            event_time,
                            view.kin.get_mass(),
                            pos_before,
                            event_diag.v_rel_before_m_s,
                            event_diag.sigma_mt_m2,
                            vel_before,
                            view.kin.vel(),
                            !active[ion_idx]);
                    }
                    elapsed += collision_dt;
                }
            }
        }
    }

    collision_runtime_stats_.add_cpu_counts(macro_attempts_local, substep_attempts_local, events_local);
    if (collision_exception.load()) {
        if (collision_exception_msg.empty()) {
            throw std::runtime_error("Collision handling failed due to an exception in parallel section");
        }
        throw std::runtime_error(collision_exception_msg);
    }
}

void SimulationEngine::perform_reactions(core::IonEnsemble& ensemble,
                                         const std::vector<double>& dt_used_per_ion,
                                         const std::vector<int>& domain_indices) {
    if (!reaction_handler_ || config_.reaction_db.reactions.empty()) {
        return;
    }

    PROFILE_SCOPE_IF_ENABLED("Reaction Handling");

    const auto active_indices = collect_active_domain_indices(ensemble, domain_indices);

    if (reaction_handler_->supports_batch()) {
        // Batch only if dt is uniform across active ions
        bool uniform_dt = !active_indices.empty();
        double dt_batch = uniform_dt ? dt_used_per_ion[active_indices.front()] : 0.0;
        for (size_t idx : active_indices) {
            if (dt_used_per_ion[idx] != dt_batch) {
                uniform_dt = false;
                break;
            }
        }
        if (uniform_dt) {
            const bool handled = reaction_handler_->handle_batch(
                ensemble,
                domain_indices,
                dt_batch,
                config_.reaction_db,
                config_.species_db,
                config_.domains,
                rng_by_ion_
            );
            if (handled) {
                return;
            }
        }
    }

    const bool use_omp = parallel_enabled_;

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
        for (int k = 0; k < static_cast<int>(active_indices.size()); ++k) {
            const size_t i = active_indices[static_cast<size_t>(k)];
            const int dom = domain_indices[i];
            auto reaction_view = ensemble.reaction_data(i);
            reaction_handler_->handle_reaction(
                reaction_view,
                dt_used_per_ion[i],
                rng_by_ion_[i],
                config_.reaction_db,
                config_.species_db,
                config_.domains[static_cast<size_t>(dom)].environment);
        }
    }
}

}  // namespace integrator
}  // namespace ICARION

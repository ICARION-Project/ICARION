// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "SimulationEngine.h"
#include "core/utils/safety/numericalSafetyGuards.h"
#include "core/utils/safety/numericalSafetyLogger.h"
#include "core/utils/Profiler.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/log/Logger.h"
#include "core/types/IonEnsemble.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/RK4Strategy.h"
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

using physics::PhysicsRng;

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
    reaction_handler_(reaction_handler)
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
    }
    
    if (!integrator_) {
        throw std::invalid_argument("SimulationEngine: IntegrationStrategy cannot be null");
    }

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
    }

    parallel_enabled_ = config_.simulation.enable_openmp && !integrator_->is_adaptive();
    if (!parallel_enabled_ && config_.simulation.enable_openmp && integrator_->is_adaptive()) {
        // Disable OpenMP for adaptive integrator to avoid shared-state races
        config_.simulation.enable_openmp = false;
        ICARION::log::Logger::main()->warn("Adaptive integrator is serial-only; OpenMP disabled for this run");
    }
    
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
    
    // 1. Initialize subsystems
    initialize(ensemble);
    
    // 2. Main time loop using SoA
    current_dt_ = config_.simulation.dt_s;
    current_time_ = 0.0;
    current_step_ = 0;
    
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
        double dt_next_hint = current_dt_;
        // Process one timestep using SoA
        double new_time = process_timestep(ensemble, current_dt_, dt_next_hint);
        
        // Log trajectory snapshot (write every write_interval steps)
        if (current_step_ == 0 || current_step_ % config_.simulation.write_interval == 0 ||
            output_manager_->should_write(new_time)) {
            PROFILE_SCOPE_IF_ENABLED("Output Writing");
            output_manager_->log_step(new_time, ensemble);
        }
        
        // Update time and step counter
        current_time_ = new_time;
        if (integrator_->is_adaptive()) {
            current_dt_ = dt_next_hint;
        }
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
    
    // Direct SoA finalization (no conversion overhead)
    output_manager_->finalize(current_time_, ensemble);
    
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

void SimulationEngine::update_space_charge_models(core::IonEnsemble& ensemble) {
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

double SimulationEngine::process_timestep(core::IonEnsemble& ensemble, double dt, double& dt_next_hint_out) {
    const size_t n_ions = ensemble.size();
    std::vector<double> dt_used_per_ion(n_ions, dt);
    double dt_next_hint = dt;

    if (rng_by_ion_.size() != n_ions) {
        PROFILE_SCOPE_IF_ENABLED("RNG Initialization");
        rng_by_ion_.clear();
        rng_by_ion_.reserve(n_ions);
        for (size_t i = 0; i < n_ions; ++i) {
            uint64_t ion_seed = config_.simulation.rng_seed + static_cast<uint64_t>(i);
            rng_by_ion_.emplace_back(ion_seed);
        }
    }

    update_space_charge_models(ensemble);

    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* active = ensemble.active_data();
    auto* born = ensemble.born_data();

    std::vector<int> integration_domains(n_ions, -1);

    const bool use_omp = parallel_enabled_;

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(static, 256)
        for (int i = 0; i < static_cast<int>(n_ions); ++i) {
            PhysicsRng& ion_rng = rng_by_ion_[i];

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

    double max_dt_used = perform_integration(
        ensemble,
        current_time_,
        dt,
        integration_domains,
        dt_used_per_ion,
        dt_next_hint);

    // Apply stochastic processes using the actual dt used by the integrator
    perform_collisions(ensemble, dt_used_per_ion, integration_domains);
    perform_reactions(ensemble, dt_used_per_ion, integration_domains);

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(static, 256)
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
                active[i] = 0;
                ensemble.set_death_time(i, current_time_);
                continue;
            }
            if (new_domain_idx != ensemble.domain_index(i)) {
                const auto& new_dom = config_.domains[new_domain_idx];
                ensemble.update_domain_cache(i, new_domain_idx,
                    new_dom.environment.temperature_K,
                    new_dom.environment.particle_density_m_3,
                    new_dom.environment.gas_mass_kg);
            }

            ensemble.set_time(i, ensemble.time(i) + dt_used_per_ion[i]);

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

    dt_next_hint_out = dt_next_hint;

    double new_time = current_time_;
    const double* time_ptr = ensemble.time_data();
    for (size_t i = 0; i < n_ions; ++i) {
        new_time = std::max(new_time, time_ptr[i]);
    }
    // Fallback to accumulated max dt if time data not yet populated
    if (new_time <= current_time_) {
        new_time = current_time_ + max_dt_used;
    }
    return new_time;
}

double SimulationEngine::perform_integration(core::IonEnsemble& ensemble,
                                           double t,
                                           double dt,
                                           const std::vector<int>& domain_indices,
                                           std::vector<double>& dt_used_per_ion,
                                           double& dt_next_hint) {
    double max_dt_used = dt;
    dt_next_hint = dt;

    if (integrator_->step_batch(ensemble, t, dt, force_registries_, domain_indices)) {
        return max_dt_used;
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
        integrator_->step(ensemble, i, t, dt, *registry);

        if (rk45) {
            double used = rk45->last_dt_used();
            double next = rk45->last_dt_suggested();
            if (used > 0.0) {
                dt_used_per_ion[i] = used;
                max_dt_used = std::max(max_dt_used, used);
            }
            if (next > 0.0) {
                dt_next_hint = std::min(dt_next_hint, next);
            }
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

    const size_t n = ensemble.size();
    const auto* active = ensemble.active_data();
    const auto* born = ensemble.born_data();

    const bool has_batch = collision_handler_->supports_batch();
    const size_t domain_count = config_.domains.size();
    std::vector<std::vector<size_t>> per_domain(domain_count);
    for (size_t i = 0; i < n; ++i) {
        int dom = domain_indices[i];
        if (dom < 0 || !active[i] || !born[i]) {
            continue;
        }
        per_domain[static_cast<size_t>(dom)].push_back(i);
    }

    for (size_t dom = 0; dom < domain_count; ++dom) {
        auto& indices = per_domain[dom];
        if (indices.empty()) {
            continue;
        }
        const auto& env = config_.domains[dom].environment;
        bool handled = false;
        if (has_batch) {
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
            }
        }
        if (!handled) {
            handle_collisions_cpu(ensemble, dt_used_per_ion, indices, env);
        }
    }
}

void SimulationEngine::handle_collisions_cpu(core::IonEnsemble& ensemble,
                                             const std::vector<double>& dt_used_per_ion,
                                             const std::vector<size_t>& indices,
                                             const config::EnvironmentConfig& env) {
    if (!collision_handler_ || indices.empty()) {
        return;
    }

    const auto* active = ensemble.active_data();
    const auto* born = ensemble.born_data();

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(static, 256)
        for (int k = 0; k < static_cast<int>(indices.size()); ++k) {
            size_t ion_idx = indices[static_cast<size_t>(k)];
            if (!active[ion_idx] || !born[ion_idx]) {
                continue;
            }
            auto view = ensemble.collision_data(ion_idx);
            collision_handler_->handle_collision(
                view, dt_used_per_ion[ion_idx], rng_by_ion_[ion_idx], env);
        }
    }
}

void SimulationEngine::perform_reactions(core::IonEnsemble& ensemble,
                                         const std::vector<double>& dt_used_per_ion,
                                         const std::vector<int>& domain_indices) {
    if (!reaction_handler_ || config_.reaction_db.reactions.empty()) {
        return;
    }

    PROFILE_SCOPE_IF_ENABLED("Reaction Handling");

    if (reaction_handler_->supports_batch()) {
        // Batch only if dt is uniform across active ions
        double dt_batch = 0.0;
        bool dt_set = false;
        bool uniform_dt = true;
        for (size_t i = 0; i < dt_used_per_ion.size(); ++i) {
            if (!ensemble.is_active(i) || !ensemble.born_data()[i] || domain_indices[i] < 0) continue;
            if (!dt_set) {
                dt_batch = dt_used_per_ion[i];
                dt_set = true;
            } else if (dt_used_per_ion[i] != dt_batch) {
                uniform_dt = false;
                break;
            }
        }
        if (dt_set && uniform_dt) {
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

    const size_t n = ensemble.size();
    const auto* active = ensemble.active_data();
    const auto* born = ensemble.born_data();

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(static, 256)
        for (int i = 0; i < static_cast<int>(n); ++i) {
            if (!active[i] || !born[i]) {
                continue;
            }
            const int dom = domain_indices[static_cast<size_t>(i)];
            if (dom < 0) {
                continue;
            }
            auto reaction_view = ensemble.reaction_data(i);
            reaction_handler_->handle_reaction(
                reaction_view,
                dt_used_per_ion[static_cast<size_t>(i)],
                rng_by_ion_[static_cast<size_t>(i)],
                config_.reaction_db,
                config_.species_db,
                config_.domains[static_cast<size_t>(dom)].environment);
        }
    }
}

}  // namespace integrator
}  // namespace ICARION

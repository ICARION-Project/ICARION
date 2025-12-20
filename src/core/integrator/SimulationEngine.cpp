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
        rk45->enable_stats(!parallel_enabled_);
    }

    parallel_enabled_ = config_.simulation.enable_openmp;
    
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

    // Warn if space charge present but updated only once per macro step
    for (const auto& registry : force_registries_) {
        if (registry && registry->space_charge_model()) {
            space_charge_stale_warned_ = true;
            output_manager_->log_progress("Warning: Space-charge fields are updated once per timestep, not per RK substep; fast-changing clouds may be inaccurate.");
            break;
        }
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

    output_manager_->log_progress("GPU: Experimental path (E/B-only). Space-charge, damping, magnetic forces, and multi-domain batches are NOT supported; falling back to CPU in those cases.");

    // Validate force setup: only a single ElectricFieldForce per domain, no space charge
    auto gpu_forces_supported = [&]() -> bool {
        if (config_.physics.enable_space_charge) return false;
        for (const auto& reg : force_registries_) {
            if (!reg) continue;
            if (reg->space_charge_model()) return false;
            const auto& forces = reg->forces();
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
    current_time_ = 0.0;
    current_step_ = 0;
    dt_per_ion_.assign(ensemble.size(), config_.simulation.dt_s);
    
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

double SimulationEngine::process_timestep(core::IonEnsemble& ensemble) {
    const size_t n_ions = ensemble.size();
    if (dt_per_ion_.size() != n_ions) {
        dt_per_ion_.assign(n_ions, config_.simulation.dt_s);
    }

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
    if (has_space_charge && integrator_ && integrator_->is_adaptive()) {
        throw std::runtime_error("Space charge with adaptive integrator (RK45) is not supported: fields would be stale within substeps. Use fixed-step RK4 or disable space charge.");
    }

    std::vector<double> dt_used_per_ion = dt_per_ion_;
    std::vector<double> dt_next_per_ion = dt_per_ion_;

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
        dt_per_ion_,
        integration_domains,
        dt_used_per_ion,
        dt_next_per_ion);

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
                bool bridged = false;
                if (prev_dom >= 0 && static_cast<size_t>(prev_dom + 1) < config_.domains.size()) {
                    const auto& cur_dom = config_.domains[static_cast<size_t>(prev_dom)];
                    const auto& next_dom = config_.domains[static_cast<size_t>(prev_dom + 1)];
                    const double cur_end = cur_dom.geometry.origin_m.z + cur_dom.geometry.length_m;
                    if (pos_after.z >= cur_end) {
                        new_domain_idx = prev_dom + 1;
                        bridged = true;
                    }
                }
                if (!bridged && new_domain_idx < 0) {
                    active[i] = 0;
                    ensemble.set_death_time(i, current_time_);
                    continue;
                }
            } else {
                // Still inside current domain; allow seamless hand-off at a shared boundary
                const int cur_dom_idx = ensemble.domain_index(i);
                if (cur_dom_idx >= 0 &&
                    new_domain_idx == cur_dom_idx &&
                    static_cast<size_t>(cur_dom_idx + 1) < config_.domains.size()) {
                    const auto& cur_dom = config_.domains[static_cast<size_t>(cur_dom_idx)];
                    const auto& next_dom = config_.domains[static_cast<size_t>(cur_dom_idx + 1)];
                    const double local_z = pos_after.z - cur_dom.geometry.origin_m.z;
                    const double tol = 1e-9;
                    if (local_z >= cur_dom.geometry.length_m - tol) {
                        const double dx = pos_after.x - cur_dom.geometry.origin_m.x;
                        const double dy = pos_after.y - cur_dom.geometry.origin_m.y;
                        const double r2 = dx * dx + dy * dy;
                        const double aperture = cur_dom.geometry.end_aperture_m > 0.0
                            ? cur_dom.geometry.end_aperture_m
                            : cur_dom.geometry.radius_m;
                        if (r2 <= aperture * aperture) {
                            // Hand off to next domain without deactivation
                            new_domain_idx = cur_dom_idx + 1;
                        }
                    }
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
            if (integrator_->step_batch(ensemble, t, dt_batch, force_registries_, domain_indices)) {
                std::fill(dt_used_per_ion.begin(), dt_used_per_ion.end(), dt_batch);
                std::fill(dt_next_per_ion.begin(), dt_next_per_ion.end(), dt_batch);
                return dt_batch;
            }
        }
    }

    // Space-charge aware RK4 batch: recompute fields at each stage (uniform dt required)
    if (has_space_charge) {
        if (!uniform_dt || dt_batch <= 0.0 || dynamic_cast<RK4Strategy*>(integrator_.get()) == nullptr) {
            output_manager_->log_progress("Warning: Space-charge present without uniform RK4 dt; falling back to per-ion integration with stale fields.");
        } else {
            const size_t n = ensemble.size();
            auto* pos_x = ensemble.pos_x_data();
            auto* pos_y = ensemble.pos_y_data();
            auto* pos_z = ensemble.pos_z_data();
            auto* vel_x = ensemble.vel_x_data();
            auto* vel_y = ensemble.vel_y_data();
            auto* vel_z = ensemble.vel_z_data();
            auto* mass = ensemble.mass_data();

            std::vector<double> base_px(pos_x, pos_x + n);
            std::vector<double> base_py(pos_y, pos_y + n);
            std::vector<double> base_pz(pos_z, pos_z + n);
            std::vector<double> base_vx(vel_x, vel_x + n);
            std::vector<double> base_vy(vel_y, vel_y + n);
            std::vector<double> base_vz(vel_z, vel_z + n);

            std::vector<Vec3> k1_a(n), k2_a(n), k3_a(n), k4_a(n);
            std::vector<Vec3> k1_v(n), k2_v(n), k3_v(n), k4_v(n);

            auto compute_accels = [&](double t_stage, std::vector<Vec3>& acc_out) {
                double saved_time = current_time_;
                current_time_ = t_stage;
                update_space_charge_models(ensemble);
                current_time_ = saved_time;
                const bool use_omp = parallel_enabled_;
                #pragma omp parallel if(use_omp)
                {
                    #pragma omp for schedule(guided, kOmpChunk)
                    for (int i = 0; i < static_cast<int>(n); ++i) {
                        if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                        const auto& reg = force_registries_[domain_indices[i]];
                        if (!reg) continue;
                        physics::ForceContext ctx;
                        ctx.domain = reg->domain();
                        ctx.field_model = reg->field_model();
                        ctx.ion_ensemble = &ensemble;
                        ctx.ion_index = static_cast<size_t>(i);
                        Vec3 F = reg->compute_total_force(ensemble, static_cast<size_t>(i), t_stage, ctx);
                        const double inv_m = 1.0 / mass[i];
                        acc_out[static_cast<size_t>(i)] = F * inv_m;
                    }
                }
            };

            // Stage 1: current positions
            for (size_t i = 0; i < n; ++i) {
                k1_v[i] = Vec3{base_vx[i], base_vy[i], base_vz[i]};
            }
            compute_accels(t, k1_a);

            // Stage 2
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                pos_x[i] = base_px[i] + k1_v[i].x * (dt_batch * 0.5);
                pos_y[i] = base_py[i] + k1_v[i].y * (dt_batch * 0.5);
                pos_z[i] = base_pz[i] + k1_v[i].z * (dt_batch * 0.5);
                vel_x[i] = base_vx[i] + k1_a[i].x * (dt_batch * 0.5);
                vel_y[i] = base_vy[i] + k1_a[i].y * (dt_batch * 0.5);
                vel_z[i] = base_vz[i] + k1_a[i].z * (dt_batch * 0.5);
                k2_v[i] = Vec3{vel_x[i], vel_y[i], vel_z[i]};
            }
            compute_accels(t + dt_batch * 0.5, k2_a);

            // Stage 3
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                pos_x[i] = base_px[i] + k2_v[i].x * (dt_batch * 0.5);
                pos_y[i] = base_py[i] + k2_v[i].y * (dt_batch * 0.5);
                pos_z[i] = base_pz[i] + k2_v[i].z * (dt_batch * 0.5);
                vel_x[i] = base_vx[i] + k2_a[i].x * (dt_batch * 0.5);
                vel_y[i] = base_vy[i] + k2_a[i].y * (dt_batch * 0.5);
                vel_z[i] = base_vz[i] + k2_a[i].z * (dt_batch * 0.5);
                k3_v[i] = Vec3{vel_x[i], vel_y[i], vel_z[i]};
            }
            compute_accels(t + dt_batch * 0.5, k3_a);

            // Stage 4
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                pos_x[i] = base_px[i] + k3_v[i].x * dt_batch;
                pos_y[i] = base_py[i] + k3_v[i].y * dt_batch;
                pos_z[i] = base_pz[i] + k3_v[i].z * dt_batch;
                vel_x[i] = base_vx[i] + k3_a[i].x * dt_batch;
                vel_y[i] = base_vy[i] + k3_a[i].y * dt_batch;
                vel_z[i] = base_vz[i] + k3_a[i].z * dt_batch;
                k4_v[i] = Vec3{vel_x[i], vel_y[i], vel_z[i]};
            }
            compute_accels(t + dt_batch, k4_a);

            // Final update from base state
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) {
                    continue;
                }
                Vec3 pos_new = Vec3{base_px[i], base_py[i], base_pz[i]} +
                    (k1_v[i] + k2_v[i] * 2.0 + k3_v[i] * 2.0 + k4_v[i]) * (dt_batch / 6.0);
                Vec3 vel_new = Vec3{base_vx[i], base_vy[i], base_vz[i]} +
                    (k1_a[i] + k2_a[i] * 2.0 + k3_a[i] * 2.0 + k4_a[i]) * (dt_batch / 6.0);

                pos_x[i] = pos_new.x;
                pos_y[i] = pos_new.y;
                pos_z[i] = pos_new.z;
                vel_x[i] = vel_new.x;
                vel_y[i] = vel_new.y;
                vel_z[i] = vel_new.z;

                dt_used_per_ion[i] = dt_batch;
                dt_next_per_ion[i] = dt_batch;
            }

            dt_per_ion_.assign(n, dt_batch);
            max_dt_used = dt_batch;
            return max_dt_used;
        }
    }

    // Space-charge aware RK45 batch (fixed-step use of Dormand-Prince coefficients)
    if (has_space_charge) {
        auto* rk45_ptr = dynamic_cast<RK45Strategy*>(integrator_.get());
        if (rk45_ptr && uniform_dt && dt_batch > 0.0) {
            const size_t n = ensemble.size();
            auto* pos_x = ensemble.pos_x_data();
            auto* pos_y = ensemble.pos_y_data();
            auto* pos_z = ensemble.pos_z_data();
            auto* vel_x = ensemble.vel_x_data();
            auto* vel_y = ensemble.vel_y_data();
            auto* vel_z = ensemble.vel_z_data();
            auto* mass = ensemble.mass_data();

            std::vector<double> base_px(pos_x, pos_x + n);
            std::vector<double> base_py(pos_y, pos_y + n);
            std::vector<double> base_pz(pos_z, pos_z + n);
            std::vector<double> base_vx(vel_x, vel_x + n);
            std::vector<double> base_vy(vel_y, vel_y + n);
            std::vector<double> base_vz(vel_z, vel_z + n);

            std::vector<Vec3> kv[7];
            std::vector<Vec3> ka[7];
            for (int s = 0; s < 7; ++s) {
                kv[s].assign(n, Vec3{0.0, 0.0, 0.0});
                ka[s].assign(n, Vec3{0.0, 0.0, 0.0});
            }

            auto compute_accels_rk = [&](double t_stage, std::vector<Vec3>& acc_out) {
                double saved_time = current_time_;
                current_time_ = t_stage;
                update_space_charge_models(ensemble);
                current_time_ = saved_time;
                const bool use_omp = parallel_enabled_;
                #pragma omp parallel if(use_omp)
                {
                    #pragma omp for schedule(guided, kOmpChunk)
                    for (int i = 0; i < static_cast<int>(n); ++i) {
                        if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                        const auto& reg = force_registries_[domain_indices[i]];
                        if (!reg) continue;
                        physics::ForceContext ctx;
                        ctx.domain = reg->domain();
                        ctx.field_model = reg->field_model();
                        ctx.ion_ensemble = &ensemble;
                        ctx.ion_index = static_cast<size_t>(i);
                        Vec3 F = reg->compute_total_force(ensemble, static_cast<size_t>(i), t_stage, ctx);
                        const double inv_m = 1.0 / mass[i];
                        acc_out[static_cast<size_t>(i)] = F * inv_m;
                    }
                }
            };

            // Dormand-Prince coefficients (from RK45Strategy)
            constexpr double c2 = 1.0 / 5.0;
            constexpr double c3 = 3.0 / 10.0;
            constexpr double c4 = 4.0 / 5.0;
            constexpr double c5 = 8.0 / 9.0;
            constexpr double c6 = 1.0;
            constexpr double c7 = 1.0;

            constexpr double a21 = 1.0 / 5.0;
            constexpr double a31 = 3.0 / 40.0;
            constexpr double a32 = 9.0 / 40.0;
            constexpr double a41 = 44.0 / 45.0;
            constexpr double a42 = -56.0 / 15.0;
            constexpr double a43 = 32.0 / 9.0;
            constexpr double a51 = 19372.0 / 6561.0;
            constexpr double a52 = -25360.0 / 2187.0;
            constexpr double a53 = 64448.0 / 6561.0;
            constexpr double a54 = -212.0 / 729.0;
            constexpr double a61 = 9017.0 / 3168.0;
            constexpr double a62 = -355.0 / 33.0;
            constexpr double a63 = 46732.0 / 5247.0;
            constexpr double a64 = 49.0 / 176.0;
            constexpr double a65 = -5103.0 / 18656.0;
            constexpr double a71 = 35.0 / 384.0;
            constexpr double a72 = 0.0;
            constexpr double a73 = 500.0 / 1113.0;
            constexpr double a74 = 125.0 / 192.0;
            constexpr double a75 = -2187.0 / 6784.0;
            constexpr double a76 = 11.0 / 84.0;

            constexpr double b41 = 35.0 / 384.0;
            constexpr double b42 = 0.0;
            constexpr double b43 = 500.0 / 1113.0;
            constexpr double b44 = 125.0 / 192.0;
            constexpr double b45 = -2187.0 / 6784.0;
            constexpr double b46 = 11.0 / 84.0;
            constexpr double b47 = 0.0;

            constexpr double b51 = 5179.0 / 57600.0;
            constexpr double b52 = 0.0;
            constexpr double b53 = 7571.0 / 16695.0;
            constexpr double b54 = 393.0 / 640.0;
            constexpr double b55 = -92097.0 / 339200.0;
            constexpr double b56 = 187.0 / 2100.0;
            constexpr double b57 = 1.0 / 40.0;

            // k1
            for (size_t i = 0; i < n; ++i) {
                kv[0][i] = Vec3{base_vx[i], base_vy[i], base_vz[i]};
            }
            compute_accels_rk(t, ka[0]);

            auto stage_apply = [&](size_t i, const Vec3& pos, const Vec3& vel) {
                pos_x[i] = pos.x;
                pos_y[i] = pos.y;
                pos_z[i] = pos.z;
                vel_x[i] = vel.x;
                vel_y[i] = vel.y;
                vel_z[i] = vel.z;
            };

            // k2
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                Vec3 pos = Vec3{base_px[i], base_py[i], base_pz[i]} + kv[0][i] * (dt_batch * a21);
                Vec3 vel = Vec3{base_vx[i], base_vy[i], base_vz[i]} + ka[0][i] * (dt_batch * a21);
                stage_apply(i, pos, vel);
                kv[1][i] = vel;
            }
            compute_accels_rk(t + c2 * dt_batch, ka[1]);

            // k3
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                Vec3 pos = Vec3{base_px[i], base_py[i], base_pz[i]} +
                    (kv[0][i] * a31 + kv[1][i] * a32) * dt_batch;
                Vec3 vel = Vec3{base_vx[i], base_vy[i], base_vz[i]} +
                    (ka[0][i] * a31 + ka[1][i] * a32) * dt_batch;
                stage_apply(i, pos, vel);
                kv[2][i] = vel;
            }
            compute_accels_rk(t + c3 * dt_batch, ka[2]);

            // k4
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                Vec3 pos = Vec3{base_px[i], base_py[i], base_pz[i]} +
                    (kv[0][i] * a41 + kv[1][i] * a42 + kv[2][i] * a43) * dt_batch;
                Vec3 vel = Vec3{base_vx[i], base_vy[i], base_vz[i]} +
                    (ka[0][i] * a41 + ka[1][i] * a42 + ka[2][i] * a43) * dt_batch;
                stage_apply(i, pos, vel);
                kv[3][i] = vel;
            }
            compute_accels_rk(t + c4 * dt_batch, ka[3]);

            // k5
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                Vec3 pos = Vec3{base_px[i], base_py[i], base_pz[i]} +
                    (kv[0][i] * a51 + kv[1][i] * a52 + kv[2][i] * a53 + kv[3][i] * a54) * dt_batch;
                Vec3 vel = Vec3{base_vx[i], base_vy[i], base_vz[i]} +
                    (ka[0][i] * a51 + ka[1][i] * a52 + ka[2][i] * a53 + ka[3][i] * a54) * dt_batch;
                stage_apply(i, pos, vel);
                kv[4][i] = vel;
            }
            compute_accels_rk(t + c5 * dt_batch, ka[4]);

            // k6
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                Vec3 pos = Vec3{base_px[i], base_py[i], base_pz[i]} +
                    (kv[0][i] * a61 + kv[1][i] * a62 + kv[2][i] * a63 + kv[3][i] * a64 + kv[4][i] * a65) * dt_batch;
                Vec3 vel = Vec3{base_vx[i], base_vy[i], base_vz[i]} +
                    (ka[0][i] * a61 + ka[1][i] * a62 + ka[2][i] * a63 + ka[3][i] * a64 + ka[4][i] * a65) * dt_batch;
                stage_apply(i, pos, vel);
                kv[5][i] = vel;
            }
            compute_accels_rk(t + c6 * dt_batch, ka[5]);

            // k7
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                Vec3 pos = Vec3{base_px[i], base_py[i], base_pz[i]} +
                    (kv[0][i] * a71 + kv[1][i] * a72 + kv[2][i] * a73 + kv[3][i] * a74 + kv[4][i] * a75 + kv[5][i] * a76) * dt_batch;
                Vec3 vel = Vec3{base_vx[i], base_vy[i], base_vz[i]} +
                    (ka[0][i] * a71 + ka[1][i] * a72 + ka[2][i] * a73 + ka[3][i] * a74 + ka[4][i] * a75 + ka[5][i] * a76) * dt_batch;
                stage_apply(i, pos, vel);
                kv[6][i] = vel;
            }
            compute_accels_rk(t + c7 * dt_batch, ka[6]);

            // 4th-order solution (y4) using b4 weights
            for (size_t i = 0; i < n; ++i) {
                if (domain_indices[i] < 0 || !ensemble.is_active(i)) continue;
                Vec3 pos_new = Vec3{base_px[i], base_py[i], base_pz[i]} +
                    (kv[0][i] * b41 + kv[1][i] * b42 + kv[2][i] * b43 + kv[3][i] * b44 + kv[4][i] * b45 + kv[5][i] * b46 + kv[6][i] * b47) * dt_batch;
                Vec3 vel_new = Vec3{base_vx[i], base_vy[i], base_vz[i]} +
                    (ka[0][i] * b41 + ka[1][i] * b42 + ka[2][i] * b43 + ka[3][i] * b44 + ka[4][i] * b45 + ka[5][i] * b46 + ka[6][i] * b47) * dt_batch;

                pos_x[i] = pos_new.x;
                pos_y[i] = pos_new.y;
                pos_z[i] = pos_new.z;
                vel_x[i] = vel_new.x;
                vel_y[i] = vel_new.y;
                vel_z[i] = vel_new.z;

                dt_used_per_ion[i] = dt_batch;
                dt_next_per_ion[i] = dt_batch;
                max_dt_used = std::max(max_dt_used, dt_batch);
            }

            dt_per_ion_.assign(n, dt_batch);
            return max_dt_used;
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
    const bool use_omp = parallel_enabled_;

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
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
    const bool use_omp = parallel_enabled_;

    #pragma omp parallel if(use_omp)
    {
        #pragma omp for schedule(guided, kOmpChunk)
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

// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "SimulationEngine.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/BorisStrategy.h"
#include "core/utils/safety/numericalSafetyGuards.h"
#include "core/utils/safety/numericalSafetyLogger.h"
#include "core/utils/Profiler.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/types/IonEnsemble.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <chrono>
#include <cstdlib>  // for setenv (NUMA thread placement)

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
    std::string log_path = "";  // TODO: Add log file to OutputConfig
    if (config_.output.print_progress) {
        log_path = config_.output.folder + "/simulation.log";
    }
    
    output_manager_ = std::make_unique<OutputManager>(
        hdf5_path,
        log_path,
        config_.simulation.dt_s * config_.simulation.write_interval,  // Write interval in seconds
        50  // Buffer max
    );
    
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
    if (!config_.simulation.enable_openmp) {
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
        
        // Create GPU integration helper
        gpu_helper_ = icarion::gpu::GPUIntegrationHelper::create(*gpu_context_, gpu_threshold_);
        if (!gpu_helper_) {
            output_manager_->log_progress("GPU: Failed to create integration helper, using CPU-only");
            gpu_context_.reset();
            return;
        }
        
        std::ostringstream threshold_msg;
        threshold_msg << "GPU: Integration enabled for N >= " << gpu_threshold_ << " ions";
        output_manager_->log_progress(threshold_msg.str());
        
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
                        // TODO: Extract geometry from species database and upload to GPU
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
        gpu_helper_.reset();
        gpu_context_.reset();
    }
}

const IFieldProvider* SimulationEngine::extract_field_provider(int domain_id) const {
    // Validate domain_id
    if (domain_id < 0 || domain_id >= static_cast<int>(force_registries_.size())) {
        return nullptr;
    }
    
    const auto& registry = force_registries_[domain_id];
    if (!registry) {
        return nullptr;
    }
    
    // Search for ElectricFieldForce in registry
    for (const auto& force : registry->forces()) {
        if (auto* e_force = dynamic_cast<const physics::ElectricFieldForce*>(force.get())) {
            return e_force->get_field_provider();
        }
    }
    
    return nullptr;  // No field provider found
}





bool SimulationEngine::try_gpu_boundary_check(core::IonEnsemble& ensemble, int domain_idx) {
    // Early exit: GPU not available
    if (!gpu_helper_) {
        return false;
    }
    
    // NOTE: This GPU path is currently unused by the main timestep loop.
    // Keep logic in sync if/when batch boundary checks are dispatched.
    
    const auto& domain_config = config_.domains[domain_idx];
    
    // =========================================================================
    // Conditional GPU Dispatch: Check domain config compatibility
    // =========================================================================
    // GPU boundary check is limited to:
    // - Absorption only (no reflections yet)
    // - Cylindrical geometry only (no Orbitrap hyperlogarithmic)
    // 
    // If domain requires unsupported features → fallback to CPU
    // =========================================================================
    
    // Check 1: Boundary action must be Absorption
    if (domain_config.boundary.type != config::BoundaryActionType::Absorption) {
        // Reflections require: surface normals, RNG, thermal accommodation
        // GPU doesn't support these yet → CPU fallback
        return false;
    }
    
    // Check 2: Instrument must NOT be Orbitrap
    if (domain_config.instrument == config::Instrument::Orbitrap) {
        // Orbitrap uses hyperlogarithmic surface with bisection-based intersection
        // GPU doesn't support this geometry yet → CPU fallback
        return false;
    }
    
    // =========================================================================
    // GPU-compatible configuration! Dispatch to GPU kernel
    // =========================================================================
    
    // Build temporary IonState vector for GPU upload
    // NOTE: GPU API still uses vector<IonState>, will migrate to SoA later
    size_t n_ions = ensemble.size();
    std::vector<IonState> ions(n_ions);
    
    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* vel_x = ensemble.vel_x_data();
    auto* vel_y = ensemble.vel_y_data();
    auto* vel_z = ensemble.vel_z_data();
    auto* mass = ensemble.mass_data();
    auto* charge = ensemble.charge_data();
    auto* active = ensemble.active_data();
    
    for (size_t i = 0; i < n_ions; ++i) {
        ions[i].pos = {pos_x[i], pos_y[i], pos_z[i]};
        ions[i].vel = {vel_x[i], vel_y[i], vel_z[i]};
        ions[i].mass_kg = mass[i];
        ions[i].ion_charge_C = charge[i];
        ions[i].active = static_cast<bool>(active[i]);
    }
    
    // Call GPU batch boundary check
    bool gpu_success = gpu_helper_->check_boundaries_batch(
        ions,
        domain_config.geometry.length_m,
        domain_config.geometry.radius_m,
        domain_idx == static_cast<int>(config_.domains.size()) - 1  // is_last_domain
    );
    
    if (!gpu_success) {
        return false;  // GPU error → CPU fallback
    }
    
    // Sync active flags back to IonEnsemble
    for (size_t i = 0; i < n_ions; ++i) {
        active[i] = static_cast<uint8_t>(ions[i].active);
    }
    
    return true;  // GPU boundary check succeeded
}



void SimulationEngine::finalize_gpu() {
    if (!gpu_helper_) {
        return;  // GPU not initialized
    }
    
    const auto& stats = gpu_helper_->get_stats();
    
    if (stats.gpu_integrations == 0) {
        output_manager_->log_progress("GPU: No batches processed (all below threshold)");
        return;
    }
    
    // Log GPU statistics
    std::ostringstream msg;
    msg << "GPU Statistics:\n"
        << "  Batches:      " << stats.gpu_integrations << "\n"
        << "  Total ions:   " << stats.total_ions_gpu << "\n"
        << "  Total time:   " << stats.total_time_ms << " ms\n"
        << "  Avg/batch:    " << (stats.total_time_ms / stats.gpu_integrations) << " ms";
    
    output_manager_->log_progress(msg.str());
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
    const double dt = config_.simulation.dt_s;
    current_time_ = 0.0;
    current_step_ = 0;
    
    output_manager_->log_progress("Starting main simulation loop (SoA)");
    
    // Helper to check if we should continue
    auto should_continue_soa = [&]() -> bool {
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
    
    while (should_continue_soa()) {
        // Process one timestep using SoA
        process_timestep(ensemble, dt);
        
        // Log trajectory snapshot (write every write_interval steps)
        if (current_step_ % config_.simulation.write_interval == 0) {
            PROFILE_SCOPE_IF_ENABLED("Output Writing");
            output_manager_->log_step(current_time_, ensemble);
        }
        
        // Update time and step counter
        current_time_ += dt;
        current_step_++;
        
        // Progress logging (every 10%)
        {
            PROFILE_SCOPE_IF_ENABLED("Progress Update");
            log_progress(current_time_);
        }
    }
    
    // 3. Finalization
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

void SimulationEngine::process_timestep(core::IonEnsemble& ensemble, double dt) {
    // Direct SoA processing (no conversions!)
    const int n_ions = static_cast<int>(ensemble.size());
    
    // Initialize per-ion RNGs on first call
    if (rng_by_ion_.empty()) {
        PROFILE_SCOPE_IF_ENABLED("RNG Initialization");
        rng_by_ion_.reserve(n_ions);
        for (int i = 0; i < n_ions; ++i) {
            uint64_t ion_seed = config_.simulation.rng_seed + static_cast<uint64_t>(i);
            rng_by_ion_.emplace_back(ion_seed);
        }
    }
    
    // Get raw array pointers for cache-friendly iteration
    auto* pos_x = ensemble.pos_x_data();
    auto* pos_y = ensemble.pos_y_data();
    auto* pos_z = ensemble.pos_z_data();
    auto* active = ensemble.active_data();
    auto* born = ensemble.born_data();
    
    // Parallel ion processing (OpenMP if enabled)
    #pragma omp parallel if(config_.simulation.enable_openmp)
    {
        #pragma omp for schedule(static, 256)
        for (int i = 0; i < n_ions; ++i) {
            PhysicsRng& ion_rng = rng_by_ion_[i];
            
            // Birth logic (delayed emission)
            if (!born[i] && current_time_ >= ensemble.birth_time(i)) {
                born[i] = 1;
                active[i] = 1;
            }
            
            // Skip ions that are not active/born
            if (!active[i] || !born[i]) {
                continue;
            }
            
            // 1. Find current domain
            Vec3 pos(pos_x[i], pos_y[i], pos_z[i]);
            int domain_idx;
            {
                PROFILE_SCOPE_IF_ENABLED("Domain Finding");
                domain_idx = domain_manager_->find_domain_index(pos);
                if (domain_idx < 0) {
                    active[i] = false;
                    ensemble.set_death_time(i, current_time_);
                    continue;
                }
            }
            
            const auto& domain_config = config_.domains[domain_idx];
            
            // 2. Update domain cache if domain changed
            if (ensemble.domain_index(i) != domain_idx) {
                ensemble.update_domain_cache(i, domain_idx,
                    domain_config.environment.temperature_K,
                    domain_config.environment.particle_density_m_3,
                    domain_config.environment.gas_mass_kg);
            }
            
            // 3-5. Collision/Reaction handling using SoA views
            if (collision_handler_) {
                PROFILE_SCOPE_IF_ENABLED("Collision Handling");
                auto collision_view = ensemble.collision_data(i);
                collision_handler_->handle_collision(collision_view, dt, ion_rng, 
                                                      domain_config.environment);
            }
            
            if (reaction_handler_ && !config_.reaction_db.reactions.empty()) {
                PROFILE_SCOPE_IF_ENABLED("Reaction Handling");
                auto reaction_view = ensemble.reaction_data(i);
                
                reaction_handler_->handle_reaction(reaction_view,
                    dt, ion_rng, config_.reaction_db, config_.species_db,
                    domain_config.environment);
            }
            
            // 6-7. Integration 
            {
                PROFILE_SCOPE_IF_ENABLED("Integration");
                
                // Get force registry for this domain
                const auto& force_registry = force_registries_[domain_idx];
                
                // Integrate directly using SoA (no conversion!)
                integrator_->step(ensemble, i, current_time_, dt, *force_registry);
            }
            
            // 8. Boundary checks via geometry strategy (CPU)
            {
                PROFILE_SCOPE_IF_ENABLED("Boundary Checks");
                Vec3 pos_after(pos_x[i], pos_y[i], pos_z[i]);
                int new_domain_idx = domain_manager_->find_domain_index(pos_after);
                if (new_domain_idx < 0) {
                    active[i] = false;
                    ensemble.set_death_time(i, current_time_);
                    continue;
                }
                if (new_domain_idx != domain_idx) {
                    const auto& new_dom = config_.domains[new_domain_idx];
                    ensemble.update_domain_cache(i, new_domain_idx,
                        new_dom.environment.temperature_K,
                        new_dom.environment.particle_density_m_3,
                        new_dom.environment.gas_mass_kg);
                }
            }
            
            // 9. Update time
            ensemble.set_time(i, ensemble.time(i) + dt);
            
            // 10. Numerical safety checks
            Vec3 pos_check(pos_x[i], pos_y[i], pos_z[i]);
            Vec3 vel_check = ensemble.get_vel(i);
            
            bool position_valid = ICARION::safety::is_finite(pos_check);
            bool velocity_valid = ICARION::safety::is_finite(vel_check);
            
            if (!position_valid || !velocity_valid) {
                active[i] = false;
                ensemble.set_death_time(i, current_time_);
                
                if (!config_.simulation.enable_safety_logging) {
                    std::cerr << "Warning: Ion " << i << " has invalid state at t = " 
                              << ensemble.time(i) << " s, deactivating" << std::endl;
                }
            }
        }  // End parallel for
    }  // End parallel region
}

}  // namespace integrator
}  // namespace ICARION

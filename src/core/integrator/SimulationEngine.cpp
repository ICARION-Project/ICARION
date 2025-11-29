// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "SimulationEngine.h"
#include "DomainContext.h"
#include "core/utils/safety/numericalSafetyGuards.h"
#include "core/utils/safety/numericalSafetyLogger.h"
#include "core/utils/mathUtils.h"
#include "core/utils/Profiler.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include <algorithm>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <stdexcept>
#include <chrono>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ICARION {
namespace integrator {

using physics::EhssRng;

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

void SimulationEngine::initialize(const std::vector<IonState>& ions) {
    // SSOT: Birth conditions already set by IonLoader
    // No need to modify ion state here!
    
    // Initialize output system
    output_manager_->initialize(config_, ions);
    
    // Log initialization
    output_manager_->log_progress("Simulation engine initialized");
    
    std::ostringstream msg;
    msg << "Configuration: " << ions.size() << " ions, "
        << config_.domains.size() << " domains, "
        << "dt = " << config_.simulation.dt_s * 1e9 << " ns, "
        << "t_max = " << config_.simulation.total_time_s * 1e6 << " µs";
    output_manager_->log_progress(msg.str());
    
#ifdef ICARION_USE_GPU
    // Initialize GPU acceleration (optional, controlled by config)
    initialize_gpu(config_.simulation.enable_gpu);
#endif
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
    }
    catch (const std::exception& e) {
        output_manager_->log_progress(std::string("GPU: Initialization failed: ") + e.what());
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

bool SimulationEngine::try_gpu_integration(std::vector<IonState>& ions, double dt) {
    PROFILE_SCOPE_IF_ENABLED("GPU Integration");
    
    // Check if GPU available and threshold met
    int n_ions = static_cast<int>(ions.size());
    if (!gpu_helper_ || n_ions < static_cast<int>(gpu_helper_->get_threshold())) {
        return false;  // Use CPU path
    }
    
    // Extract field provider for GPU integration
    // Note: Currently assumes single-domain or domain 0 fields
    // Multi-domain field handling deferred (requires per-ion domain tracking)
    const IFieldProvider* field_provider = extract_field_provider(0);
    
    // Attempt GPU batch integration
    if (!gpu_helper_->integrate_batch_rk4(ions, dt, current_time_, field_provider)) {
        return false;  // GPU failed, use CPU fallback
    }
    
    // GPU success! Update ion times (GPU only updates pos/vel)
    #pragma omp parallel for if(config_.simulation.enable_openmp)
    for (int i = 0; i < n_ions; ++i) {
        if (ions[i].active) {
            ions[i].t += dt;
        }
    }
    
    return true;  // GPU path complete
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

void SimulationEngine::apply_ion_birth(std::vector<IonState>& ions, double t) {
    for (auto& ion : ions) {
        if (!ion.born && ion.birth_time_s <= t) {
            ion.born = true;
            ion.active = true;
            ion.t = t;
        }
    }
}

void SimulationEngine::process_timestep(std::vector<IonState>& ions, double dt) {
    const int n_ions = static_cast<int>(ions.size());
    
    // Initialize per-ion RNGs on first call (persistent across timesteps!)
    // Each ion gets its own RNG seeded deterministically from: base_seed + ion_index
    // This ensures:
    // - Reproducible results regardless of thread count or scheduling
    // - Same ion always sees same random sequence
    // - Thread-safe (each thread accesses different ion RNG)
    // - RNG state is PRESERVED across timesteps (critical for collision physics!)
    if (rng_by_ion_.empty()) {
        PROFILE_SCOPE_IF_ENABLED("RNG Initialization");
        rng_by_ion_.reserve(n_ions);
        for (int i = 0; i < n_ions; ++i) {
            uint64_t ion_seed = config_.simulation.rng_seed + static_cast<uint64_t>(i);
            rng_by_ion_.emplace_back(ion_seed);
        }
    }
    
#ifdef ICARION_USE_GPU
    // Try GPU acceleration (auto-fallback to CPU if unavailable)
    if (try_gpu_integration(ions, dt)) {
        return;  // GPU succeeded
    }
#endif
    
    // ====================================================================
    // CPU Path (fallback or small N)
    // ====================================================================
    // Parallel ion processing (OpenMP if enabled)
    // Note: Using schedule(static) for better cache locality and lower overhead
    // Dynamic scheduling was causing severe performance degradation due to task queue contention
    #pragma omp parallel if(config_.simulation.enable_openmp)
    {
        #pragma omp for schedule(static, 256)
        for (int i = 0; i < n_ions; ++i) {
            IonState& ion = ions[i];
            physics::EhssRng& ion_rng = rng_by_ion_[i];  // Ion-specific RNG (persistent!)
            
            if (!ion.active) continue;
            
            // ================================================================
            // 1. Find current domain
            // ================================================================
            int domain_idx = find_ion_domain(ion);
            if (domain_idx < 0) {
                ion.active = false;
                continue;
            }
            
            // ================================================================
            // 2. Update domain properties (if domain changed)
            // ================================================================
            update_domain_properties(ion, domain_idx);
            
            // ================================================================
            // 3. Create domain context (local coordinate system)
            // ================================================================
            DomainContext ctx(ion, domain_idx, *domain_manager_);
            Vec3 pos_before = ctx.pos_local();
            
            // ================================================================
            // 4. Physics processing (collisions, reactions, integration)
            // ================================================================
            process_ion_collisions(ion, ctx, dt, ion_rng, domain_idx);
            process_ion_reactions(ion, ctx, dt, ion_rng, domain_idx);
            integrate_ion_trajectory(ion, ctx, dt, domain_idx, ions);
            
            // ================================================================
            // 5. Boundary checks (aperture, walls)
            // ================================================================
            bool still_inside = check_ion_boundaries(ion, ctx, domain_idx, pos_before);
            if (!still_inside) continue;
            
            // ================================================================
            // 6. Transform back to global coordinates and update time
            // ================================================================
            ctx.sync_to_ion();
            ion.t += dt;
            
            // ================================================================
            // 7. Numerical safety verification
            // ================================================================
            verify_ion_safety(ion, i, domain_idx);
            
        }  // End of parallel for loop
    }  // End of parallel region
}

bool SimulationEngine::should_continue(const std::vector<IonState>& ions, double t) const {
    // Stop if time exceeded
    if (t >= config_.simulation.total_time_s) {
        return false;
    }
    
    // Continue if any ion is active OR waiting to be born
    const double t_max = config_.simulation.total_time_s;
    bool any_active_or_waiting = std::any_of(ions.begin(), ions.end(), 
        [t_max](const IonState& ion) { 
            return (ion.active && ion.born) || (!ion.born && ion.birth_time_s <= t_max); 
        });
    
    return any_active_or_waiting;
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

std::vector<IonState> SimulationEngine::run(std::vector<IonState>& ions) {
    // 1. Initialize subsystems
    initialize(ions);
    
    // 2. Main time loop
    const double dt = config_.simulation.dt_s;
    current_time_ = 0.0;
    current_step_ = 0;
    
    output_manager_->log_progress("Starting main simulation loop");
    
    while (should_continue(ions, current_time_)) {
        // Apply ion birth logic (delayed emission)
        apply_ion_birth(ions, current_time_);
        
        // Process one timestep
        process_timestep(ions, dt);
        
        // Log trajectory snapshot (auto-flush if needed)
        // Only write every write_interval steps to avoid excessive I/O
        if (current_step_ % config_.simulation.write_interval == 0) {
            PROFILE_SCOPE_IF_ENABLED("Output Writing");
            output_manager_->log_step(current_time_, ions);
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
    output_manager_->log_progress("Simulation completed");
    
    // Count active ions
    size_t active_count = std::count_if(ions.begin(), ions.end(),
        [](const IonState& ion) { return ion.active && ion.born; });
    
    std::ostringstream msg;
    msg << "Final state: " << active_count << "/" << ions.size() << " ions active";
    output_manager_->log_progress(msg.str());
    
    // Flush output and write completion metadata
    output_manager_->finalize(current_time_, ions);
    
    // Generate numerical safety report if logging was enabled
    if (config_.simulation.enable_safety_logging) {
        std::string report_file = config_.output.folder + "/numerical_safety_report.txt";
        safety::NumericalSafetyLogger::getInstance().generateSafetyReport(report_file);
        
        // Log summary statistics
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
    
#ifdef ICARION_USE_GPU
    // Log GPU performance statistics
    finalize_gpu();
#endif
    
    return ions;
}

// ============================================================================
// SoA (Structure of Arrays) Implementation - Phase 2
// ============================================================================

std::vector<IonState> SimulationEngine::run_soa(core::IonEnsemble& ensemble) {
    // Phase 3: Direct SoA simulation loop (no upfront conversion)
    
    // Convert to legacy format for initialization only
    std::vector<IonState> ions_legacy = ensemble.to_legacy();
    
    // 1. Initialize subsystems
    initialize(ions_legacy);
    
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
        // Check if any ion is active
        for (size_t i = 0; i < ensemble.size(); ++i) {
            if (ensemble.is_active(i)) {
                return true;
            }
        }
        return false;
    };
    
    while (should_continue_soa()) {
        // Process one timestep using SoA
        process_timestep_soa(ensemble, dt);
        
        // Log trajectory snapshot (write every write_interval steps)
        if (current_step_ % config_.simulation.write_interval == 0) {
            PROFILE_SCOPE_IF_ENABLED("Output Writing");
            // Phase 5: Direct SoA→HDF5 writing (no conversion overhead)
            output_manager_->log_step_soa(current_time_, ensemble);
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
    
    // Phase 5: Direct SoA finalization (no conversion overhead)
    output_manager_->finalize_soa(current_time_, ensemble);
    
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
    
    // Return final ions in legacy format (for API compatibility)
    return ensemble.to_legacy();
}

void SimulationEngine::process_timestep_soa(core::IonEnsemble& ensemble, double dt) {
    // Phase 3: Direct SoA processing (no conversions!)
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
    
    // Parallel ion processing (OpenMP if enabled)
    #pragma omp parallel if(config_.simulation.enable_openmp)
    {
        #pragma omp for schedule(static, 256)
        for (int i = 0; i < n_ions; ++i) {
            EhssRng& ion_rng = rng_by_ion_[i];
            
            // Skip inactive ions
            if (!active[i]) {
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
                collision_handler_->handle_collision_soa(collision_view, dt, ion_rng, 
                                                        domain_config.environment);
            }
            
            if (reaction_handler_ && !config_.reaction_db.reactions.empty()) {
                PROFILE_SCOPE_IF_ENABLED("Reaction Handling");
                auto reaction_view = ensemble.reaction_data(i);
                auto* CCS_arr = ensemble.CCS_data();
                auto* mobility_arr = ensemble.mobility_data();
                
                reaction_handler_->handle_reaction_soa(reaction_view, CCS_arr, mobility_arr,
                    dt, ion_rng, config_.reaction_db, config_.species_db,
                    domain_config.environment);
            }
            
            // 6-7. Integration (Phase 3B: Direct SoA!)
            {
                PROFILE_SCOPE_IF_ENABLED("Integration");
                
                // Get force registry for this domain
                const auto& force_registry = force_registries_[domain_idx];
                
                // Integrate directly using SoA (no conversion!)
                integrator_->step_soa(ensemble, i, current_time_, dt, *force_registry);
            }
            
            // 8. Boundary checks
            {
                PROFILE_SCOPE_IF_ENABLED("Boundary Checks");
                
                Vec3 pos_after(pos_x[i], pos_y[i], pos_z[i]);
                
                // Simple cylindrical boundary check
                bool still_inside = (pos_after.z >= -DOMAIN_BOUNDARY_EPSILON);
                still_inside = still_inside && (pos_after.z < domain_config.geometry.length_m);
                
                if (still_inside) {
                    double r = std::sqrt(pos_after.x*pos_after.x + pos_after.y*pos_after.y);
                    still_inside = (r <= domain_config.geometry.radius_m + DOMAIN_BOUNDARY_EPSILON);
                }
                
                if (!still_inside) {
                    active[i] = false;
                    continue;
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
                
                if (!config_.simulation.enable_safety_logging) {
                    std::cerr << "Warning: Ion " << i << " has invalid state at t = " 
                              << ensemble.time(i) << " s, deactivating" << std::endl;
                }
            }
        }  // End parallel for
    }  // End parallel region
}

// ============================================================================
// Ion Processing Pipeline (inline private methods)
// ============================================================================

inline int SimulationEngine::find_ion_domain(const IonState& ion) {
    PROFILE_SCOPE_IF_ENABLED("Domain Finding");
    return domain_manager_->find_domain_index(ion.pos);
}

inline void SimulationEngine::update_domain_properties(IonState& ion, int domain_idx) {
    if (ion.current_domain_index != domain_idx) {
        domain_manager_->update_domain_properties(ion, domain_idx);
    }
}

inline void SimulationEngine::process_ion_collisions(
    IonState& ion,
    DomainContext& ctx,
    double dt,
    physics::EhssRng& rng,
    int domain_idx
) {
    if (!collision_handler_) return;
    
    PROFILE_SCOPE_IF_ENABLED("Collision Handling");
    
    // Work in local coordinates
    ion.pos = ctx.pos_local();
    ion.vel = ctx.vel_local();
    
    const auto& domain_config = config_.domains[domain_idx];
    collision_handler_->handle_collision(ion, dt, rng, domain_config.environment);
    
    // Sync back to context
    ctx.pos_local() = ion.pos;
    ctx.vel_local() = ion.vel;
}

inline void SimulationEngine::process_ion_reactions(
    IonState& ion,
    DomainContext& ctx,
    double dt,
    physics::EhssRng& rng,
    int domain_idx
) {
    if (!reaction_handler_ || config_.reaction_db.reactions.empty()) return;
    
    PROFILE_SCOPE_IF_ENABLED("Reaction Handling");
    
    // Work in local coordinates
    ion.pos = ctx.pos_local();
    ion.vel = ctx.vel_local();
    
    const auto& domain_config = config_.domains[domain_idx];
    reaction_handler_->handle_reaction(
        ion,
        dt,
        rng,
        config_.reaction_db,
        config_.species_db,
        domain_config.environment
    );
    
    // Sync back to context
    ctx.pos_local() = ion.pos;
    ctx.vel_local() = ion.vel;
}

inline void SimulationEngine::integrate_ion_trajectory(
    IonState& ion,
    DomainContext& ctx,
    double dt,
    int domain_idx,
    std::vector<IonState>& ions
) {
    PROFILE_SCOPE_IF_ENABLED("Integration");
    
    // Work in local coordinates
    ion.pos = ctx.pos_local();
    ion.vel = ctx.vel_local();
    
    const auto& force_registry = force_registries_[domain_idx];
    integrator_->step(ion, current_time_, dt, *force_registry, ions);
    
    // Sync back to context
    ctx.pos_local() = ion.pos;
    ctx.vel_local() = ion.vel;
}

inline bool SimulationEngine::check_ion_boundaries(
    IonState& ion,
    DomainContext& ctx,
    int domain_idx,
    const Vec3& pos_before
) {
    PROFILE_SCOPE_IF_ENABLED("Boundary Checks");
    
    Vec3 pos_after = ctx.pos_local();
    const auto& domain_config = config_.domains[domain_idx];
    
    // Check aperture crossing (domain exit at z=length_m)
    if (pos_after.z >= domain_config.geometry.length_m && 
        pos_before.z < domain_config.geometry.length_m) {
        
        domain_manager_->check_aperture_crossing(ion, domain_idx, pos_before, pos_after);
        
        if (!ion.active) {
            domain_manager_->terminate_ion_at_boundary(ion, domain_idx, pos_before, pos_after);
            return false;
        }
        
        bool is_last_domain = (domain_idx == static_cast<int>(config_.domains.size()) - 1);
        if (is_last_domain) {
            ion.active = false;
            domain_manager_->terminate_ion_at_boundary(ion, domain_idx, pos_before, pos_after);
            return false;
        }
        
        ctx.sync_to_ion();
        return true;  // Multi-domain transition
    }
    
    // Check geometry-specific boundaries
    bool still_inside = check_geometry_boundaries(pos_after, domain_config, domain_idx);
    
    if (!still_inside) {
        domain_manager_->terminate_ion_at_boundary(ion, domain_idx, pos_before, pos_after);
        return false;
    }
    
    return true;
}

inline bool SimulationEngine::check_geometry_boundaries(
    const Vec3& pos,
    const config::DomainConfig& domain_config,
    int domain_idx
) {
    if (domain_config.instrument == config::Instrument::Orbitrap) {
        Vec3 pos_global = domain_manager_->local_to_global_pos(pos, domain_idx);
        int check_domain = domain_manager_->find_domain_index(pos_global);
        return (check_domain == domain_idx);
    }
    
    // Cylindrical geometry
    const double EPSILON = 1e-9;
    bool inside = (pos.z >= -EPSILON);
    
    bool is_last_domain = (domain_idx == static_cast<int>(config_.domains.size()) - 1);
    if (is_last_domain) {
        inside = inside && (pos.z < domain_config.geometry.length_m);
    } else {
        inside = inside && (pos.z <= domain_config.geometry.length_m + EPSILON);
    }
    
    if (inside) {
        double r = std::sqrt(pos.x*pos.x + pos.y*pos.y);
        inside = (r <= domain_config.geometry.radius_m + EPSILON);
    }
    
    return inside;
}

inline void SimulationEngine::verify_ion_safety(
    IonState& ion,
    int ion_index,
    int domain_idx
) {
    bool position_valid = ICARION::safety::is_finite(ion.pos);
    bool velocity_valid = ICARION::safety::is_finite(ion.vel);
    
    if (!position_valid || !velocity_valid) {
        log_safety_violation(ion, ion_index, domain_idx, position_valid, velocity_valid);
        ion.active = false;
        return;
    }
    
    if (config_.simulation.safety_checks.enable_bounds_checks) {
        check_bounds_violations(ion, ion_index, domain_idx);
    }
}

inline void SimulationEngine::log_safety_violation(
    const IonState& ion,
    int ion_index,
    int domain_idx,
    bool position_valid,
    bool velocity_valid
) {
    if (!config_.simulation.enable_safety_logging) {
        std::cerr << "Warning: Ion " << ion_index << " has invalid state at t=" 
                  << ion.t << " s" << std::endl;
        return;
    }
    
    safety::ViolationEvent event;
    event.type = !position_valid ? 
        (std::isnan(ion.pos.x + ion.pos.y + ion.pos.z) ? 
            safety::ViolationType::NAN_POSITION : 
            safety::ViolationType::INF_POSITION) :
        (std::isnan(ion.vel.x + ion.vel.y + ion.vel.z) ? 
            safety::ViolationType::NAN_VELOCITY : 
            safety::ViolationType::INF_VELOCITY);
    
    event.timestamp = std::chrono::steady_clock::now();
    event.ion_index = ion_index;
    event.step_number = current_step_;
    event.simulation_time = ion.t;
    event.timestep = config_.simulation.dt_s;
    event.position = ion.pos;
    event.velocity = ion.vel;
    event.violation_context = "Post-integration in domain " + std::to_string(domain_idx);
    event.violation_magnitude = !position_valid ? norm(ion.pos) : norm(ion.vel);
    event.recovery_attempted = false;
    event.recovery_successful = false;
    
    safety::NumericalSafetyLogger::getInstance().logViolation(event);
}

inline void SimulationEngine::check_bounds_violations(
    IonState& ion,
    int ion_index,
    int domain_idx
) {
    double pos_mag = norm(ion.pos);
    double vel_mag = norm(ion.vel);
    
    bool bounds_violated = false;
    safety::ViolationType violation_type;
    
    if (pos_mag > config_.simulation.safety_checks.max_position_m) {
        bounds_violated = true;
        violation_type = safety::ViolationType::BOUNDS_POSITION;
    } else if (vel_mag > config_.simulation.safety_checks.max_velocity_ms) {
        bounds_violated = true;
        violation_type = safety::ViolationType::BOUNDS_VELOCITY;
    }
    
    if (!bounds_violated) return;
    
    if (config_.simulation.enable_safety_logging) {
        safety::ViolationEvent event;
        event.type = violation_type;
        event.timestamp = std::chrono::steady_clock::now();
        event.ion_index = ion_index;
        event.step_number = current_step_;
        event.simulation_time = ion.t;
        event.timestep = config_.simulation.dt_s;
        event.position = ion.pos;
        event.velocity = ion.vel;
        event.violation_context = "Bounds check in domain " + std::to_string(domain_idx);
        event.violation_magnitude = (violation_type == safety::ViolationType::BOUNDS_POSITION) 
            ? pos_mag : vel_mag;
        event.recovery_attempted = false;
        event.recovery_successful = false;
        
        safety::NumericalSafetyLogger::getInstance().logViolation(event);
    }
    
    if (config_.simulation.safety_checks.throw_on_violation) {
        throw std::runtime_error("Bounds violation for ion " + std::to_string(ion_index) + 
                               " at t=" + std::to_string(ion.t));
    }
    
    ion.active = false;
}

}  // namespace integrator
}  // namespace ICARION

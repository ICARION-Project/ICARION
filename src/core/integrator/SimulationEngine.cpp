// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "SimulationEngine.h"
#include "DomainContext.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/BorisStrategy.h"
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

bool SimulationEngine::try_gpu_integration(std::vector<IonState>& ions, double dt) {
    PROFILE_SCOPE_IF_ENABLED("GPU Integration");
    
    // Early exit: GPU not available
    if (!gpu_helper_) {
        return false;
    }
    
    // =========================================================================
    // Smart threshold: Adjust based on integrator complexity
    // =========================================================================
    // Cache integrator type to avoid repeated dynamic_cast (expensive!)
    if (!integrator_type_cached_) {
        if (dynamic_cast<RK4Strategy*>(integrator_.get())) {
            integrator_type_ = IntegratorType::RK4;
        } else if (dynamic_cast<RK45Strategy*>(integrator_.get())) {
            integrator_type_ = IntegratorType::RK45;
        } else if (dynamic_cast<BorisStrategy*>(integrator_.get())) {
            integrator_type_ = IntegratorType::Boris;
        } else {
            integrator_type_ = IntegratorType::Unknown;
        }
        integrator_type_cached_ = true;
    }
    
    // Early exit: Unknown integrator
    if (integrator_type_ == IntegratorType::Unknown) {
        return false;
    }
    
    // Dynamic threshold based on integrator cost:
    // - Boris: 1 force eval → lower threshold (GPU beneficial at ~2000 ions)
    // - RK4:   4 force evals → standard threshold (5000 ions)
    // - RK45:  6 force evals → standard threshold (5000 ions)
    size_t effective_threshold = gpu_threshold_;
    if (integrator_type_ == IntegratorType::Boris) {
        effective_threshold = gpu_threshold_ / 2;  // Boris cheaper → lower threshold
    }
    
    int n_ions = static_cast<int>(ions.size());
    if (n_ions < static_cast<int>(effective_threshold)) {
        return false;  // Below threshold, use CPU
    }
    
    // =========================================================================
    // Extract field provider (assumes single-domain or domain 0)
    // =========================================================================
    const IFieldProvider* field_provider = extract_field_provider(0);
    
    // =========================================================================
    // Dispatch to appropriate GPU kernel (no dynamic_cast overhead!)
    // =========================================================================
    bool gpu_success = false;
    
    switch (integrator_type_) {
        case IntegratorType::RK4:
            gpu_success = gpu_helper_->integrate_batch_rk4(ions, dt, current_time_, field_provider);
            break;
            
        case IntegratorType::RK45: {
            // Get RK45 tolerance parameters
            auto* rk45 = static_cast<RK45Strategy*>(integrator_.get());  // Safe: type cached
            const auto& config = rk45->get_config();
            gpu_success = gpu_helper_->integrate_batch_rk45(
                ions, dt, current_time_, field_provider,
                config.atol, config.rtol
            );
            break;
        }
            
        case IntegratorType::Boris:
            gpu_success = gpu_helper_->integrate_batch_boris(ions, dt, current_time_, field_provider);
            break;
            
        default:
            return false;  // Should never reach here
    }
    
    // GPU integration failed → fallback to CPU
    if (!gpu_success) {
        return false;
    }
    
    // =========================================================================
    // GPU success! Update ion times (GPU only updates pos/vel, not time)
    // =========================================================================
    #pragma omp parallel for if(config_.simulation.enable_openmp)
    for (int i = 0; i < n_ions; ++i) {
        if (ions[i].active) {
            ions[i].t += dt;
        }
    }
    
    return true;  // GPU path complete
}

bool SimulationEngine::try_gpu_collisions(std::vector<IonState>& ions, double dt) {
    PROFILE_SCOPE_IF_ENABLED("GPU Collision Processing");
    
    // Early exit: GPU collision helper not available
    if (!gpu_collision_helper_) {
        return false;
    }
    
    // Count active ions for threshold check
    int n_active = 0;
    for (const auto& ion : ions) {
        if (ion.active) n_active++;
    }
    
    // Early exit: Below threshold (use active ion count, not total)
    if (n_active < static_cast<int>(gpu_collision_threshold_)) {
        return false;
    }
    
    // Need to get environment config from current domain
    // For simplicity, use domain 0's environment (assumes single-domain or uniform env)
    // TODO: Support multi-domain with different environments
    if (config_.domains.empty()) {
        return false;
    }
    
    const auto& env = config_.domains[0].environment;
    
    // Dispatch to GPU
    bool gpu_success = gpu_collision_helper_->process_collisions_batch(ions, dt, env);
    
    if (!gpu_success) {
        return false;  // GPU failed, fallback to CPU
    }
    
    return true;  // GPU collision processing complete
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

bool SimulationEngine::try_gpu_space_charge(const std::vector<IonState>& ions, std::vector<Vec3>& E_fields) {
    PROFILE_SCOPE_IF_ENABLED("GPU Space Charge");
    
    // Early exit: GPU context not available
    if (!gpu_context_) {
        return false;
    }
    
    // NOTE: This helper is not yet invoked from the main integration loop.
    // When wired in, caller must supply E_fields and handle CPU fallback.
    
    // =========================================================================
    // Threshold check: P³M only beneficial above ~1000 ions
    // =========================================================================
    // For small N, direct summation O(N²) is faster than P³M O(N log N)
    // because FFT overhead dominates. P³M is not dispatched from SimulationEngine yet.
    // =========================================================================
    
    size_t n_active = 0;
    for (const auto& ion : ions) {
        if (ion.active) ++n_active;
    }
    
    if (n_active < gpu_space_charge_threshold_) {
        return false;  // Below threshold → CPU direct summation faster
    }
    
    // =========================================================================
    // Lazy initialization: Create P³M solver on first use
    // =========================================================================
    if (!gpu_space_charge_) {
        // Auto-configure based on first domain's geometry
        // (Multi-domain space charge requires more sophisticated handling)
        if (config_.domains.empty()) {
            return false;  // No domains configured
        }
        
        const auto& domain = config_.domains[0];
        
        // Estimate domain bounds from geometry
        // TODO: Use actual SpaceChargeConfig when available
        double L = domain.geometry.length_m;
        double R = domain.geometry.radius_m;
        
        Vec3 domain_min, domain_max;
        if (L > 0 && R > 0) {
            // Use cylindrical geometry bounds
            domain_min = Vec3{-R, -R, 0.0};
            domain_max = Vec3{R, R, L};
        } else {
            // Default fallback (1cm cube)
            domain_min = Vec3{-0.01, -0.01, -0.01};
            domain_max = Vec3{0.01, 0.01, 0.01};
        }
        
        icarion::gpu::GPUSpaceChargeP3M::Config p3m_config;
        
        // Grid dimensions (adaptive based on domain size)
        // Target: ~30 µm cell size for good accuracy
        Vec3 domain_size = domain_max - domain_min;
        double target_cell_size = 3e-5;  // 30 µm
        
        p3m_config.grid_nx = std::max(32, std::min(256, static_cast<int>(domain_size.x / target_cell_size)));
        p3m_config.grid_ny = std::max(32, std::min(256, static_cast<int>(domain_size.y / target_cell_size)));
        p3m_config.grid_nz = std::max(32, std::min(256, static_cast<int>(domain_size.z / target_cell_size)));
        
        p3m_config.domain_min = domain_min;
        p3m_config.domain_max = domain_max;
        p3m_config.epsilon_0 = 8.854187817e-12;  // F/m
        
        gpu_space_charge_ = icarion::gpu::GPUSpaceChargeP3M::create(*gpu_context_, p3m_config);
        
        if (!gpu_space_charge_) {
            output_manager_->log_progress("GPU: Space Charge P³M initialization failed, using CPU");
            return false;
        }
        
        std::ostringstream msg;
        msg << "GPU: Space Charge P³M initialized (grid: " 
            << p3m_config.grid_nx << "×" << p3m_config.grid_ny << "×" << p3m_config.grid_nz
            << ", N >= " << gpu_space_charge_threshold_ << " ions)";
        output_manager_->log_progress(msg.str());
    }
    
    // =========================================================================
    // Compute space charge field on GPU
    // =========================================================================
    bool success = gpu_space_charge_->compute_space_charge_field(ions, E_fields);
    
    if (!success) {
        // GPU computation failed (out of bounds, CUDA error, etc.)
        // Fall back to CPU
        return false;
    }
    
    return true;  // GPU space charge succeeded
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
    // Try GPU acceleration - integration and collisions can succeed independently
    bool gpu_integration_done = try_gpu_integration(ions, dt);
    bool gpu_collisions_done = false;
    
    // Try GPU collisions even if integration failed (they are independent operations)
    if (collision_handler_) {
        gpu_collisions_done = try_gpu_collisions(ions, dt);
    }
    
    // If both GPU operations succeeded, we're done
    if (gpu_integration_done && (!collision_handler_ || gpu_collisions_done)) {
        return;  // Fully processed on GPU
    }
    
    // Partial GPU: Some work done on GPU, rest needs CPU
    // This is fine - the CPU loop below will skip already-completed work
#else
    bool gpu_integration_done = false;
    bool gpu_collisions_done = false;
#endif
    
    // ====================================================================
    // CPU Path (fallback or small N)
    // ====================================================================
    // Cache-blocking optimization:
    // - Process ions in fixed blocks (256) tuned for L2/cache friendliness
    // - Smaller block size reduces memory bandwidth contention
    // 
    // Thread scaling notes:
    // - Memory bandwidth saturates at ~4-8 threads (typical system)
    // - Block size chosen for cache locality (ions are independent → no false sharing)
    // - Static scheduling with blocks reduces task queue overhead
    constexpr int CACHE_BLOCK_SIZE = 256;  // Tuned for L2 cache locality
    
    // Parallel ion processing (OpenMP if enabled)
    #pragma omp parallel if(config_.simulation.enable_openmp)
    {
        #pragma omp for schedule(static, CACHE_BLOCK_SIZE)
        for (int i = 0; i < n_ions; ++i) {
            IonState& ion = ions[i];
            physics::PhysicsRng& ion_rng = rng_by_ion_[i];  // Ion-specific RNG (persistent!)
            
            if (!ion.active) continue;
            
            // ================================================================
            // 1. Find current domain
            // ================================================================
            int domain_idx = find_ion_domain(ion);
            if (domain_idx < 0) {
                ion.death_time_s = ion.t;
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
            // Skip parts already done on GPU
            // ================================================================
            if (!gpu_collisions_done) {
                process_ion_collisions(ion, ctx, dt, ion_rng, domain_idx);
            }
            process_ion_reactions(ion, ctx, dt, ion_rng, domain_idx);
            
            if (!gpu_integration_done) {
                integrate_ion_trajectory(ion, ctx, dt, domain_idx, ions);
            }
            
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
// SoA (Structure of Arrays) Implementation 
// ============================================================================

std::vector<IonState> SimulationEngine::run_soa(core::IonEnsemble& ensemble) {
    
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
    
    // Direct SoA finalization (no conversion overhead)
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
    
    // Parallel ion processing (OpenMP if enabled)
    #pragma omp parallel if(config_.simulation.enable_openmp)
    {
        #pragma omp for schedule(static, 256)
        for (int i = 0; i < n_ions; ++i) {
            PhysicsRng& ion_rng = rng_by_ion_[i];
            
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
            
            // 6-7. Integration 
            {
                PROFILE_SCOPE_IF_ENABLED("Integration");
                
                // Get force registry for this domain
                const auto& force_registry = force_registries_[domain_idx];
                
                // Integrate directly using SoA (no conversion!)
                integrator_->step_soa(ensemble, i, current_time_, dt, *force_registry);
            }
            
            // 8. Boundary checks via geometry strategy (CPU)
            {
                PROFILE_SCOPE_IF_ENABLED("Boundary Checks");
                Vec3 pos_after(pos_x[i], pos_y[i], pos_z[i]);
                if (!domain_manager_->is_inside_domain(domain_config, pos_after)) {
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
    physics::PhysicsRng& rng,
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
    physics::PhysicsRng& rng,
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
            ion.death_time_s = ion.t;  // Will be overwritten by terminate_ion_at_boundary
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
        ion.death_time_s = ion.t;
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
    
    ion.death_time_s = ion.t;
    ion.active = false;
}

}  // namespace integrator
}  // namespace ICARION

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

void SimulationEngine::initialize_soa(const core::IonEnsemble& ensemble) {
    // Single conversion at startup for metadata/output initialization
    auto ions_legacy = ensemble.to_legacy();
    initialize(ions_legacy);
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
    // Thin wrapper: convert AoS to SoA, run SoA path, return AoS
    core::IonEnsemble ensemble = core::IonEnsemble::from_legacy(ions);
    return run(ensemble);
}

std::vector<IonState> SimulationEngine::run(core::IonEnsemble& ensemble) {
    
    // 1. Initialize subsystems
    initialize_soa(ensemble);
    
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
                int new_domain_idx = domain_manager_->find_domain_index(pos_after);
                if (new_domain_idx < 0) {
                    active[i] = false;
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

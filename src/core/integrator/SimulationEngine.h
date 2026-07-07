// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file SimulationEngine.h
 * @brief Main simulation orchestrator (SSOT: Phase 5A)
 * 
 * **Responsibilities:**
 * - Main simulation loop (ion birth, integration, output)
 * - Physics module coordination (forces, collisions, reactions)
 * - Multi-domain management (DomainManager)
 * - Output management (OutputManager)
 * - Progress logging and error handling
 * 
 * **Design Principles:**
 * - SSOT: Uses config::FullConfig exclusively (no GlobalParams)
 * - Dependency injection: Physics modules passed via constructor
 * - Modular: Clear separation of concerns (domain management, output, physics)
 * - Testable: Isolated from legacy code, mockable interfaces
 * 
 * **Replaces:** integrate_trajectory() in integrator.cpp (legacy)
 */

#pragma once

#include "core/config/types/FullConfig.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"  // NEW: SoA data structure
#include "core/physics/forces/ForceRegistry.h"
#include "core/integrator/strategies/IIntegrationStrategy.h"
#include "core/physics/collisions/ICollisionHandler.h"
#include "core/types/CollisionTypes.h"  // PhysicsRng
#include "core/physics/reactions/IReactionHandler.h"
#include "core/integrator/DomainManager.h"
#include "core/integrator/OutputManager.h"
#include "core/integrator/DeepCollisionDiagnosticsTracker.h"
#include "core/integrator/SimulationEngineUtils.h"
#include <cstdint>
#include <vector>
#include <memory>

#ifdef ICARION_USE_GPU
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/spacecharge/GPUSpaceChargeP3M.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#endif

namespace ICARION {
namespace integrator {

/**
 * @class SimulationEngine
 * @brief Main trajectory integration orchestrator
 * 
 * Coordinates all simulation subsystems:
 * - Physics: ForceRegistry, IntegrationStrategy, CollisionHandler, ReactionHandler
 * - Geometry: DomainManager (multi-domain transitions)
 * - Output: OutputManager (HDF5 + text logging)
 * 
 * **Main Loop:**
 * 1. Initialize (create output files, log parameters)
 * 2. For each timestep:
 *    - Apply ion birth logic
 *    - Parallel ion processing:
 *      * Find domain (DomainManager)
 *      * Transform to local coordinates
 *      * Compute forces (ForceRegistry)
 *      * Handle collisions (ICollisionHandler)
 *      * Handle reactions (IReactionHandler)
 *      * Integrate trajectory (IIntegrationStrategy)
 *      * Check aperture crossings (DomainManager)
 *      * Transform back to global coordinates
 *    - Log trajectory snapshot (OutputManager)
 *    - Update progress logging
 * 3. Finalize (flush output, write completion metadata)
 * 
 * **Thread Safety:**
 * - Parallel ion loop (OpenMP if enabled)
 * - Thread-local RNG states
 * - Synchronized output writes (OutputManager handles locking)
 */
class SimulationEngine {
public:
    /**
     * @brief Construct simulation engine
     * @param config Simulation configuration (SSOT)
     * @param force_registries Force computation systems (one per domain)
     * @param integrator Trajectory integration strategy (RK4, RK45, Boris)
     * @param collision_handler Collision physics (optional, nullptr = no collisions)
     * @param reaction_handler Reaction chemistry (optional, nullptr = no reactions)
     * 
     * **Per-Domain ForceRegistry (Phase 12 Enhancement):**
     * - Each domain has its own ForceRegistry with domain-specific context
     * - Eliminates need to pass domain config through integration methods
     * - Improves SSOT compliance (domain config stored once in registry)
     * 
     * Physics modules are injected to enable:
     * - Testing with mock implementations
     * - Runtime algorithm selection (RK4 vs RK45 vs Boris)
     * - Modular collision model swapping (EHSS vs HSS vs OU)
     * 
     * @note force_registries.size() must equal config.domains.size()
     */
    SimulationEngine(
        const config::FullConfig& config,
        std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries,
        std::shared_ptr<IIntegrationStrategy> integrator,
        std::shared_ptr<physics::ICollisionHandler> collision_handler = nullptr,
        std::shared_ptr<physics::IReactionHandler> reaction_handler = nullptr
    );
    
    /**
     * @brief Run simulation using SoA (Structure of Arrays) data layout
     * @param ensemble Initial ion ensemble (SoA)
     * @return Final ion states (SoA)
     * 
     * **Performance Notes:**
     * - SoA loop is used internally for per-timestep work.
     * - Entry from legacy AoS calls converts once and reuses SoA thereafter.
     * - OpenMP friendliness remains (no false sharing in the hot arrays).
     * 
     * **Note:** Uses process_timestep() internally for bulk operations
     */
    core::IonEnsemble run(core::IonEnsemble& ensemble);
    
    /**
     * @brief Get simulation configuration
     */
    const config::FullConfig& get_config() const { return config_; }
    
    /**
     * @brief Get domain manager (for testing)
     */
    const DomainManager& get_domain_manager() const { return *domain_manager_; }
    
    /**
     * @brief Get output manager (for testing)
     */
    const OutputManager& get_output_manager() const { return *output_manager_; }

    uint64_t collision_events_total() const { return collision_runtime_stats_.events_total; }
    uint64_t collision_macro_attempts_total() const { return collision_runtime_stats_.macro_attempts_total; }
    uint64_t collision_substep_attempts_total() const { return collision_runtime_stats_.substep_attempts_total; }
    bool collision_monitor_complete() const { return collision_runtime_stats_.monitor_complete; }
    int steps_completed() const { return current_step_; }
    double last_wall_runtime_s() const { return last_wall_runtime_s_; }

    double mean_collisions_per_step() const {
        return collision_runtime_stats_.mean_events_per_step(current_step_);
    }

    double collision_event_fraction_per_ion_step() const {
        return collision_runtime_stats_.event_fraction_per_ion_step();
    }

    double collision_event_fraction_per_substep() const {
        return collision_runtime_stats_.event_fraction_per_substep();
    }
    
private:
    // Configuration (SSOT)
    config::FullConfig config_;
    
    // Physics modules (dependency injection)
    std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries_;  // One per domain
    std::shared_ptr<IIntegrationStrategy> integrator_;
    std::shared_ptr<physics::ICollisionHandler> collision_handler_;
    std::shared_ptr<physics::IReactionHandler> reaction_handler_;
    
    // Subsystems (owned)
    std::unique_ptr<DomainManager> domain_manager_;
    std::unique_ptr<OutputManager> output_manager_;
    bool parallel_enabled_ = false;  ///< OpenMP allowed (disabled for adaptive integrators)
    DeepCollisionDiagnosticsTracker deep_collision_diagnostics_;
    
    // Simulation state
    double current_time_ = 0.0;
    int current_step_ = 0;
    std::vector<double> dt_per_ion_;
    
    // Per-ion RNG states (persistent across timesteps!)
    std::vector<physics::PhysicsRng> rng_by_ion_;
    std::vector<uint64_t> rng_fingerprints_;
    bool space_charge_stale_warned_ = false;
    CollisionRuntimeStats collision_runtime_stats_;
    double last_wall_runtime_s_ = 0.0;
    
#ifdef ICARION_USE_GPU
    // GPU acceleration (optional)
    std::unique_ptr<icarion::gpu::GPUContext> gpu_context_;
    std::unique_ptr<icarion::gpu::GPUCollisionHelper> gpu_collision_helper_;
    std::unique_ptr<icarion::gpu::GPUSpaceChargeP3M> gpu_space_charge_;  ///< P³M space charge solver
    size_t gpu_collision_threshold_ = 5000;  ///< Minimum ions for GPU collision dispatch
    size_t gpu_space_charge_threshold_ = 1000;  ///< Minimum ions for GPU space charge
#endif
    
    /**
     * @brief Initialize simulation subsystems
     * @param ions Initial ion ensemble (for metadata)
     * 
     * Creates DomainManager and OutputManager, writes initial HDF5 metadata.
     */
    void initialize(const core::IonEnsemble& ensemble);
    
    /**
     * @brief Initialize OpenMP thread settings and NUMA awareness
     * 
     * Configures:
     * - Thread count (uses omp_get_max_threads; honors OMP_NUM_THREADS)
     * - NUMA-aware thread placement (OMP_PLACES=cores, OMP_PROC_BIND=close)
     * - Thread affinity to prevent migration
     * 
     * Called automatically during construction.
     */
    void initialize_openmp_settings();
    
#ifdef ICARION_USE_GPU
    /**
     * @brief Initialize GPU acceleration (if available and enabled)
     * @param enable_gpu Whether GPU is enabled in config
     * 
     * Attempts to create GPUContext and helper components. If GPU unavailable,
     * disabled, or initialization fails, continues with CPU-only.
     */
    void initialize_gpu(bool enable_gpu);
#endif

    /**
     * @brief Update all registered space-charge models (if any).
     *
     * Ensures each unique model recomputes its field cache once per timestep.
     */
    void update_space_charge_models(core::IonEnsemble& ensemble);

    /**
     * @brief Refresh time-varying environment parameters (e.g., pressure waveform).
     */
    void update_dynamic_environments(double t);
    
    /**
     * @brief Process one timestep using SoA (Structure of Arrays)
     * @param ensemble Ion ensemble (SoA)
     * @param dt Timestep [s]
     * @param dt_next_hint_out Suggested timestep for the next iteration (adaptive integrators)
     * @return New simulation time (max over ions)
     * 
     * **Cache-Optimized Processing:**
     * - Bulk position/velocity updates (SIMD-friendly)
     * - Zero-copy view access (IonKinematics for integration)
     * - Reduced cache misses (hot data packed together)
     * - Compatibility layer: converts to IonState for collision/reaction handlers
     * 
     */
    double process_timestep(core::IonEnsemble& ensemble);

    double perform_integration(core::IonEnsemble& ensemble,
                             double t,
                             const std::vector<double>& dt_per_ion,
                             const std::vector<int>& domain_indices,
                             std::vector<double>& dt_used_per_ion,
                             std::vector<double>& dt_next_per_ion);

    void perform_collisions(core::IonEnsemble& ensemble,
                             const std::vector<double>& dt_used_per_ion,
                             const std::vector<int>& domain_indices);

    void handle_collisions_cpu(core::IonEnsemble& ensemble,
                               const std::vector<double>& dt_used_per_ion,
                               const std::vector<size_t>& indices,
                               const config::EnvironmentConfig& env,
                               int domain_index);

    void perform_reactions(core::IonEnsemble& ensemble,
                           const std::vector<double>& dt_used_per_ion,
                           const std::vector<int>& domain_indices);

    void log_collision_runtime_stats();
    
    /**
     * @brief Log progress message (every 10%)
     * @param t Current time [s]
     * 
     * Uses OutputManager for consistent logging.
     */
    void log_progress(double t);
    
    // ========================================================================
    // Ion Processing Pipeline (private inline for performance)
    // ========================================================================
};

}  // namespace integrator
}  // namespace ICARION

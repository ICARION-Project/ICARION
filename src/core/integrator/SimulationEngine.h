// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

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
#include "core/integrator/DomainContext.h"
#include "core/integrator/OutputManager.h"
#include <vector>
#include <memory>

#ifdef ICARION_USE_GPU
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/core/GPUIntegrationHelper.h"
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
     * @brief Run simulation (legacy AoS interface, internally SoA)
     * @param ions Initial ion ensemble
     * @return Final ion states
     */
    std::vector<IonState> run(std::vector<IonState>& ions);
    
    /**
     * @brief Run simulation using SoA (Structure of Arrays) data layout
     * @param ensemble Initial ion ensemble (SoA)
     * @return Final ion states (converted back to AoS for compatibility)
     * 
     * **Performance Notes:**
     * - SoA loop is used internally for per-timestep work.
     * - Current entry/exit still convert to AoS for initialization/output; end-to-end
     *   speedup depends on overriding step_soa() without AoS copies.
     * - OpenMP friendliness remains (no false sharing in the hot arrays).
     * 
     * **Note:** Uses process_timestep_soa() internally for bulk operations
     */
    std::vector<IonState> run_soa(core::IonEnsemble& ensemble);
    
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
    
    // Simulation state
    double current_time_ = 0.0;
    int current_step_ = 0;
    
    // Per-ion RNG states (persistent across timesteps!)
    std::vector<physics::PhysicsRng> rng_by_ion_;
    
#ifdef ICARION_USE_GPU
    // GPU acceleration (optional)
    std::unique_ptr<icarion::gpu::GPUContext> gpu_context_;
    std::unique_ptr<icarion::gpu::GPUIntegrationHelper> gpu_helper_;
    std::unique_ptr<icarion::gpu::GPUCollisionHelper> gpu_collision_helper_;
    std::unique_ptr<icarion::gpu::GPUSpaceChargeP3M> gpu_space_charge_;  ///< P³M space charge solver
    size_t gpu_threshold_ = 5000;  ///< Minimum ions for GPU dispatch
    size_t gpu_collision_threshold_ = 5000;  ///< Minimum ions for GPU collision dispatch
    size_t gpu_space_charge_threshold_ = 1000;  ///< Minimum ions for GPU space charge
    
    // GPU dispatch cache (avoid repeated dynamic_cast)
    enum class IntegratorType { RK4, RK45, Boris, Unknown };
    IntegratorType integrator_type_ = IntegratorType::Unknown;
    bool integrator_type_cached_ = false;
#endif
    
    /**
     * @brief Initialize simulation subsystems
     * @param ions Initial ion ensemble (for metadata)
     * 
     * Creates DomainManager and OutputManager, writes initial HDF5 metadata.
     */
    void initialize(const std::vector<IonState>& ions);
    
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
     * Attempts to create GPUContext and GPUIntegrationHelper.
     * If GPU unavailable, disabled, or initialization fails, continues with CPU-only.
     * 
     * Called automatically during initialize().
     */
    void initialize_gpu(bool enable_gpu);
    
    /**
     * @brief Try GPU batch integration (if enabled and above threshold)
     * @param ions Ion ensemble
     * @param dt Timestep [s]
     * @return true if GPU integration succeeded, false if CPU fallback needed
     * 
     * Automatically falls back to CPU if GPU unavailable or N < threshold.
     */
    bool try_gpu_integration(std::vector<IonState>& ions, double dt);
    
    /**
     * @brief Try GPU batch collision processing (if enabled and above threshold)
     * @param ions Ion ensemble
     * @param dt Timestep [s]
     * @return true if GPU collision processing succeeded, false if CPU fallback needed
     * 
     * Automatically falls back to CPU if GPU unavailable or N < threshold.
     */
    bool try_gpu_collisions(std::vector<IonState>& ions, double dt);
    
    /**
     * @brief Try GPU space charge field computation (P³M algorithm)
     * @param ions Ion ensemble
     * @param E_fields Output: E-field at each ion position [V/m]
     * @return true if GPU space charge succeeded, false if CPU fallback needed
     * 
     * **Dispatch Logic:**
     * - Requires N >= gpu_space_charge_threshold_ (default: 1000 ions)
     * - Uses P³M algorithm (O(N log N) via FFT)
     * - SimulationEngine does not dispatch this helper yet; caller must invoke manually.
     * 
     * **Automatically falls back to CPU if:**
     * - GPU not available
     * - N < threshold (direct summation O(N²) faster for small N)
     * - Space charge config missing or disabled
     */
    bool try_gpu_space_charge(const std::vector<IonState>& ions, std::vector<Vec3>& E_fields);
    
    /**
     * @brief Extract field provider from force registry
     * @param domain_id Domain index
     * @return Pointer to field provider or nullptr if not available
     * 
     * Searches force registry for ElectricFieldForce and extracts its field provider.
     * Used by GPU integration to upload fields to texture memory.
     */
    const ::IFieldProvider* extract_field_provider(int domain_id) const;
    
    /**
     * @brief Try GPU batch boundary checking (conditional on domain config)
     * @param ensemble Ion ensemble (SoA)
     * @param domain_idx Domain index
     * @return true if GPU boundary check succeeded, false if CPU fallback needed
     * 
     * **Conditional Dispatch:**
     * - GPU boundary check only used if:
     *   1. Boundary type is Absorption (GPU doesn't support reflections yet)
     *   2. Instrument is NOT Orbitrap (GPU doesn't support hyperlogarithmic boundaries)
     * - Falls back to CPU for:
     *   - Specular/Diffuse/Thermal Reflection (requires surface normals + RNG)
     *   - Orbitrap (requires bisection-based intersection)
     * 
     * **GPU Limitations (Phase 11):**
     * - Geometry: Cylindrical only (no Orbitrap hyperlogarithmic surface)
     * - Action: Absorption only (no reflection, no thermal accommodation)
     * - Future: Phase 12 will add reflection + Orbitrap support
     * - Note: Helper is defined but not wired into the timestep loop yet.
     */
    bool try_gpu_boundary_check(core::IonEnsemble& ensemble, int domain_idx);
    
    /**
     * @brief Log GPU performance statistics
     * 
     * Called at simulation end to report GPU usage and speedup.
     */
    void finalize_gpu();
#endif
    
    /**
     * @brief Process one timestep for all ions (legacy AoS)
     * @param ions Ion ensemble
     * @param dt Timestep [s]
     * 
     * Parallel ion loop with domain management, physics, and integration.
     */
    void process_timestep(std::vector<IonState>& ions, double dt);
    
    /**
     * @brief Process one timestep using SoA (Structure of Arrays)
     * @param ensemble Ion ensemble (SoA)
     * @param dt Timestep [s]
     * 
     * **Cache-Optimized Processing:**
     * - Bulk position/velocity updates (SIMD-friendly)
     * - Zero-copy view access (IonKinematics for integration)
     * - Reduced cache misses (hot data packed together)
     * - Compatibility layer: converts to IonState for collision/reaction handlers
     * 
     * **Performance:** Gains materialize when integrators override step_soa() and
     * avoid AoS conversions; default wrappers may limit gains.
     */
    void process_timestep_soa(core::IonEnsemble& ensemble, double dt);
    
    /**
     * @brief Apply ion birth logic (delayed emission)
     * @param ions Ion ensemble
     * @param t Current simulation time [s]
     * 
     * Activates ions with birth_time_s <= t.
     */
    void apply_ion_birth(std::vector<IonState>& ions, double t);
    
    /**
     * @brief Check if simulation should continue
     * @param ions Ion ensemble
     * @param t Current time [s]
     * @return true if simulation should continue
     * 
     * Stops if all ions inactive or time exceeded.
     */
    bool should_continue(const std::vector<IonState>& ions, double t) const;
    
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
    
    /**
     * @brief Find domain index for ion position
     * @return Domain index, or -1 if outside all domains
     */
    inline int find_ion_domain(const IonState& ion);
    
    /**
     * @brief Update ion properties when entering new domain
     */
    inline void update_domain_properties(IonState& ion, int domain_idx);
    
    /**
     * @brief Apply collision effects to single ion
     * @param ctx DomainContext (manages coordinate transforms)
     * @param rng Ion-specific RNG (thread-safe)
     */
    inline void process_ion_collisions(
        IonState& ion,
        DomainContext& ctx,
        double dt,
        physics::PhysicsRng& rng,
        int domain_idx
    );
    
    /**
     * @brief Apply reaction effects to single ion
     * @param ctx DomainContext (manages coordinate transforms)
     * @param rng Ion-specific RNG (thread-safe)
     */
    inline void process_ion_reactions(
        IonState& ion,
        DomainContext& ctx,
        double dt,
        physics::PhysicsRng& rng,
        int domain_idx
    );
    
    /**
     * @brief Integrate ion trajectory (RK4/RK45/Boris)
     * @param ctx DomainContext (manages coordinate transforms)
     * @param ions Full ion ensemble (needed by integrator for space charge)
     */
    inline void integrate_ion_trajectory(
        IonState& ion,
        DomainContext& ctx,
        double dt,
        int domain_idx,
        std::vector<IonState>& ions
    );
    
    /**
     * @brief Check if ion crossed boundaries (aperture, walls)
     * @param pos_before Position before integration
     * @return true if ion is still inside domain
     */
    inline bool check_ion_boundaries(
        IonState& ion,
        DomainContext& ctx,
        int domain_idx,
        const Vec3& pos_before
    );
    
    /**
     * @brief Check geometry-specific boundaries (cylindrical vs Orbitrap)
     */
    inline bool check_geometry_boundaries(
        const Vec3& pos,
        const config::DomainConfig& domain_config,
        int domain_idx
    );
    
    /**
     * @brief Verify numerical safety (NaN, Inf, bounds)
     */
    inline void verify_ion_safety(
        IonState& ion,
        int ion_index,
        int domain_idx
    );
    
    /**
     * @brief Log safety violation event
     */
    inline void log_safety_violation(
        const IonState& ion,
        int ion_index,
        int domain_idx,
        bool position_valid,
        bool velocity_valid
    );
    
    /**
     * @brief Check bounds violations (position/velocity magnitude)
     */
    inline void check_bounds_violations(
        IonState& ion,
        int ion_index,
        int domain_idx
    );
};

}  // namespace integrator
}  // namespace ICARION

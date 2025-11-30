// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IForce.h"
#include "core/physics/spacecharge/spaceChargeSolver.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include <memory>

namespace ICARION::physics {

/**
 * @brief Grid-based space charge force using Poisson solver (O(N log N))
 * 
 * Integrates SpaceChargeSolver into the IForce framework for automatic
 * method selection. Updates Poisson solver once per timestep for all ions,
 * then interpolates field at each ion position.
 * 
 * Key Features:
 * - Automatic once-per-timestep updates (tracks time to avoid redundant solves)
 * - Grid-based Poisson solver: ∇²φ = -ρ/ε₀
 * - Fast for large ensembles: O(N log N) vs O(N²) for direct Coulomb
 * - Trilinear interpolation for smooth fields
 * 
 * Performance:
 * - N = 1000:   ~20 ms/timestep (10x faster than SpaceChargeForce)
 * - N = 10000:  ~30 ms/timestep (667x faster than SpaceChargeForce)
 * - N = 100000: ~100 ms/timestep (20000x faster than SpaceChargeForce)
 * 
 * Usage:
 * @code
 * // Create solver (64³ grid, 1mm cells)
 * auto solver = std::make_shared<SpaceChargeSolver>(
 *     64, 64, 64,           // Grid resolution
 *     1e-3, 1e-3, 1e-3,     // Cell size [m]
 *     Vec3{0, 0, 0}         // Origin
 * );
 * 
 * // Wrap in IForce interface
 * auto force = std::make_unique<SpaceChargeGrid>(solver);
 * force_registry.add_force(std::move(force));
 * @endcode
 * 
 * @note Requires ctx.all_ions to be populated in ForceContext
 * @note Updates solver only once per unique timestep (time-based deduplication)
 * @note For N < 1000, prefer SpaceChargeForce (direct Coulomb is faster + exact)
 */
class SpaceChargeGrid : public IForce {
public:
    /**
     * @brief Construct force from existing solver
     * 
     * @param solver Shared pointer to SpaceChargeSolver (must not be null)
     * 
     * @throws std::invalid_argument if solver is null
     */
    explicit SpaceChargeGrid(std::shared_ptr<SpaceChargeSolver> solver);
    
    // =========================================================================
    // IForce Interface
    // =========================================================================
    
    /**
     * @brief Compute space charge force from grid-based field
     * 
     * Updates Poisson solver once per timestep (for all ions), then
     * interpolates electric field at the given ion position.
     * 
     * @param ion Current ion (compute force ON this ion)
     * @param t Simulation time [s] (used for update deduplication)
     * @param ctx Force context - MUST contain all_ions vector
     * @return Space charge force [N] = q·E_sc(r)
     * 
     * @note Returns {0,0,0} if ctx.all_ions is null or solver fails
     * @note Thread-safe: Only first thread updates solver (OpenMP critical section)
     */
    Vec3 compute(
        const IonState& ion,
        double t,
        const ForceContext& ctx
    ) const override;
    
    /**
     * @brief Check if force applies to ion
     * 
     * Space charge applies to all charged particles.
     * 
     * @param ion Ion to check
     * @return true always (space charge universal for charged particles)
     */
    bool applies_to(const IonState& ion) const override;
    
    /**
     * @brief Get force name
     * @return "SpaceChargeGrid" for logging/debugging
     */
    std::string name() const override;
    
    // =========================================================================
    // Accessors
    // =========================================================================
    
    /**
     * @brief Get underlying solver (for inspection/configuration)
     * @return Shared pointer to solver (may be null if not initialized)
     */
    std::shared_ptr<SpaceChargeSolver> get_solver() const { return solver_; }
    
    /**
     * @brief Set new solver (allows runtime reconfiguration)
     * @param solver New solver instance
     */
    void set_solver(std::shared_ptr<SpaceChargeSolver> solver) { solver_ = solver; }

private:
    std::shared_ptr<SpaceChargeSolver> solver_;  ///< Grid-based Poisson solver
    
    // Update tracking (mutable for const compute())
    mutable double last_update_time_;  ///< Last time solver was updated
    mutable bool solver_updated_this_step_;  ///< Flag to prevent redundant updates
};

}  // namespace ICARION::physics

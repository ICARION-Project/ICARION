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
 * @brief Grid-based space charge force using Poisson solver (CPU box grid)
 * 
 * Wraps the CPU SpaceChargeSolver (rectangular grid) as an IForce. Updates the
 * solver once per timestep and interpolates the resulting field at ion positions.
 * Geometry awareness is limited to the box grid; not suitable for cylindrical/Orbitrap.
 * No validated performance numbers; use for rough estimates only.
 * 
 * @note Requires ctx.all_ions to be populated. Updates are deduped by timestep.
 * @note For small N, the direct O(N²) force is exact and may be faster.
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

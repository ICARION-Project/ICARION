// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#ifdef ICARION_USE_GPU

#include "IForce.h"
#include "core/gpu/spacecharge/GPUSpaceChargeP3M.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include <memory>
#include <vector>

namespace ICARION::physics {

/**
 * @brief GPU space charge force wrapper (experimental P³M)
 * 
 * Wraps GPUSpaceChargeP3M as an IForce. GPU solver is experimental, rectangular-only,
 * not validated, and not geometry-aware (no cylindrical/Orbitrap). Use only for
 * exploratory runs; CPU direct remains the exact reference.
 */
class SpaceChargeGPU : public IForce {
public:
    /**
     * @brief Construct GPU space charge force from P³M solver
     * 
     * @param solver Shared pointer to GPUSpaceChargeP3M (must not be null)
     * 
     * @throws std::invalid_argument if solver is null
     */
    explicit SpaceChargeGPU(std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> solver);
    
    // =========================================================================
    // IForce Interface
    // =========================================================================
    
    /**
     * @brief Compute space charge force from GPU P³M field
     * 
     * Updates P³M solver once per timestep (for all ions), then
     * returns precomputed electric field at the given ion position.
     * 
     * @param ion Current ion (compute force ON this ion)
     * @param t Simulation time [s] (used for update deduplication)
     * @param ctx Force context - MUST contain all_ions vector
     * @return Space charge force [N] = q·E_sc(r)
     * 
     * @note Returns {0,0,0} if ctx.all_ions is null or solver fails
     * @note Thread-safe: Only first thread updates solver (OpenMP critical section)
     * @note E-field is cached after update, so subsequent calls are O(1)
     */
    Vec3 compute(
        const IonState& ion,
        double t,
        const ForceContext& ctx
    ) const override;

    /**
     * @brief SoA-aware computation (uses ctx.ion_ensemble if provided)
     */
    Vec3 compute_soa(
        const core::IonEnsemble& ensemble,
        size_t ion_idx,
        double t,
        const ForceContext& ctx
    ) const override;
    
    /**
     * @brief Check if force applies to this ion
     * 
     * Space charge applies to all ions (Coulomb interaction is universal).
     * Individual ion activity/birth checks are handled in compute().
     * 
     * @param ion Ion to check
     * @return Always true (all charged particles feel space charge)
     */
    bool applies_to(const IonState& ion) const override;
    
    /**
     * @brief Get force name for logging/debugging
     * @return "SpaceChargeGPU"
     */
    std::string name() const override;
    
    // =========================================================================
    // GPU-specific API
    // =========================================================================
    
    /**
     * @brief Get underlying P³M solver (for diagnostics/tuning)
     * @return Shared pointer to GPU solver
     */
    std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> get_solver() const { return solver_; }
    
    /**
     * @brief Check if solver was updated this timestep
     * @return True if solver.compute_space_charge_field() was called at current time
     */
    bool was_updated_this_step() const { return solver_updated_this_step_; }

private:
    std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> solver_;  ///< GPU P³M solver
    
    // Timestep tracking (for once-per-timestep updates)
    mutable double last_update_time_;        ///< Last time solver was updated [s]
    mutable bool solver_updated_this_step_;  ///< Flag: was solver updated this timestep?
    
    // Cached E-fields (updated once per timestep, indexed by ion position in all_ions)
    mutable std::vector<Vec3> cached_E_fields_;  ///< E-field at each ion position [V/m]

    /**
     * @brief Ensure the GPU solver and cache are updated for the given timestep
     */
    void update_fields_if_needed(
        const core::IonEnsemble& ensemble,
        double t
    ) const;
};

}  // namespace ICARION::physics

#endif  // ICARION_USE_GPU

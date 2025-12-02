// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IForce.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "utils/constants.h"

namespace ICARION::physics {

/**
 * @brief Space charge force via direct N-body Coulomb interactions (CPU)
 * 
 * Computes electrostatic interactions in O(N²) over the ensemble with optional
 * softening to avoid singularities. No geometry masking or cutoffs are applied.
 * Intended for small N or verification; use with caution for cylindrical/Orbitrap
 * setups where this point-charge model does not respect boundaries.
 */
class SpaceChargeDirect : public IForce {
public:
    /**
     * @brief Construct space charge force
     * 
     * @param softening_m Softening length scale [m] to prevent 1/r² divergence
     *                    Typical: 1e-10 m (0.1 nm) - much smaller than ion spacing
     *                    Set to 0.0 to disable softening (use with caution)
     */
    explicit SpaceChargeDirect(double softening_m = 1e-10);
    
    // =========================================================================
    // IForce Interface
    // =========================================================================
    
    /**
     * @brief Compute total Coulomb force from all other ions
     * 
     * Iterates over all ions in ctx.all_ions and sums pairwise forces.
     * Excludes self-interaction (ion.id == other_ion.id).
     * 
     * @param ion Current ion (computes force on THIS ion)
     * @param t Simulation time [s] (unused - Coulomb force is time-independent)
     * @param ctx Force context - MUST contain all_ions vector
     * @return Total Coulomb force [N] from all other ions
     * 
     * @note Returns {0,0,0} if ctx.all_ions is null or empty
     * @note Returns {0,0,0} if ion_ensemble.size() == 1 (no other ions)
     */
    Vec3 compute(
        const IonState& ion,
        double t,
        const ForceContext& ctx
    ) const override;
    
    /**
     * @brief Check if space charge applies
     * 
     * Space charge only applies if there are other ions (N-body interactions).
     * Always returns true here (actual check happens in compute()).
     * 
     * @param ion Ion to check (unused)
     * @return true (space charge can apply to any ion)
     */
    bool applies_to(const IonState& ion) const override;
    
    /**
     * @brief Get force name
     * @return "SpaceCharge" for logging/debugging
     */
    std::string name() const override;
    
    // =========================================================================
    // Accessors
    // =========================================================================
    
    /**
     * @brief Get softening parameter
     * @return Softening length [m]
     */
    double get_softening() const { return softening_m_; }
    
    /**
     * @brief Set softening parameter
     * @param softening_m New softening length [m]
     * 
     * @note For typical ion clouds (spacing ~1 µm), use ε = 0.01-0.1 nm
     * @note Setting ε = 0 disables softening (use with caution - can diverge)
     */
    void set_softening(double softening_m) { softening_m_ = softening_m; }

private:
    // =========================================================================
    // Helper Methods
    // =========================================================================
    
    /**
     * @brief Compute pairwise Coulomb force
     * 
     * F = k_e · q₁ · q₂ · r̂ / (r² + ε²)^(3/2)
     * 
     * Softening modification:
     *   - Standard Coulomb: F ∝ 1/r²
     *   - Softened Coulomb: F ∝ 1/(r² + ε²)^(3/2)
     *   - As r → 0: F → k·q₁·q₂/ε³ (finite, not infinite)
     *   - As r → ∞: F → k·q₁·q₂/r² (recovers standard Coulomb)
     * 
     * @param ion1 First ion (force computed ON this ion)
     * @param ion2 Second ion (force computed FROM this ion)
     * @return Force [N] on ion1 due to ion2
     * 
     * @note Returns {0,0,0} if ions overlap exactly (r = 0)
     */
    Vec3 compute_pairwise_force(
        const IonState& ion1,
        const IonState& ion2
    ) const;
    
    // =========================================================================
    // Member Variables
    // =========================================================================
    
    /// Softening length [m] to prevent 1/r² divergence at small distances
    double softening_m_;
};

} // namespace ICARION::physics

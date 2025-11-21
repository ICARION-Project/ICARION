// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#pragma once

#include "IForce.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"

namespace ICARION::physics {

/**
 * @brief Space charge force (N-body Coulomb interactions)
 * 
 * Computes electrostatic repulsion/attraction between ions in the ensemble.
 * This is the most computationally expensive force (O(N²) for N ions).
 * 
 * Physics:
 *   Coulomb force between two point charges:
 *   F = k_e · q₁ · q₂ · r̂ / r²
 * 
 *   where:
 *   - k_e = 8.987551787e9 N·m²/C² (Coulomb constant)
 *   - q₁, q₂ = charges [C]
 *   - r = distance between ions [m]
 *   - r̂ = unit vector from q₂ to q₁
 * 
 * Key Features:
 *   - Self-interaction exclusion (ion doesn't interact with itself)
 *   - Softening parameter for numerical stability (avoids 1/0 at r→0)
 *   - Newton's 3rd law: F₁₂ = -F₂₁ (symmetric interactions)
 *   - Optional spatial partitioning for performance (future: O(N log N))
 * 
 * Usage:
 * @code
 * // Create space charge force with softening
 * auto sc_force = std::make_unique<SpaceChargeForce>(1e-10);  // 0.1 nm softening
 * 
 * // Compute force requires all_ions in context
 * ForceContext ctx;
 * ctx.all_ions = &ion_ensemble;
 * Vec3 force = sc_force->compute(ion, t, ctx);
 * @endcode
 * 
 * Performance Considerations:
 *   - O(N²) scaling: 1000 ions → 500k pairs → ~1ms per timestep
 *   - Spatial partitioning (future): reduces to O(N log N) or O(N)
 *   - OpenMP parallelization recommended for N > 100
 *   - Consider cutoff radius for long-range interactions (future)
 * 
 * Numerical Stability:
 *   - Softening parameter prevents divergence at r → 0
 *   - Modified Coulomb: r_eff = √(r² + ε²)
 *   - Typical ε: 0.01-0.1 nm (much smaller than ion spacing)
 * 
 * @note Requires ctx.all_ions to be non-null, otherwise returns zero force
 * @note For single-ion systems, always returns zero force (no pairs)
 */
class SpaceChargeForce : public IForce {
public:
    /**
     * @brief Construct space charge force
     * 
     * @param softening_m Softening length scale [m] to prevent 1/r² divergence
     *                    Typical: 1e-10 m (0.1 nm) - much smaller than ion spacing
     *                    Set to 0.0 to disable softening (use with caution)
     */
    explicit SpaceChargeForce(double softening_m = 1e-10);
    
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
    // Constants
    // =========================================================================
    
    /// Coulomb constant: k_e = 1/(4πε₀) = 8.987551787e9 N·m²/C²
    static constexpr double k_coulomb_ = 8.987551787e9;
    
    // =========================================================================
    // Member Variables
    // =========================================================================
    
    /// Softening length [m] to prevent 1/r² divergence at small distances
    double softening_m_;
};

} // namespace ICARION::physics

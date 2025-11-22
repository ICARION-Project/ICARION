// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "SpaceChargeForce.h"
#include "ForceContext.h"
#include "core/utils/mathUtils.h"
#include "utils/constants.h"

#include <cmath>
#include <algorithm>

namespace ICARION::physics {

// ============================================================================
// Constructor
// ============================================================================

SpaceChargeForce::SpaceChargeForce(double softening_m)
    : softening_m_(softening_m)
{
    if (softening_m_ < 0.0) {
        throw std::invalid_argument(
            "SpaceChargeForce: Softening parameter must be non-negative. "
            "Got: " + std::to_string(softening_m_) + " m"
        );
    }
}

// ============================================================================
// IForce Interface Implementation
// ============================================================================

Vec3 SpaceChargeForce::compute(
    const IonState& ion,
    double t,
    const ForceContext& ctx
) const {
    (void)t;  // Time-independent force
    
    // Early exit: no ion ensemble provided
    if (!ctx.all_ions || ctx.all_ions->empty()) {
        return Vec3{0.0, 0.0, 0.0};
    }
    
    // Early exit: only one ion (no pairwise interactions)
    if (ctx.all_ions->size() == 1) {
        return Vec3{0.0, 0.0, 0.0};
    }
    
    // Sum pairwise forces from all other ions
    Vec3 total_force{0.0, 0.0, 0.0};
    
    for (const IonState& other_ion : *ctx.all_ions) {
        // Skip self-interaction (compare by address)
        if (&ion == &other_ion) {
            continue;
        }
        
        // Compute and accumulate pairwise force
        Vec3 pairwise_force = compute_pairwise_force(ion, other_ion);
        total_force = total_force + pairwise_force;
    }
    
    return total_force;
}

bool SpaceChargeForce::applies_to(const IonState& ion) const {
    (void)ion;  // Applies to all ions
    return true;
}

std::string SpaceChargeForce::name() const {
    return "SpaceCharge";
}

// ============================================================================
// Helper Methods
// ============================================================================

Vec3 SpaceChargeForce::compute_pairwise_force(
    const IonState& ion1,
    const IonState& ion2
) const {
    // Vector from ion2 to ion1: r̂ = (r₁ - r₂) / |r₁ - r₂|
    const Vec3 r_vec = ion1.pos - ion2.pos;
    const double r_sq = r_vec.x * r_vec.x 
                      + r_vec.y * r_vec.y 
                      + r_vec.z * r_vec.z;
    
    // Handle exact overlap (r = 0)
    // Use softening threshold: if r < softening, skip (ions too close)
    // Minimum threshold: 1% of softening, or 1e-30 m² (numerical safety for softening=0)
    const double min_distance_sq = std::max(softening_m_ * softening_m_ * 1e-4, 1e-30);
    if (r_sq < min_distance_sq) {
        return Vec3{0.0, 0.0, 0.0};  // Skip if inside minimum safe distance
    }
    
    // Softened distance: r_eff = √(r² + ε²)
    const double epsilon_sq = softening_m_ * softening_m_;
    const double r_eff_sq = r_sq + epsilon_sq;
    const double r_eff = std::sqrt(r_eff_sq);
    
    // Unit vector: r̂ = r_vec / r
    const double r = std::sqrt(r_sq);
    const Vec3 r_hat = r_vec / r;
    
    // Coulomb force magnitude:
    // Standard Coulomb: F = k · q₁ · q₂ / r²
    //
    // With softening (ε > 0):
    //   Softened potential: φ(r) = k·q / √(r² + ε²)
    //   Softened force: F = -∇φ · r̂ = k·q₁·q₂ / (r² + ε²) · r̂
    //
    //   Note: NOT (r²+ε²)^(3/2)! The gradient of 1/√(r²+ε²) gives 1/(r²+ε²), not 1/(r²+ε²)^(3/2)
    //
    // This ensures:
    //   - r → ∞: F → k·q₁·q₂/r² (recovers standard Coulomb)
    //   - r → 0:  F → k·q₁·q₂/ε² (finite, not infinite)
    //   Coulomb constant: COULOMB_CONST from utils/constants.h defined as k = 1/(4πε₀)
    const double force_magnitude = COULOMB_CONST * ion1.ion_charge_C * ion2.ion_charge_C 
                                 / r_eff_sq;
    
    // Force vector: F⃗ = F_mag · r̂
    return r_hat * force_magnitude;
}

} // namespace ICARION::physics

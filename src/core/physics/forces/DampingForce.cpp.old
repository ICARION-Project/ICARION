// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "DampingForce.h"
#include "utils/constants.h"

#include <cmath>

namespace ICARION {
namespace physics {

// ============================================================================
// Constructor
// ============================================================================

DampingForce::DampingForce(const DampingParams& params)
    : params_(params)
    , rng_(params.random_seed)
    , normal_dist_(0.0, 1.0)  // Standard normal N(0,1)
{
    // No validation - all params can be zero (disabled force)
}

// ============================================================================
// IForce Interface Implementation
// ============================================================================

Vec3 DampingForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    (void)t;  // Time-independent (except for Langevin random kicks)
    
    if (params_.model == DampingModel::None) {
        return Vec3{0.0, 0.0, 0.0};
    }
    
    // Friction term: F_friction = -γ·v
    Vec3 F_friction = ion.vel * (-params_.damping_coefficient);
    
    if (params_.model == DampingModel::Friction) {
        return F_friction;
    }
    
    // Langevin: F = -γ·v + ξ(t)
    if (params_.model == DampingModel::Langevin) {
        // Use temperature from context if available, otherwise from params
        double T = ctx.temperature_K > 0.0 ? ctx.temperature_K : params_.temperature_K;
        
        // Estimate time step from context (if available)
        // TODO: Pass dt explicitly in ForceContext for proper noise scaling
        double dt = 1e-9;  // Default: 1 ns (placeholder)
        
        Vec3 F_random = compute_random_force(ion.mass_kg, T, dt);
        
        return F_friction + F_random;
    }
    
    return Vec3{0.0, 0.0, 0.0};
}

std::string DampingForce::name() const {
    switch (params_.model) {
        case DampingModel::Friction: return "Damping(Friction)";
        case DampingModel::Langevin: return "Damping(Langevin)";
        default:                     return "Damping(None)";
    }
}

// ============================================================================
// Langevin Random Force
// ============================================================================

Vec3 DampingForce::compute_random_force(double ion_mass, double temperature, double dt) const {
    // Fluctuation-dissipation theorem:
    // ⟨F_random²⟩ = 2·γ·k_B·T·m / Δt
    // 
    // For each component: F_i = σ·ξ_i, where ξ_i ~ N(0,1)
    // σ = sqrt(2·γ·k_B·T·m / Δt)
    
    const double variance = 2.0 * params_.damping_coefficient * BOLTZMANN_CONSTANT 
                          * temperature * ion_mass / dt;
    const double sigma = std::sqrt(variance);
    
    // Generate three independent Gaussian random numbers
    Vec3 random_force;
    random_force.x = sigma * normal_dist_(rng_);
    random_force.y = sigma * normal_dist_(rng_);
    random_force.z = sigma * normal_dist_(rng_);
    
    return random_force;
}

} // namespace physics
} // namespace ICARION

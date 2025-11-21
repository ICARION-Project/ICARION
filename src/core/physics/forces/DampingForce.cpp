// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "DampingForce.h"
#include "utils/constants.h"

#include <cmath>
#include <algorithm>

namespace ICARION {
namespace physics {

// ============================================================================
// Constructor
// ============================================================================

DampingForce::DampingForce(const DampingParams& params)
    : params_(params)
{
    // No validation - all params can be zero (disabled force)
}

// ============================================================================
// IForce Interface Implementation
// ============================================================================

Vec3 DampingForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    (void)t;  // Time-independent (deterministic damping)
    
    if (params_.model == DampingModel::None) {
        return Vec3{0.0, 0.0, 0.0};
    }
    
    // Calculate damping coefficient γ [1/s]
    double gamma = calculate_gamma(ion, ctx);
    
    if (gamma <= 0.0) {
        return Vec3{0.0, 0.0, 0.0};
    }
    
    // Damping force: F = -γ·m·v [N]
    // (Note: Force, not acceleration - will be divided by mass in integrator)
    return ion.vel * (-gamma * ion.mass_kg);
}

std::string DampingForce::name() const {
    switch (params_.model) {
        case DampingModel::HardSphere: return "Damping(HardSphere)";
        case DampingModel::Langevin:   return "Damping(Langevin)";
        case DampingModel::Friction:   return "Damping(Friction)";
        default:                       return "Damping(None)";
    }
}

// ============================================================================
// Damping Coefficient Calculation
// ============================================================================

double DampingForce::calculate_gamma(const IonState& ion, const ForceContext& ctx) const {
    // If explicit gamma provided, use it directly
    if (params_.gamma_coefficient > 0.0) {
        return params_.gamma_coefficient;
    }
    
    // Extract gas properties from context (domain environment)
    const double gas_density = params_.gas_density_m3 > 0.0 
                             ? params_.gas_density_m3 
                             : ctx.gas_density_m3;
    
    const double v_th = params_.mean_thermal_velocity_m_s > 0.0
                      ? params_.mean_thermal_velocity_m_s
                      : ctx.mean_thermal_velocity_m_s;
    
    const double m_neutral = params_.neutral_mass_kg > 0.0
                           ? params_.neutral_mass_kg
                           : ctx.neutral_mass_kg;
    
    // ========================================================================
    // Model-specific damping coefficient calculation
    // ========================================================================
    
    switch (params_.model) {
        case DampingModel::HardSphere: {
            // Hard-sphere collision frequency:
            // γ = ν_collision = n·σ·v_th·(m_i/(m_n+m_i))
            // Note: Legacy uses m_i in numerator (ion mass factor)
            
            const double CCS = params_.CCS_m2 > 0.0 ? params_.CCS_m2 : ion.CCS_m2;
            const double m_ion = ion.mass_kg;
            
            if (CCS <= 0.0 || gas_density <= 0.0 || v_th <= 0.0 || m_neutral <= 0.0) {
                return 0.0;  // Missing parameters
            }
            
            // gamma = n·σ·v_th · m_i/(m_n+m_i)
            const double gamma = gas_density * CCS * v_th 
                                                * m_ion / (m_neutral + m_ion);
            
            return gamma;  // γ = ν_collision [1/s]
        }
        
        case DampingModel::Langevin: {
            // Langevin collision frequency (velocity-dependent):
            // γ = ν_Langevin = n·σ_L(v)·v_th·m_reduced/m_ion
            // where σ_L = π·q·√(α/(4πε₀·m_reduced))/|v|
            
            const double alpha = params_.neutral_polarizability_m3 > 0.0
                               ? params_.neutral_polarizability_m3
                               : ctx.neutral_polarizability_m3;
            
            const double m_ion = ion.mass_kg;
            const double q = ion.ion_charge_C;
            
            if (alpha <= 0.0 || gas_density <= 0.0 || v_th <= 0.0 || m_neutral <= 0.0) {
                return 0.0;  // Missing parameters
            }
            
            // Ion velocity magnitude (guard against division by zero)
            const double v_mag = std::max(1e-6, norm(ion.vel));
            
            // Reduced mass
            const double m_reduced = (m_ion * m_neutral) / (m_ion + m_neutral);
            
            // Langevin cross-section (velocity-dependent)
            // σ_L = π·q·√(α/(4πε₀·m_reduced))/v
            const double epsilon_factor = 4.0 * M_PI * EPSILON_0;
            const double cs = M_PI * q * std::sqrt(alpha / (epsilon_factor * m_reduced)) / v_mag;
            
            // Collision frequency
            const double gamma = gas_density * cs * v_th * m_neutral / (m_neutral + m_ion);
            
            return gamma;  // γ = ν_Langevin [1/s]
        }
        
        case DampingModel::Friction: {
            // Mobility-based friction:
            // γ = q/(K₀·m_ion) where K₀ = reduced mobility
            
            const double K0_cm2_Vs = params_.reduced_mobility_cm2_Vs > 0.0
                                   ? params_.reduced_mobility_cm2_Vs
                                   : ion.reduced_mobility_cm2_Vs;
            
            const double m_ion = ion.mass_kg;
            const double q = ion.ion_charge_C;
            
            if (K0_cm2_Vs <= 0.0 || m_ion <= 0.0) {
                return 0.0;  // Missing parameters
            }
            
            // Convert K₀ from cm²/(V·s) to m²/(V·s)
            const double K0_m2_Vs = K0_cm2_Vs * 1e-4;
            
            // Ion mobility at current gas density
            // K = K₀·(n/n₀) where n₀ = Loschmidt constant
            const double ion_mobility = K0_m2_Vs * LOSCHMIDT_CONSTANT / gas_density;
            
            // Friction coefficient: γ = q/(K·m)
            const double gamma = q / (ion_mobility * m_ion);
            
            return gamma;  // γ [1/s]
        }
        
        default:
            return 0.0;
    }
}

} // namespace physics
} // namespace ICARION

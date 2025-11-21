// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file DampingForce.h
 * @brief Deterministic collision damping forces (HardSphere, Langevin, Friction)
 * 
 * Computes velocity-dependent damping forces F = -γ·v for ions in background gas.
 * All models are DETERMINISTIC - they compute continuous friction opposing motion.
 * 
 * The damping coefficient γ differs by collision model:
 * - **HardSphere**: γ = ν_collision = n·σ·v_th·(m_n/(m_i+m_n))
 * - **Langevin**: γ = ν_Langevin(v) (velocity-dependent, polarization)
 * - **Friction**: γ = q/K₀ (mobility-based)
 * 
 * @note RANDOM THERMAL KICKS ARE NOT COMPUTED HERE!
 * Stochastic diffusion (Ornstein-Uhlenbeck process) is handled separately
 * by CollisionEngine via apply_ou_velocity_kick() when enable_ou_thermalization=true.
 * This separation matches the legacy architecture:
 *   1. ODE solver uses deterministic damping forces
 *   2. After ODE step, apply_ou_velocity_kick() adds thermal noise
 */

#pragma once

#include "IForce.h"
#include "ForceContext.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"

namespace ICARION {
namespace physics {

/**
 * @enum DampingModel
 * @brief Collision model for damping force calculation
 * 
 * Each model computes F = -γ·m·v with different formulas for γ:
 * - HardSphere: Elastic collisions with momentum transfer rate
 * - Langevin: Long-range ion-neutral polarization interactions
 * - Friction: Mobility-based drag (Mason-Schamp equation)
 */
enum class DampingModel {
    None,       ///< No damping
    HardSphere, ///< Hard-sphere elastic collisions
    Langevin,   ///< Langevin polarization model
    Friction    ///< Mobility-based friction
};

// ============================================================================
// ⚠️ DEPRECATED: DampingParams violates SSOT principle!
// ============================================================================
// This struct duplicates parameters from FullConfig → DomainConfig → EnvironmentConfig.
// 
// **TODO (Phase 2):** Replace with direct EnvironmentConfig reference:
//   DampingForce(const EnvironmentConfig& env, DampingModel model)
// **KEPT FOR NOW:** To avoid breaking changes during Phase 1.
// ============================================================================

/**
 * @brief Parameters for damping force calculation
 * 
 * @deprecated Violates SSOT. Use EnvironmentConfig directly in Phase 2.
 * 
 * Contains collision parameters needed to compute damping coefficient γ.
 * Different models use different subsets of these parameters.
 */
struct DampingParams {
    DampingModel model = DampingModel::None;
    
    // --- Explicit damping coefficient (if > 0, overrides model calculation) ---
    double gamma_coefficient = 0.0;  ///< Friction coefficient γ [1/s] (F = -γ·m·v)
    
    // --- HardSphere model parameters ---
    double gas_density_m3 = 0.0;              ///< Neutral gas number density [1/m³]
    double mean_thermal_velocity_m_s = 0.0;   ///< Mean thermal velocity √(8kT/πm_n) [m/s]
    double neutral_mass_kg = 0.0;             ///< Neutral molecule mass [kg]
    double CCS_m2 = 0.0;                      ///< Collision cross-section [m²]
    
    // --- Langevin model parameters ---
    double neutral_polarizability_m3 = 0.0;   ///< Neutral polarizability [m³]
    
    // --- Friction model parameters ---
    double reduced_mobility_cm2_Vs = 0.0;     ///< Reduced mobility K₀ [cm²/(V·s)]
};

/**
 * @class DampingForce
 * @brief Computes deterministic collision damping forces
 * 
 * **All models compute:** F = -γ·m·v (friction opposing velocity)
 * 
 * **Damping coefficient γ [1/s] by model:**
 * 
 * 1. **HardSphere** (elastic collisions):
 *    γ = ν_collision = n·σ·v_th·m_reduced/m_ion
 *    where m_reduced = m_n/(m_i+m_n)
 * 
 * 2. **Langevin** (polarization interactions):
 *    γ = ν_Langevin = n·σ_Langevin(v)·v_th·m_reduced/m_ion
 *    where σ_Langevin = π·q·√(α/(4πε₀·m_reduced))/|v|
 * 
 * 3. **Friction** (mobility-based):
 *    γ = q/(K₀·m_ion) where K₀ = reduced mobility
 * 
 * **Usage:**
 * ```cpp
 * // From ion/domain state (automatic model selection)
 * DampingParams params;
 * params.model = DampingModel::Langevin;
 * // Parameters extracted from ForceContext (ion.CCS_m2, domain.env, etc.)
 * auto force = std::make_unique<DampingForce>(params);
 * 
 * // Explicit damping coefficient (overrides model)
 * DampingParams explicit_params;
 * explicit_params.model = DampingModel::Friction;
 * explicit_params.gamma_coefficient = 1e6;  // 1/s
 * auto force = std::make_unique<DampingForce>(explicit_params);
 * ```
 * 
 * @note Random thermal kicks (OU process) are handled by CollisionEngine, NOT here
 * @note Matches legacy defineCollisionForces.cpp behavior (deterministic only)
 */
class DampingForce : public IForce {
public:
    /**
     * @brief Construct damping force from parameters
     * @param params Damping model configuration
     */
    explicit DampingForce(const DampingParams& params);
    
    /**
     * @brief Compute damping force F = -γ·m·v
     * 
     * @param ion Ion state (velocity, mass, charge, CCS, mobility)
     * @param t Current simulation time [s] (unused for damping)
     * @param ctx Force context (domain environment properties)
     * @return Force vector [N]
     */
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override;
    
    /**
     * @brief Get force name
     * @return "Damping(HardSphere)", "Damping(Langevin)", or "Damping(Friction)"
     */
    std::string name() const override;

private:
    DampingParams params_;
    
    /**
     * @brief Calculate damping coefficient γ [1/s] based on collision model
     * 
     * @param ion Ion state (mass, charge, CCS, mobility)
     * @param ctx Force context (domain gas properties)
     * @return γ such that F = -γ·m·v [1/s]
     * 
     * Returns params_.gamma_coefficient if > 0 (explicit override).
     * Otherwise computes from model:
     * - HardSphere: ν_collision from momentum transfer rate
     * - Langevin: ν_Langevin from polarization cross-section
     * - Friction: q/(K₀·m) from reduced mobility
     */
    double calculate_gamma(const IonState& ion, const ForceContext& ctx) const;
};

} // namespace physics
} // namespace ICARION

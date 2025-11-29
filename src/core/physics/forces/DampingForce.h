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
#include <unordered_set>

namespace ICARION {
namespace config {
    struct EnvironmentConfig;
    struct SpeciesDatabase;
}

namespace physics {

/**
 * @enum DampingModel
 * @brief Collision model for damping force calculation
 * 
 * Each model computes F = -γ·m·v with different formulas for γ:
 * - **Friction**: Mobility-based (Mason-Schamp) - **RECOMMENDED**
 * - **HardSphere**: Kinetic theory collisions - **EXPERIMENTAL**
 * - **Langevin**: Polarization interactions - **EXPERIMENTAL** (polar molecules only)
 * 
 * @warning HardSphere and Langevin are experimental! Use Friction for production.
 */
enum class DampingModel {
    None,       ///< No damping
    HardSphere, ///< **EXPERIMENTAL:** Hard-sphere elastic collisions (use with caution)
    Langevin,   ///< **EXPERIMENTAL:** Langevin polarization (only for polar molecules!)
    Friction    ///< **RECOMMENDED:** Mobility-based friction (Mason-Schamp equation)
};

/**
 * @class DampingForce
 * @brief Computes deterministic collision damping forces
 * 
 * **All models compute:** F = -γ·m·v (friction opposing velocity)
 * 
 * **Damping coefficient γ [1/s] by model:**
 * 
 * 1. **Friction** (mobility-based) - **RECOMMENDED**:
 *    γ = q/(K·m_ion) where K = K₀·(n₀/n) is actual mobility
 *    - Validated: H3O+ in N2 gives exact mobility match (356 m/s)
 *    - Based on Mason-Schamp equation (experimentally validated)
 * 
 * 2. **HardSphere** (kinetic theory) - **EXPERIMENTAL**:
 *    γ = n·σ·v_th
 *    - Requires accurate CCS for target gas (not literature values!)
 *    - Works for non-polar gases (N2, He, Ar)
 *    - γ typically 10-20% lower than Friction model
 * 
 * 3. **Langevin** (polarization) - **EXPERIMENTAL**:
 *    γ = n·σ_Langevin(v)·v_th·m_reduced/m_ion
 *    where σ_Langevin = π·q·√(α/(4πε₀·m_reduced))/|v|
 *    - **Only valid for polar molecules!**
 *    - Overpredicts damping for non-polar gases (N2: 76x too strong!)
 *    - Use only for ion-polar molecule interactions
 * 
 * @warning HardSphere and Langevin are experimental. For production simulations,
 *          use **Friction model** which has been extensively validated.
 * 
 * **Usage (SSOT):**
 * ```cpp
 * // Recommended: Friction model
 * const auto& env = domain.environment;
 * auto force = std::make_unique<DampingForce>(env, DampingModel::Friction);
 * ```
 * 
 * @note Random thermal kicks (OU process) are handled by CollisionEngine, NOT here
 * @note Reads gas properties from EnvironmentConfig (SSOT compliance)
 */
class DampingForce : public IForce {
public:
    /**
     * @brief Construct damping force from environment config
     * @param env Environment configuration (SSOT reference)
     * @param model Damping model selection
     */
    DampingForce(
        const ICARION::config::EnvironmentConfig& env,
        DampingModel model,
        const ICARION::config::SpeciesDatabase* species_db = nullptr
    );
    
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
    const ICARION::config::EnvironmentConfig* env_;
    DampingModel model_;
    const ICARION::config::SpeciesDatabase* species_db_ = nullptr;
    mutable std::unordered_set<std::string> warned_missing_sigma_;
    
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

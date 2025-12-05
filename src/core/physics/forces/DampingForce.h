// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file DampingForce.h
 * @brief Deterministic collision damping forces (HardSphere, Langevin, Friction)
 * 
 * Computes velocity-dependent damping forces F = -γ·v for ions in background gas.
 * Thermal kicks are not applied here; stochastic diffusion (OU) is handled in the
 * collision pipeline when enabled.
 */

#pragma once

#include "IForce.h"
#include "ForceContext.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
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
 * - Friction: mobility-based (Mason-Schamp) – recommended
 * - HardSphere: kinetic-theory collisions – experimental
 * - Langevin: polarization interactions – experimental (polar molecules only)
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
 * 1. Friction (mobility-based) – recommended:
 *    γ = q/(K·m_ion) where K = K₀·(n₀/n) is actual mobility
 * 
 * 2. HardSphere (kinetic theory) – experimental:
 *    γ = n·σ·v_th; sensitive to CCS choice and target gas
 * 
 * 3. Langevin (polarization) – experimental, polar molecules only:
 *    γ = n·σ_Langevin(v)·v_th·m_reduced/m_ion
 *    where σ_Langevin = π·q·√(α/(4πε₀·m_reduced))/|v|
 * 
 * @warning HardSphere and Langevin are experimental; prefer Friction for production.
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
    Vec3 compute_soa(const core::IonEnsemble& ensemble, size_t ion_idx, double t,
                     const ForceContext& ctx) const override;
    
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

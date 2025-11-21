// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file DampingForce.h
 * @brief Damping force implementations (friction, Langevin)
 * 
 * Computes velocity-dependent damping forces for ions in background gas:
 * - Friction damping: F = -γ·v (continuous drag)
 * - Langevin damping: F = -γ·v + random thermal kicks (Brownian motion)
 */

#pragma once

#include "IForce.h"
#include "ForceContext.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"

#include <random>

namespace ICARION {
namespace physics {

/**
 * @enum DampingModel
 * @brief Type of damping force model
 */
enum class DampingModel {
    None,      ///< No damping
    Friction,  ///< Deterministic friction: F = -γ·v
    Langevin   ///< Stochastic Langevin: F = -γ·v + ξ(t) (Brownian)
};

/**
 * @brief Parameters for damping force configuration
 */
struct DampingParams {
    DampingModel model = DampingModel::None;
    double damping_coefficient = 0.0;  ///< Damping coefficient γ [kg/s]
    double temperature_K = 300.0;      ///< Gas temperature [K] (for Langevin)
    unsigned int random_seed = 42;     ///< RNG seed for reproducibility
};

/**
 * @class DampingForce
 * @brief Computes velocity-dependent damping forces
 * 
 * **Friction Damping:**
 * - F = -γ·v (deterministic drag)
 * - Energy dissipation without thermal fluctuations
 * - Suitable for high-pressure regime or phenomenological modeling
 * 
 * **Langevin Damping:**
 * - F = -γ·v + ξ(t) (stochastic force)
 * - ξ(t): Gaussian white noise with ⟨ξ⟩ = 0
 * - Variance: σ² = 2·γ·k_B·T·m/Δt (fluctuation-dissipation theorem)
 * - Models Brownian motion in dilute gas
 * 
 * **Usage:**
 * ```cpp
 * // Friction damping (deterministic)
 * DampingParams params;
 * params.model = DampingModel::Friction;
 * params.damping_coefficient = 1e-15;  // kg/s
 * auto force = std::make_unique<DampingForce>(params);
 * 
 * // Langevin damping (stochastic, temperature-dependent)
 * DampingParams langevin;
 * langevin.model = DampingModel::Langevin;
 * langevin.damping_coefficient = 1e-15;
 * langevin.temperature_K = 300.0;
 * langevin.random_seed = 12345;
 * auto force = std::make_unique<DampingForce>(langevin);
 * ```
 * 
 * @note Damping coefficient γ depends on ion-gas collision cross-section
 * @note For Langevin: requires dt from context for proper noise scaling
 */
class DampingForce : public IForce {
public:
    /**
     * @brief Construct damping force from parameters
     * @param params Damping model configuration
     */
    explicit DampingForce(const DampingParams& params);
    
    /**
     * @brief Compute damping force
     * 
     * @param ion Ion state (velocity required!)
     * @param t Current simulation time [s]
     * @param ctx Force context (temperature from environment)
     * @return Force vector [N]
     * 
     * Friction: F = -γ·v
     * Langevin: F = -γ·v + ξ(t), with ξ ~ N(0, σ²)
     */
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override;
    
    /**
     * @brief Get force name
     * @return "Damping(Friction)" or "Damping(Langevin)"
     */
    std::string name() const override;

private:
    /**
     * @brief Compute random thermal force for Langevin dynamics
     * @param ion_mass Ion mass [kg]
     * @param temperature Temperature [K]
     * @param dt Time step [s]
     * @return Random force vector [N]
     */
    Vec3 compute_random_force(double ion_mass, double temperature, double dt) const;
    
    DampingParams params_;
    
    // Mutable RNG for thread-safe random number generation
    mutable std::mt19937 rng_;
    mutable std::normal_distribution<double> normal_dist_;
};

} // namespace physics
} // namespace ICARION

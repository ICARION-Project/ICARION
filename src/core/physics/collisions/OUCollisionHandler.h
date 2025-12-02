// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file OUCollisionHandler.h
 * @brief Ornstein-Uhlenbeck (OU) thermal kick handler
 * 
 * Applies stochastic thermal velocity kicks to maintain Boltzmann distribution.
 * Used as add-on for deterministic damping models (Friction, Langevin, HardSphere).
 * 
 * **Physics:**
 * - Based on Ornstein-Uhlenbeck process (Langevin equation solution)
 * - Ensures correct thermal equilibrium without explicit collision events
 * - Satisfies Fluctuation-Dissipation Theorem
 * 
 * **ONLY VALID WITH DETERMINISTIC DAMPING:**
 * - Friction + OU: Mobility-based damping + thermal noise
 * - Langevin + OU: Velocity-dependent damping + thermal noise
 * - HardSphere + OU: Collision frequency damping + thermal noise
 * 
 * **NOT compatible with:**
 * - EHSS (already has thermal scattering)
 * - HSS (already has thermal scattering)
 * - NoCollisions (no damping to balance)
 * 
 * **SSOT Design:**
 * - Reads temperature directly from `EnvironmentConfig`
 * - Gamma coefficient must match DampingForce setting
 * 
 * @date 2025-11-21
 * @version 1.0
 */

#pragma once

#include "ICollisionHandler.h"

namespace ICARION::physics {

/**
 * @brief OU (Ornstein-Uhlenbeck) thermal kick handler
 * 
 * Applies stochastic velocity kicks to ions undergoing deterministic damping.
 * Ensures correct equilibrium temperature distribution.
 * 
 * **Physical basis:**
 * Langevin equation: m·dv/dt = -γ·m·v + ξ(t)
 * where ξ(t) is white noise with <ξ(t)·ξ(t')> = 2·γ·m·kB·T·δ(t-t')
 * 
 * Solution gives velocity kick:
 * Δv = v_th · √(1 - exp(-2·γ·dt)) · N(0,1)
 * where v_th = √(kB·T/m)
 * 
 * **Use cases:**
 * - Friction model + OU → correct temperature equilibrium
 * - Langevin model + OU → velocity-dependent damping + thermal noise
 * - Hard-sphere model + OU → collision frequency + thermal noise
 * 
 * **SSOT Pattern:**
 * ```cpp
 * // Correct: Gamma must match DampingForce setting
 * double gamma = 1e6;  // [1/s] from DampingForce
 * OUCollisionHandler handler(gamma);
 * 
 * // Handler reads temperature from EnvironmentConfig
 * config::EnvironmentConfig env;
 * env.temperature_K = 300.0;
 * handler.handle_collision(ion, dt, rng, env);  // SSOT!
 * ```
 * 
 * Gamma coefficient should match the damping coefficient used in DampingForce;
 * otherwise equilibrium temperature will be off.
 * 
 * @see DampingForce for deterministic damping models
 * @see apply_ou_velocity_kick() for underlying implementation
 */
class OUCollisionHandler : public ICollisionHandler {
public:
    /**
     * @brief Construct OU handler
     * 
     * @param gamma_coefficient Damping coefficient [1/s] - MUST match DampingForce!
     * @param apply_damping If true, applies full OU (damping + thermalization).
     *                      If false, applies only thermal kicks (for use with DampingForce).
     * 
     *  **IMPORTANT:** This gamma must equal the damping coefficient used in
     * DampingForce, otherwise equilibrium temperature will be wrong!
     * 
     * **Usage:**
     * - With DampingForce: apply_damping=false (only thermal kicks)
     * - Without DampingForce: apply_damping=true (full OU process)
     * 
     * @throws std::invalid_argument if gamma_coefficient <= 0
     */
    explicit OUCollisionHandler(double gamma_coefficient, bool apply_damping = true);
    
    /**
     * @brief Apply OU thermal kick for single timestep
     * 
     * **Algorithm:**
     * 1. Compute thermal velocity: v_th = √(kB·T/m)
     * 2. Compute kick amplitude: A = v_th · √(1 - exp(-2·γ·dt))
     * 3. Sample 3D Gaussian noise: N(0,1)
     * 4. Apply kick: Δv = A · N(0,1)
     * 
     * **SSOT:** Reads temperature directly from `env`:
     * - env.temperature_K → thermal velocity scale
     * 
     * @param[in,out] ion Ion state (velocity modified by thermal kick)
     * @param[in] dt Timestep [s]
     * @param[in,out] rng Random number generator
     * @param[in] env Environment configuration (SSOT!)
     * 
     * @return true (always "collides" - continuous process)
     * 
     * @note Unlike EHSS/HSS, this is a continuous stochastic process,
     *       not a discrete collision event. Return value always true.
     */
    bool handle_collision(
        IonState& ion,
        double dt,
        EhssRng& rng,
        const config::EnvironmentConfig& env
    ) override;
    
    std::string name() const override { return "OU"; }
    
    /**
     * @brief Get damping coefficient
     * @return Gamma coefficient [1/s]
     */
    double gamma() const { return gamma_; }
    
private:
    double gamma_;          ///< Damping coefficient [1/s] - must match DampingForce!
    bool apply_damping_;    ///< If false, only thermal kicks (use with DampingForce)
};

} // namespace ICARION::physics

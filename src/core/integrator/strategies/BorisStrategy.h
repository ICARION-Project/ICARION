// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   ------------------------------------------------
 *   Modular framework for simulating ion trajectories in custom
 *   electric fields and background gas environments.
 *
 *   @file       BorisStrategy.h
 *   @brief      Boris pusher integration strategy for magnetic fields
 *
 *   @details
 *   Symplectic integration method optimized for charged particle motion
 *   in electromagnetic fields. Second-order accurate and preserves energy
 *   in pure magnetic fields.
 *
 *   **Algorithm (Velocity Verlet + Boris rotation):**
 *   1. Half-step velocity advance with electric field: v^- = v^n + (q/m)*E*dt/2
 *   2. Magnetic rotation: v^+ = v^- rotated by q*B*dt/m
 *   3. Half-step velocity advance: v^(n+1) = v^+ + (q/m)*E*dt/2
 *   4. Position update: x^(n+1) = x^n + v^(n+1)*dt
 *
 *   **Properties:**
 *   - Order: 2 (explicit symplectic)
 *   - Timestep: Fixed
 *   - Stability: Excellent for magnetic fields (phase error only, no amplitude growth)
 *   - Cost: 1 force evaluation per step
 *   - Energy: Conserved in pure B-field, bounded error in E+B fields
 *
 *   **When to Use:**
 *   - Strong magnetic fields (ωc*dt ~ 0.1-1.0)
 *   - ICR (Ion Cyclotron Resonance) simulations
 *   - Penning traps, Orbitraps
 *   - Long-term tracking (symplectic → no energy drift)
 *
 *   **When NOT to Use:**
 *   - No magnetic field (use RK4 instead)
 *   - Highly non-uniform E-fields (symplectic advantage lost)
 *   - Adaptive timestep needed (Boris is fixed-step)
 *
 *   **Reference:**
 *   - Boris, J.P. (1970): "Relativistic plasma simulation-optimization of a hybrid code"
 *   - Qin et al. (2013): "Why is Boris algorithm so good?" Physics of Plasmas 20, 084503
 *
 *   @date       2025-11-22
 *   @version    1.0.0
 *   @authors    ICARION Development Team
 *
 * =====================================================================
 */
#pragma once

#include "IIntegrationStrategy.h"

namespace ICARION {
namespace integrator {

/**
 * @brief Boris pusher integration strategy
 * 
 * Symplectic method for charged particle dynamics in electromagnetic fields.
 * Optimized for magnetic field-dominated systems (ICR, Orbitrap, Penning trap).
 * 
 * **Algorithm Details:**
 * ```
 * // Half-step electric kick
 * v_minus = v_n + (q/m) * E * (dt/2)
 * 
 * // Magnetic rotation (Boris rotation)
 * t = (q/m) * B * (dt/2)
 * s = 2*t / (1 + t²)
 * v_prime = v_minus + v_minus × t
 * v_plus = v_minus + v_prime × s
 * 
 * // Half-step electric kick
 * v_n+1 = v_plus + (q/m) * E * (dt/2)
 * 
 * // Position update
 * x_n+1 = x_n + v_n+1 * dt
 * ```
 * 
 * **Stability:**
 * - Cyclotron motion: Stable for all ωc*dt (phase error only)
 * - E×B drift: Accurate to O(dt²)
 * - Energy: Bounded error (symplectic property)
 * 
 * **Optimal Timestep:**
 * - ωc*dt ~ 0.1-0.5 for <1% phase error per orbit
 * - For ωc = 2π*1MHz, use dt ~ 1-5 ns
 * 
 * **SSOT Compliance:**
 * - Uses ForceRegistry (not compute_accelerations())
 * - Uses DomainConfig (not GlobalParams)
 * - Zero-copy references
 */
class BorisStrategy : public IIntegrationStrategy {
public:
    /**
     * @brief Constructor
     * 
     * No configuration needed for fixed-timestep Boris.
     */
    BorisStrategy() = default;
    
    /**
     * @brief Advance ion by one timestep using Boris pusher
     * 
     * @param ion Ion state (updated in-place)
     * @param t Current time [s]
     * @param dt Timestep [s]
     * @param force_registry Force computation engine
     * @param domain Domain configuration (SSOT!)
     * @param all_ions All ions (for space charge)
     * 
     * **Implementation Details:**
     * 1. Compute electric and magnetic forces at current position
     * 2. Half-step electric acceleration
     * 3. Boris rotation in magnetic field
     * 4. Half-step electric acceleration
     * 5. Position update with new velocity
     * 
     * **Force Decomposition:**
     * - Electric force: F_E = q*E → contributes to v^- and v^+
     * - Magnetic force: F_B = q*v×B → rotation operator
     * - Other forces (drag, collisions): treated as electric-like
     * 
     * **Numerical Stability:**
     * - Unconditionally stable for pure B-field
     * - Stable for E-field if dt < 2/ω_plasma
     * - Typical safety: dt < 0.1 * T_cyclotron
     * 
     * **Thread Safety:**
     * - Read-only access to force_registry, domain, all_ions
     * - Modifies only `ion` parameter
     * - Can parallelize over ions (no shared state)
     */
    void step(
        IonState& ion,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry,
        const config::DomainConfig& domain,
        const std::vector<IonState>& all_ions
    ) override;
    
    std::string name() const override { return "Boris"; }
    bool is_adaptive() const override { return false; }
};

} // namespace integrator
} // namespace ICARION

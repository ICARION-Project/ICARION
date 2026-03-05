// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
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
     * @brief SoA-aware Boris step (avoids AoS reconstruction)
     */
    void step(
        core::IonEnsemble& ensemble,
        size_t ion_idx,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry
    ) override;
    
    std::string name() const override { return "Boris"; }
    bool is_adaptive() const override { return false; }
};

} // namespace integrator
} // namespace ICARION

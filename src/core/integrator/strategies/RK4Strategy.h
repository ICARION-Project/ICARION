// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IIntegrationStrategy.h"

namespace ICARION {
namespace integrator {

/**
 * @brief 4th-order Runge-Kutta integration strategy
 * 
 * Standard RK4 method for solving ODEs: dy/dt = f(t, y)
 * 
 * **Use Cases:**
 * - General-purpose integration (most common choice)
 * - Smooth, deterministic forces
 * - Fixed timestep requirements
 * 
 * **Not Recommended For:**
 * - Stiff systems (use implicit methods)
 * - Highly oscillatory systems (use symplectic methods like Boris)
 * - Adaptive timestep requirements (use RK45Strategy)
 * 
 * **Performance:**
 * - 4 force evaluations per step
 * - ~2x slower than Euler, ~2x faster than RK45
 * - Memory: O(1) per ion (no state history)
 * 
 * **SSOT Compliance:**
 * - Uses ForceRegistry (not compute_accelerations())
 * - Uses DomainConfig (not GlobalParams)
 * - Zero-copy references
 */
class RK4Strategy : public IIntegrationStrategy {
public:
    /**
     * @brief Constructor
     * 
     * No configuration needed for fixed-timestep RK4.
     */
    RK4Strategy() = default;
    
    /**
     * @brief Advance ion by one timestep using RK4
     * 
     * @param ion Ion state (updated in-place)
     * @param t Current time [s]
     * @param dt Timestep [s]
     * @param force_registry Force computation engine
     * @param domain Domain configuration (SSOT!)
     * @param all_ions All ions (for space charge)
     * 
     * **Implementation Details:**
     * - Computes 4 intermediate stages (k1, k2, k3, k4)
     * - Each stage requires force evaluation at different (t, y)
     * - Final update uses weighted average of stages
     * 
     * **Numerical Stability:**
     * - Stable for dt < 2.78/|λ_max| (linear systems)
     * - For nonlinear systems, empirical testing required
     * - Typical safe choice: dt ~ T_oscillation / 20
     * 
     * **Thread Safety:**
     * - Read-only access to force_registry, all_ions
     * - Modifies only `ion` parameter
     * - Can parallelize over multiple ions (no shared state)
     */
    void step(
        IonState& ion,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry,
        const std::vector<IonState>& all_ions
    ) override;
    
    /**
     * @brief RK4 integration using SoA (Phase 3B - optimized)
     * 
     * Direct array access for cache efficiency.
     * No IonState construction overhead.
     */
    void step_soa(
        core::IonEnsemble& ensemble,
        size_t ion_idx,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry
    ) override;
    
    std::string name() const override { return "RK4"; }
    bool is_adaptive() const override { return false; }
};

} // namespace integrator
} // namespace ICARION

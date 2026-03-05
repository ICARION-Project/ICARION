// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
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
     * @brief RK4 integration using SoA layout
     * 
     * Direct array access for cache efficiency.
     * No IonState construction overhead.
     */
    void step(
        core::IonEnsemble& ensemble,
        size_t ion_idx,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry
    ) override;

    bool supports_batch() const override { return true; }

    bool step_batch(
        core::IonEnsemble& ensemble,
        double t,
        double dt,
        const std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
        const std::vector<int>& domain_indices
    ) override;

    void set_parallel_enabled(bool enabled) override { parallel_enabled_ = enabled; }
    
    std::string name() const override { return "RK4"; }
    bool is_adaptive() const override { return false; }

private:
    bool parallel_enabled_ = false;
};

} // namespace integrator
} // namespace ICARION

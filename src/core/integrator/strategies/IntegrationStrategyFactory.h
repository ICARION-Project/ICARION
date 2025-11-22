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
 *   @file       IntegrationStrategyFactory.h
 *   @brief      Factory for creating integration strategies
 *
 *   @details
 *   Factory pattern for strategy selection.
 *   Supports: RK4, RK45 (future), Boris (future)
 *
 *   @date       2025-11-22
 *   @version    1.0.0
 *   @authors    ICARION Development Team
 *
 * =====================================================================
 */
#pragma once

#include <memory>
#include <string>
#include <stdexcept>

#include "IIntegrationStrategy.h"
#include "RK4Strategy.h"

namespace ICARION {
namespace integrator {

/**
 * @brief Factory for creating integration strategies
 * 
 * **Supported Strategies:**
 * - "RK4" (default): 4th-order Runge-Kutta (fixed timestep)
 * - "RK45" (Phase 4B): Dormand-Prince (adaptive timestep)
 * - "Boris" (Phase 4B): Boris pusher (magnetic fields)
 * 
 * **Example Usage:**
 * ```cpp
 * // Create from string identifier
 * auto strategy = IntegrationStrategyFactory::create("RK4");
 * 
 * // Use in simulation
 * strategy->step(ion, t, dt, force_registry, domain, all_ions);
 * ```
 */
class IntegrationStrategyFactory {
public:
    /**
     * @brief Create integration strategy from identifier
     * 
     * @param strategy_name Strategy identifier ("RK4", "RK45", "Boris")
     * @return Unique pointer to strategy instance
     * @throws std::invalid_argument if strategy_name unknown
     * 
     * **SSOT Compliance:**
     * - Reads strategy_name from config.simulation.integrator
     * - No hardcoded defaults (caller must specify)
     * 
     * **Supported Values:**
     * - "RK4": Classic 4th-order Runge-Kutta (Phase 4A)
     * - "RK45": Adaptive Dormand-Prince (Phase 4B, not yet implemented)
     * - "Boris": Boris pusher for magnetic fields (Phase 4B, not yet implemented)
     * 
     * **Thread Safety:**
     * - Stateless factory (thread-safe)
     * - Returns new instance per call (no sharing)
     */
    static std::unique_ptr<IIntegrationStrategy> create(const std::string& strategy_name) {
        if (strategy_name == "RK4") {
            return std::make_unique<RK4Strategy>();
        }
        else if (strategy_name == "RK45") {
            // Phase 4B: Not yet implemented
            throw std::invalid_argument(
                "IntegrationStrategyFactory: RK45 strategy not yet implemented (Phase 4B). "
                "Use 'RK4' for now."
            );
        }
        else if (strategy_name == "Boris") {
            // Phase 4B: Not yet implemented
            throw std::invalid_argument(
                "IntegrationStrategyFactory: Boris strategy not yet implemented (Phase 4B). "
                "Use 'RK4' for now."
            );
        }
        else {
            throw std::invalid_argument(
                "IntegrationStrategyFactory: Unknown strategy '" + strategy_name + "'. "
                "Supported: 'RK4', 'RK45' (Phase 4B), 'Boris' (Phase 4B)."
            );
        }
    }
    
    /**
     * @brief Get list of supported strategies
     * 
     * @return Vector of strategy identifiers
     * 
     * Useful for validation and CLI help text.
     */
    static std::vector<std::string> supported_strategies() {
        return {
            "RK4"     // Phase 4A (implemented)
            // "RK45",   // Phase 4B (not yet implemented)
            // "Boris"   // Phase 4B (not yet implemented)
        };
    }
};

} // namespace integrator
} // namespace ICARION

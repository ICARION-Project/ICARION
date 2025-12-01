// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <memory>
#include <string>
#include <stdexcept>

#include "IIntegrationStrategy.h"
#include "RK4Strategy.h"
#include "RK45Strategy.h"
#include "BorisStrategy.h"

namespace ICARION {
namespace integrator {

/**
 * @brief Factory for creating integration strategies
 * 
 * **Supported Strategies:**
 * - "RK4" (default): 4th-order Runge-Kutta (fixed timestep)
 * - "RK45": Dormand-Prince (adaptive timestep)
 * - "Boris": Boris pusher (magnetic fields)
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
     * - "RK4": Classic 4th-order Runge-Kutta
     * - "RK45": Adaptive Dormand-Prince
     * - "Boris": Boris pusher for magnetic fields
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
            // Adaptive Dormand-Prince with default config
            return std::make_unique<RK45Strategy>();
        }
        else if (strategy_name == "Boris") {
            // Boris pusher for magnetic fields
            return std::make_unique<BorisStrategy>();
        }
        else {
            throw std::invalid_argument(
                "IntegrationStrategyFactory: Unknown strategy '" + strategy_name + "'. "
                "Supported: 'RK4', 'RK45', 'Boris'."
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
            "RK4",     // Classic 4th-order Runge-Kutta
            "RK45",    // Adaptive Dormand-Prince
            "Boris"    // Boris pusher
        };
    }
};

} // namespace integrator
} // namespace ICARION

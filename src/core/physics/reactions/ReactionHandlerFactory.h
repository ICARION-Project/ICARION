// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IReactionHandler.h"
#include "core/config/types/PhysicsConfig.h"
#include <memory>

namespace ICARION {
namespace physics {

/**
 * @brief Factory for creating reaction handlers
 * 
 * **SSOT Design:** Factory reads enable_reactions directly from PhysicsConfig.
 * No intermediate conversion or parameter copies.
 * 
 * **Design Pattern:** Factory Pattern
 * - Encapsulates handler creation logic
 * - Easy to extend with new reaction models
 * - Consistent with CollisionHandlerFactory
 * 
 * **Example Usage:**
 * @code
 * // In integrate_trajectory():
 * auto reaction_handler = ReactionHandlerFactory::create(
 *     full_config.physics,  // Direct config reference (SSOT!)
 *     enable_logging
 * );
 * 
 * if (reaction_handler) {
 *     logger->log("Reaction handler created: " + reaction_handler->name());
 * }
 * @endcode
 * 
 * @see CollisionHandlerFactory (similar pattern)
 */
class ReactionHandlerFactory {
public:
    /**
     * @brief Create reaction handler from config
     * 
     * @param config Physics configuration (contains enable_reactions flag)
     * @param enable_logging Enable debug logging (optional)
     * 
     * @return Reaction handler instance (nullptr if reactions disabled)
     * 
     * **SSOT:** Reads `config.enable_reactions` directly (no conversion!).
     * 
 * **Behavior:**
 * - If `enable_reactions == false`: returns nullptr
 * - If `enable_reactions == true`: returns StochasticReactionHandler
     */
    static std::unique_ptr<IReactionHandler> create(
        const config::PhysicsConfig& config,
        bool enable_logging = false
    );
};

} // namespace physics
} // namespace ICARION

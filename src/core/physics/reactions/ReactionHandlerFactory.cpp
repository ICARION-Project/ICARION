// ReactionHandlerFactory.cpp
// Implementation of reaction handler factory
//
// Created: 2025-11-22 (Phase 3 Refactor)

#include "ReactionHandlerFactory.h"
#include "StochasticReactionHandler.h"
#include "NoReactionHandler.h"

namespace ICARION {
namespace physics {

std::unique_ptr<IReactionHandler> ReactionHandlerFactory::create(
    const config::PhysicsConfig& config,
    bool enable_logging
) {
    // SSOT: Read enable_reactions directly from config
    if (config.enable_reactions) {
        return std::make_unique<StochasticReactionHandler>(enable_logging);
    } else {
        // Null Object Pattern: Return no-op handler (avoids null checks in caller)
        return std::make_unique<NoReactionHandler>();
    }
}

} // namespace physics
} // namespace ICARION

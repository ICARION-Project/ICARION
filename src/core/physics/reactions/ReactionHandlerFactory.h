// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
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
     * @param enable_gpu Request GPU backend (requires CUDA build)
     * @param gpu_seed RNG seed for GPU helper
     * @param gpu_threshold Minimum ensemble size for GPU dispatch
     *
     * @return Reaction handler instance (NoReactionHandler if disabled)
     *
     * **SSOT:** Reads `config.enable_reactions` directly (no conversion!).
     */
    static std::unique_ptr<IReactionHandler> create(
        const config::PhysicsConfig& config,
        bool enable_logging = false,
        bool enable_gpu = false,
        unsigned long long gpu_seed = 0,
        size_t gpu_threshold = 2000
    );
};

} // namespace physics
} // namespace ICARION

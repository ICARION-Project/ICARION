// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once

#include "IReactionHandler.h"

namespace ICARION {
namespace physics {

/**
 * @brief Null object pattern for disabled reactions
 * 
 * Used when config.physics.enable_reactions = false.
 * Does nothing (no-op), avoiding null pointer checks in integration loop.
 * 
 * **Design Pattern:** Null Object Pattern
 * - Provides valid interface without functionality
 * - Eliminates need for `if (handler != nullptr)` checks
 * - Zero performance overhead (compiler optimizes away empty function)
 * 
 * **Example Usage:**
 * @code
 * auto handler = (config.physics.enable_reactions)
 *     ? std::make_unique<StochasticReactionHandler>()
 *     : std::make_unique<NoReactionHandler>();
 * 
 * // Safe to call without null check:
 * handler->handle_reaction(...);  // No-op if NoReactionHandler
 * @endcode
 */
class NoReactionHandler : public IReactionHandler {
public:
    bool handle_reaction(
        core::IonReactionData&,
        double,
        PhysicsRng&,
        const config::ReactionDatabase&,
        const config::SpeciesDatabase&,
        const config::EnvironmentConfig&
    ) override {
        // No-op: reactions disabled
        return false;
    }
    
    std::string name() const override { return "None"; }
};

} // namespace physics
} // namespace ICARION

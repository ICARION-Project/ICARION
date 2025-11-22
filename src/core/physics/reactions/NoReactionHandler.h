// NoReactionHandler.h
// Null object pattern for disabled reactions
//
// Created: 2025-11-22 (Phase 3 Refactor)

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
        IonState& ion,
        double dt,
        EhssRng& rng,
        const config::ReactionDatabase& reaction_db,
        const config::SpeciesDatabase& species_db,
        const config::EnvironmentConfig& env
    ) override {
        // No-op: reactions disabled
        return false;
    }
    
    std::string name() const override { return "None"; }
};

} // namespace physics
} // namespace ICARION

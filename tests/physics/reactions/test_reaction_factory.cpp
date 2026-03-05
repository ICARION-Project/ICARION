// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
// test_reaction_factory.cpp
// Unit tests for ReactionHandlerFactory
//
// Tests:
// - Factory creation from PhysicsConfig
// - SSOT compliance (direct config reference)
// - Correct handler type returned
// - Null handling when reactions disabled
//
// Created: 2025-11-22 (Phase 3 Refactor)

#include <catch2/catch_test_macros.hpp>
#include "core/physics/reactions/ReactionHandlerFactory.h"
#include "core/physics/reactions/StochasticReactionHandler.h"
#include "core/physics/reactions/IReactionHandler.h"
#include "core/config/types/PhysicsConfig.h"

using namespace ICARION;
using namespace ICARION::physics;
using namespace ICARION::config;

TEST_CASE("ReactionHandlerFactory: Create StochasticReactionHandler", "[reaction][factory]") {
    PhysicsConfig config;
    config.enable_reactions = true;  // Enable reactions
    config.collision_model = CollisionModel::NoCollisions;
    config.enable_space_charge = false;
    
    // SSOT: Pass config directly
    auto handler = ReactionHandlerFactory::create(config, false);
    
    // Verify: Handler created
    REQUIRE(handler != nullptr);
    REQUIRE(handler->name() == "Stochastic");
}

TEST_CASE("ReactionHandlerFactory: Reactions disabled", "[reaction][factory]") {
    PhysicsConfig config;
    config.enable_reactions = false;  // ✅ Disable reactions
    config.collision_model = CollisionModel::NoCollisions;
    config.enable_space_charge = false;
    
    // ✅ SSOT: Pass config directly
    auto handler = ReactionHandlerFactory::create(config, false);
    
    // Verify: NoReactionHandler returned (Null Object Pattern)
    REQUIRE(handler != nullptr);
    REQUIRE(handler->name() == "None");
}

TEST_CASE("ReactionHandlerFactory: Enable logging", "[reaction][factory]") {
    PhysicsConfig config;
    config.enable_reactions = true;
    
    // Create with logging enabled
    auto handler = ReactionHandlerFactory::create(config, true);
    
    REQUIRE(handler != nullptr);
    REQUIRE(handler->name() == "Stochastic");
    // Note: Logging behavior tested in integration tests
}

TEST_CASE("ReactionHandlerFactory: SSOT compliance", "[reaction][factory][ssot]") {
    PhysicsConfig config;
    config.enable_reactions = true;
    
    // SSOT: Factory reads enable_reactions directly from PhysicsConfig
    // No intermediate conversion, no parameter copies
    auto handler = ReactionHandlerFactory::create(config, false);
    
    REQUIRE(handler != nullptr);
    
    // Change config (handler creation is decoupled)
    config.enable_reactions = false;
    
    // Handler still valid (created from previous config state)
    REQUIRE(handler->name() == "Stochastic");
}

TEST_CASE("ReactionHandlerFactory: GPU request without CUDA falls back", "[reaction][factory][gpu]") {
    PhysicsConfig config;
    config.enable_reactions = true;

    auto handler = ReactionHandlerFactory::create(
        config,
        false,
        true,   // enable_gpu
        42,
        16
    );

    REQUIRE(handler != nullptr);
#ifdef ICARION_USE_GPU
    CHECK(handler->supports_batch());
#else
    CHECK_FALSE(handler->supports_batch());
#endif
}

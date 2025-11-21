// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file test_collision_factory.cpp
 * @brief Unit tests for CollisionHandlerFactory
 * 
 * Tests factory creation logic, validation, and SSOT compliance.
 */

#include "core/physics/collisions/CollisionHandlerFactory.h"
#include "core/physics/collisions/EHSSCollisionHandler.h"
#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/physics/collisions/OUCollisionHandler.h"
#include <catch2/catch_test_macros.hpp>
#include <memory>

using namespace ICARION::physics;
using namespace ICARION::config;

// =====================================================================
// Test Fixture Helper
// =====================================================================

struct CollisionFactoryTestFixture {
    CollisionFactoryTestFixture() {
        // Create minimal geometry map for EHSS tests
        geometry_map_["test_species"] = GeometryData{
            {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}},  // positions
            {1.0, 1.0}                            // radii
        };
    }

    PhysicsConfig config_;
    GeometryMap geometry_map_;
    const double gamma_test_ = 1e5;  // Test gamma coefficient
};

// =====================================================================
// Test Cases
// =====================================================================

TEST_CASE("CollisionHandlerFactory: EHSS handler creation", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::EHSS;
    
    auto handler = CollisionHandlerFactory::create(
        fixture.config_,
        &fixture.geometry_map_,
        0.0,
        false
    );
    
    REQUIRE(handler != nullptr);
    
    // Verify it's actually an EHSS handler (dynamic cast)
    auto* ehss = dynamic_cast<EHSSCollisionHandler*>(handler.get());
    REQUIRE(ehss != nullptr);
}

TEST_CASE("CollisionHandlerFactory: EHSS requires geometry", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::EHSS;
    
    REQUIRE_THROWS_AS(
        CollisionHandlerFactory::create(fixture.config_, nullptr, 0.0, false),
        std::invalid_argument
    );
}

TEST_CASE("CollisionHandlerFactory: HSS handler creation", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::HSS;
    
    auto handler = CollisionHandlerFactory::create(
        fixture.config_,
        nullptr,  // HSS doesn't need geometry
        0.0,
        false
    );
    
    REQUIRE(handler != nullptr);
    
    auto* hss = dynamic_cast<HSSCollisionHandler*>(handler.get());
    REQUIRE(hss != nullptr);
}

TEST_CASE("CollisionHandlerFactory: OU handler with Friction", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::Friction;
    fixture.config_.enable_ou_thermalization = true;
    
    auto handler = CollisionHandlerFactory::create(
        fixture.config_,
        nullptr,
        fixture.gamma_test_,
        false
    );
    
    REQUIRE(handler != nullptr);
    
    auto* ou = dynamic_cast<OUCollisionHandler*>(handler.get());
    REQUIRE(ou != nullptr);
}

TEST_CASE("CollisionHandlerFactory: OU requires positive gamma", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::Langevin;
    fixture.config_.enable_ou_thermalization = true;
    
    // gamma <= 0 should throw
    REQUIRE_THROWS_AS(
        CollisionHandlerFactory::create(fixture.config_, nullptr, 0.0, false),
        std::invalid_argument
    );
    
    REQUIRE_THROWS_AS(
        CollisionHandlerFactory::create(fixture.config_, nullptr, -1.0, false),
        std::invalid_argument
    );
}

TEST_CASE("CollisionHandlerFactory: Deterministic without OU returns null", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::Friction;
    fixture.config_.enable_ou_thermalization = false;
    
    auto handler = CollisionHandlerFactory::create(
        fixture.config_,
        nullptr,
        0.0,
        false
    );
    
    REQUIRE(handler == nullptr);
}

TEST_CASE("CollisionHandlerFactory: HSD without OU returns null", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::HSD;
    fixture.config_.enable_ou_thermalization = false;
    
    auto handler = CollisionHandlerFactory::create(
        fixture.config_,
        nullptr,
        0.0,
        false
    );
    
    REQUIRE(handler == nullptr);
}

TEST_CASE("CollisionHandlerFactory: NoCollisions returns null", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::NoCollisions;
    
    auto handler = CollisionHandlerFactory::create(
        fixture.config_,
        nullptr,
        0.0,
        false
    );
    
    REQUIRE(handler == nullptr);
}

TEST_CASE("CollisionHandlerFactory: Unknown model throws", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::UnknownCollisionModel;
    
    REQUIRE_THROWS_AS(
        CollisionHandlerFactory::create(fixture.config_, nullptr, 0.0, false),
        std::invalid_argument
    );
}

TEST_CASE("CollisionHandlerFactory: OU with HSD is valid", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::HSD;
    fixture.config_.enable_ou_thermalization = true;
    
    auto handler = CollisionHandlerFactory::create(
        fixture.config_,
        nullptr,
        fixture.gamma_test_,
        false
    );
    
    REQUIRE(handler != nullptr);
    
    auto* ou = dynamic_cast<OUCollisionHandler*>(handler.get());
    REQUIRE(ou != nullptr);
}

TEST_CASE("CollisionHandlerFactory: Logging enabled works", "[collision][factory]") {
    CollisionFactoryTestFixture fixture;
    fixture.config_.collision_model = ICARION::config::CollisionModel::HSS;
    
    // Should not throw with logging enabled
    auto handler = CollisionHandlerFactory::create(
        fixture.config_,
        nullptr,
        0.0,
        true  // enable_logging
    );
    
    REQUIRE(handler != nullptr);
}

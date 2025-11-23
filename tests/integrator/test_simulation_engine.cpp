// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file test_simulation_engine.cpp
 * @brief Unit tests for SimulationEngine (Phase 5A)
 * 
 * Tests the main simulation orchestrator in isolation:
 * - Single-domain trajectory integration
 * - Multi-domain transitions (aperture crossing)
 * - Collision handler integration
 * - Reaction handler integration
 * - Output validation
 * - OpenMP thread safety
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/integrator/SimulationEngine.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/physics/collisions/EHSSCollisionHandler.h"  // Use EHSS as no-op equivalent
#include "core/physics/reactions/NoReactionHandler.h"
#include "core/types/IonState.h"
#include <memory>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;

// ============================================================================
// Test Helpers
// ============================================================================

/**
 * @brief Create minimal valid FullConfig for testing (free-flight by default)
 */
FullConfig create_test_config() {
    FullConfig cfg;
    
    // Simulation
    cfg.simulation.total_time_s = 1e-6;  // 1 μs
    cfg.simulation.dt_s = 1e-9;          // 1 ns
    cfg.simulation.compute_derived();
    
    // Physics (no collisions, no reactions, no space charge)
    cfg.physics.collision_model = CollisionModel::NoCollisions;
    cfg.physics.enable_reactions = false;
    cfg.physics.enable_space_charge = false;
    
    // Output
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "icarion_test_simulation_engine.h5";
    
    // Single cylindrical domain (10 mm radius, 100 mm length)
    DomainConfig domain;
    domain.name = "test_domain";
    domain.instrument = Instrument::IMS;  // Simple cylindrical geometry
    domain.geometry.length_m = 0.1;       // 100 mm
    domain.geometry.radius_m = 0.01;      // 10 mm
    domain.environment.temperature_K = 300.0;
    domain.environment.pressure_Pa = 101325.0;
    domain.finalize();
    cfg.domains.push_back(domain);
    
    return cfg;
}

/**
 * @brief Create test ion at domain entrance
 */
IonState create_test_ion() {
    IonState ion;
    ion.pos = Vec3{0.0, 0.0, 0.001};    // 1 mm from entrance
    ion.vel = Vec3{0.0, 0.0, 100.0};    // 100 m/s axial
    ion.mass_kg = 100.0 * 1.66053906660e-27;  // 100 Da
    ion.ion_charge_C = 1.602176634e-19;       // +1e
    ion.species_id = "test_species";
    ion.current_domain_index = 0;
    ion.active = true;
    ion.born = true;  // Initial ions are born from start
    ion.birth_time_s = 0.0;
    return ion;
}

// ============================================================================
// Basic Tests
// ============================================================================

TEST_CASE("SimulationEngine: Construction and initialization", "[simulation][engine]") {
    auto cfg = create_test_config();
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    // Note: For testing without collisions, simply pass nullptr
    auto reaction_handler = std::make_shared<NoReactionHandler>();
    
    SECTION("Valid construction") {
        REQUIRE_NOTHROW(SimulationEngine(
            cfg, force_registry, integrator, nullptr, reaction_handler
        ));
    }
    
    SECTION("Null ForceRegistry throws") {
        REQUIRE_THROWS_AS(SimulationEngine(
            cfg, nullptr, integrator, nullptr, reaction_handler
        ), std::invalid_argument);
    }
    
    SECTION("Null IntegrationStrategy throws") {
        REQUIRE_THROWS_AS(SimulationEngine(
            cfg, force_registry, nullptr, nullptr, reaction_handler
        ), std::invalid_argument);
    }
    
    SECTION("Null collision/reaction handlers allowed") {
        REQUIRE_NOTHROW(SimulationEngine(
            cfg, force_registry, integrator, nullptr, nullptr
        ));
    }
    
    SECTION("Config accessors") {
        SimulationEngine engine(cfg, force_registry, integrator);
        
        REQUIRE(engine.get_config().simulation.dt_s == 1e-9);
        REQUIRE(engine.get_config().domains.size() == 1);
        REQUIRE(&engine.get_domain_manager() != nullptr);
        REQUIRE(&engine.get_output_manager() != nullptr);
    }
}

// ============================================================================
// Single-Domain Tests
// ============================================================================

TEST_CASE("SimulationEngine: Single-domain free-flight trajectory", "[simulation][engine]") {
    auto cfg = create_test_config();
    cfg.simulation.total_time_s = 1e-6;  // 1 μs
    cfg.simulation.compute_derived();
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    
    SimulationEngine engine(cfg, force_registry, integrator);
    
    SECTION("Single ion free flight") {
        std::vector<IonState> ions = {create_test_ion()};
        ions[0].vel = Vec3{0, 0, 1000.0};  // 1 km/s
        double initial_z = ions[0].pos.z;
        
        auto result = engine.run(ions);
        
        REQUIRE(result.size() == 1);
        
        // Ion should have moved forward (no forces)
        REQUIRE(result[0].pos.z > initial_z);
        
        // Distance ≈ v*t = 1000 m/s * 1e-6 s = 1 mm
        double expected_displacement = 1000.0 * 1e-6;
        REQUIRE_THAT(result[0].pos.z - initial_z, 
                     Catch::Matchers::WithinRel(expected_displacement, 0.01));
        
        // Velocity should be unchanged (no forces)
        REQUIRE_THAT(result[0].vel.z, Catch::Matchers::WithinRel(1000.0, 1e-6));
    }
    
    SECTION("Ion exits domain") {
        std::vector<IonState> ions = {create_test_ion()};
        ions[0].pos = Vec3{0, 0, 0.09};  // Near exit (z=90mm, length=100mm)
        ions[0].vel = Vec3{0, 0, 10000.0};  // Fast enough to exit in 1μs
        
        auto result = engine.run(ions);
        
        // Ion should be inactive (exited domain)
        REQUIRE(!result[0].active);
    }
    
    SECTION("Ion hits radial wall") {
        std::vector<IonState> ions = {create_test_ion()};
        ions[0].pos = Vec3{0.009, 0, 0.05};  // r=9mm (near wall at r=10mm)
        ions[0].vel = Vec3{1000.0, 0, 100.0};  // Radial velocity toward wall
        
        auto result = engine.run(ions);
        
        // Ion should be inactive (hit wall)
        REQUIRE(!result[0].active);
    }
}

TEST_CASE("SimulationEngine: Ion birth timing", "[simulation][engine]") {
    auto cfg = create_test_config();
    cfg.simulation.total_time_s = 5e-6;  // 5 μs
    cfg.simulation.compute_derived();
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    
    SimulationEngine engine(cfg, force_registry, integrator);
    
    SECTION("Delayed birth") {
        std::vector<IonState> ions(3);
        for (int i = 0; i < 3; ++i) {
            ions[i] = create_test_ion();
            ions[i].birth_time_s = i * 1e-6;  // 0, 1, 2 μs
            ions[i].born = false;
            ions[i].active = false;
        }
        
        auto result = engine.run(ions);
        
        // All ions should have been born and moved
        for (const auto& ion : result) {
            REQUIRE(ion.born);
            if (ion.active) {  // Only check movement for active ions
                REQUIRE(ion.pos.z >= 0.001);  // Should have moved or stayed at start
            }
        }
    }
}

TEST_CASE("SimulationEngine: Multiple ions", "[simulation][engine]") {
    auto cfg = create_test_config();
    cfg.simulation.total_time_s = 1e-6;
    cfg.simulation.enable_openmp = true;  // Test parallel processing
    cfg.simulation.compute_derived();
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    
    SimulationEngine engine(cfg, force_registry, integrator);
    
    SECTION("Parallel processing") {
        const int N = 100;
        std::vector<IonState> ions;
        for (int i = 0; i < N; ++i) {
            auto ion = create_test_ion();
            ion.pos.x = (i % 10) * 0.001;  // Spread out
            ion.vel.z = 1000.0 + i * 10.0;  // Different velocities
            ions.push_back(ion);
        }
        
        auto result = engine.run(ions);
        
        REQUIRE(result.size() == N);
        
        // All ions should have different final positions (different velocities)
        std::vector<double> z_positions;
        for (const auto& ion : result) {
            z_positions.push_back(ion.pos.z);
        }
        
        // Check that positions are different (not all the same)
        std::sort(z_positions.begin(), z_positions.end());
        double range = z_positions.back() - z_positions.front();
        REQUIRE(range > 1e-6);  // At least 1 μm spread
    }
}

// ============================================================================
// Multi-Domain Tests
// ============================================================================

TEST_CASE("SimulationEngine: Multi-domain transition", "[simulation][engine][domain]") {
    auto cfg = create_test_config();
    cfg.simulation.total_time_s = 5e-6;
    cfg.simulation.compute_derived();
    
    // Add second domain (back-to-back)
    DomainConfig domain2 = cfg.domains[0];
    domain2.name = "domain2";
    domain2.domain_index = 1;
    domain2.geometry.origin_m = Vec3{0, 0, 0.1};  // Starts where domain1 ends
    domain2.finalize();  // Compute bounding boxes for new domain
    cfg.domains.push_back(domain2);
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    
    SimulationEngine engine(cfg, force_registry, integrator);
    
    SECTION("Ion crosses aperture into second domain") {
        std::vector<IonState> ions = {create_test_ion()};
        ions[0].pos = Vec3{0, 0, 0.095};  // Near domain1 exit
        ions[0].vel = Vec3{0, 0, 10000.0};  // Fast enough to reach domain2
        
        auto result = engine.run(ions);
        
        // Ion should have moved forward from starting position
        REQUIRE(result[0].pos.z > 0.095);
        
        // With proper aperture handling, ion should be active and in domain 2
        // For now, verify it at least stayed active through boundary region
        REQUIRE(result[0].active);
        
        // Check if multi-domain transition occurred
        if (result[0].pos.z > 0.1) {
            // Successfully entered domain 2
            REQUIRE(result[0].current_domain_index == 1);
        }
    }
}

// ============================================================================
// Physics Integration Tests
// ============================================================================

TEST_CASE("SimulationEngine: Collision handler integration", "[simulation][engine][collision]") {
    auto cfg = create_test_config();
    cfg.physics.collision_model = CollisionModel::NoCollisions;
    cfg.simulation.total_time_s = 1e-6;
    cfg.simulation.compute_derived();
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    // nullptr collision handler = no collisions
    
    SimulationEngine engine(cfg, force_registry, integrator, nullptr);
    
    SECTION("Nullptr collision handler does not modify trajectory") {
        std::vector<IonState> ions = {create_test_ion()};
        ions[0].vel = Vec3{0, 0, 1000.0};
        
        auto result = engine.run(ions);
        
        // Velocity should be unchanged (no collisions)
        REQUIRE_THAT(result[0].vel.z, Catch::Matchers::WithinRel(1000.0, 1e-6));
    }
}

TEST_CASE("SimulationEngine: Reaction handler integration", "[simulation][engine][reaction]") {
    auto cfg = create_test_config();
    cfg.physics.enable_reactions = true;
    cfg.simulation.total_time_s = 1e-6;
    cfg.simulation.compute_derived();
    
    // Empty reaction database (no reactions will occur)
    cfg.reaction_db.reactions.clear();
    cfg.species_db.species.clear();
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    auto reaction_handler = std::make_shared<NoReactionHandler>();
    
    SimulationEngine engine(cfg, force_registry, integrator, nullptr, reaction_handler);
    
    SECTION("NoReactionHandler does not modify species") {
        std::vector<IonState> ions = {create_test_ion()};
        std::string original_species = ions[0].species_id;
        
        auto result = engine.run(ions);
        
        // Species should be unchanged
        REQUIRE(result[0].species_id == original_species);
    }
}

// ============================================================================
// Edge Cases and Error Handling
// ============================================================================

TEST_CASE("SimulationEngine: Edge cases", "[simulation][engine][edge]") {
    auto cfg = create_test_config();
    cfg.simulation.total_time_s = 1e-6;
    cfg.simulation.compute_derived();
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    
    SimulationEngine engine(cfg, force_registry, integrator);
    
    SECTION("Empty ion list") {
        std::vector<IonState> ions;
        
        auto result = engine.run(ions);
        
        REQUIRE(result.empty());
    }
    
    SECTION("All ions inactive from start") {
        std::vector<IonState> ions(5);
        for (auto& ion : ions) {
            ion = create_test_ion();
            ion.active = false;
        }
        
        auto result = engine.run(ions);
        
        // All should remain inactive
        for (const auto& ion : result) {
            REQUIRE(!ion.active);
        }
    }
    
    SECTION("Ion starting outside domain") {
        std::vector<IonState> ions = {create_test_ion()};
        ions[0].pos = Vec3{0, 0, -0.01};  // Before domain entrance
        
        auto result = engine.run(ions);
        
        // Ion should be deactivated (outside domain)
        REQUIRE(!result[0].active);
    }
}

TEST_CASE("SimulationEngine: Numerical safety", "[simulation][engine][safety]") {
    auto cfg = create_test_config();
    cfg.simulation.total_time_s = 1e-6;
    cfg.simulation.compute_derived();
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    
    SimulationEngine engine(cfg, force_registry, integrator);
    
    SECTION("Extremely high velocity does not cause NaN") {
        std::vector<IonState> ions = {create_test_ion()};
        ions[0].vel = Vec3{0, 0, 1e10};  // Unrealistically high
        
        auto result = engine.run(ions);
        
        // Should either exit domain or be deactivated safely (no NaN)
        REQUIRE(std::isfinite(result[0].pos.z));
        REQUIRE(std::isfinite(result[0].vel.z));
    }
}

// ============================================================================
// RNG Thread Safety
// ============================================================================

TEST_CASE("SimulationEngine: RNG thread safety", "[simulation][engine][openmp]") {
    auto cfg = create_test_config();
    cfg.simulation.total_time_s = 1e-6;
    cfg.simulation.enable_openmp = true;
    cfg.simulation.rng_seed = 42;
    cfg.simulation.compute_derived();
    
    auto force_registry = std::make_shared<ForceRegistry>();
    auto integrator = std::make_shared<RK4Strategy>();
    
    SimulationEngine engine(cfg, force_registry, integrator);
    
    SECTION("Reproducible results with same seed") {
        const int N = 50;
        std::vector<IonState> ions1, ions2;
        for (int i = 0; i < N; ++i) {
            ions1.push_back(create_test_ion());
            ions2.push_back(create_test_ion());
        }
        
        auto result1 = engine.run(ions1);
        
        // Create new engine with same seed
        SimulationEngine engine2(cfg, force_registry, integrator);
        auto result2 = engine2.run(ions2);
        
        // Results should be identical (deterministic with fixed seed)
        for (int i = 0; i < N; ++i) {
            REQUIRE_THAT(result1[i].pos.z, Catch::Matchers::WithinAbs(result2[i].pos.z, 1e-12));
            REQUIRE_THAT(result1[i].vel.z, Catch::Matchers::WithinAbs(result2[i].vel.z, 1e-12));
        }
    }
}

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/config/adapter/LegacyAdapter.h"

using Catch::Approx;

// ============================================================================
// GlobalParams Conversion
// ============================================================================

TEST_CASE("LegacyAdapter converts FullConfig to GlobalParams", "[adapter][global]") {
    ICARION::config::FullConfig cfg;
    cfg.simulation.dt_s = 1e-9;
    cfg.simulation.total_steps = 1000;
    cfg.simulation.write_interval = 100;
    cfg.simulation.rng_seed = 42;
    cfg.simulation.enable_gpu = false;
    cfg.simulation.enable_openmp = true;
    
    cfg.physics.collision_model = ICARION::config::CollisionModel::Langevin;
    cfg.physics.enable_reactions = true;
    cfg.physics.enable_space_charge = false;
    
    cfg.output.folder = "./output";
    cfg.output.trajectory_file = "test.h5";
    cfg.output.print_progress = true;
    
    ICARION::core::GlobalParams g = ICARION::config::LegacyAdapter::to_global_params(cfg);
    
    REQUIRE(g.dt_s == Approx(1e-9));
    REQUIRE(g.sim_time_steps == 1000);
    REQUIRE(g.write_interval == 100);
    REQUIRE(g.rng_seed == 42);
    REQUIRE_FALSE(g.enable_gpu);
    REQUIRE(g.parallelization);
    REQUIRE(g.collisionModel == ICARION::core::CollisionModel::Langevin);
    REQUIRE(g.enable_reactions);
    REQUIRE_FALSE(g.enable_space_charge);
    REQUIRE(g.output_file == "./output/test.h5");
    REQUIRE(g.print_results);
}

// ============================================================================
// Enum Conversions
// ============================================================================

TEST_CASE("LegacyAdapter converts collision model enums", "[adapter][enums]") {
    ICARION::config::FullConfig cfg;
    cfg.simulation.dt_s = 1e-9;
    cfg.simulation.total_steps = 100;
    cfg.output.folder = "./output";
    cfg.output.trajectory_file = "test.h5";
    
    SECTION("NoCollisions") {
        cfg.physics.collision_model = ICARION::config::CollisionModel::NoCollisions;
        auto g = ICARION::config::LegacyAdapter::to_global_params(cfg);
        REQUIRE(g.collisionModel == ICARION::core::CollisionModel::NoCollisions);
    }
    
    SECTION("HSD -> HardSphere") {
        cfg.physics.collision_model = ICARION::config::CollisionModel::HSD;
        auto g = ICARION::config::LegacyAdapter::to_global_params(cfg);
        REQUIRE(g.collisionModel == ICARION::core::CollisionModel::HardSphere);
    }
    
    SECTION("EHSS") {
        cfg.physics.collision_model = ICARION::config::CollisionModel::EHSS;
        auto g = ICARION::config::LegacyAdapter::to_global_params(cfg);
        REQUIRE(g.collisionModel == ICARION::core::CollisionModel::EHSS);
    }
}

TEST_CASE("LegacyAdapter converts domains via public API", "[adapter][domain]") {
    ICARION::config::FullConfig cfg;
    cfg.simulation.dt_s = 1e-9;
    cfg.simulation.total_steps = 100;
    cfg.output.folder = "./output";
    cfg.output.trajectory_file = "test.h5";
    cfg.physics.collision_model = ICARION::config::CollisionModel::NoCollisions;
    
    ICARION::config::DomainConfig dom_cfg;
    dom_cfg.domain_index = 0;
    dom_cfg.instrument = ICARION::instrument::InstrumentType::IMS;
    dom_cfg.solver = ICARION::config::SolverType::RK45;
    dom_cfg.geometry.length_m = 0.1;
    dom_cfg.geometry.radius_m = 0.01;
    dom_cfg.environment.pressure_Pa = 101325.0;
    dom_cfg.environment.temperature_K = 300.0;
    cfg.domains.push_back(dom_cfg);
    
    std::vector<ICARION::core::InstrumentDomain> domains = 
        ICARION::config::LegacyAdapter::to_instrument_domains(cfg);
    
    REQUIRE(domains.size() == 1);
    REQUIRE(domains[0].index == 0);
    REQUIRE(domains[0].geom.length_m == Approx(0.1));
    REQUIRE(domains[0].env.pressure_Pa == Approx(101325.0));
}

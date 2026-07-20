// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file test_main_integration.cpp
 * @brief Integration test for main.cpp workflow with SimulationEngine
 * 
 * Tests the complete pipeline using real example configs
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/config/loader/ConfigLoader.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/integrator/SimulationEngine.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/BorisStrategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/collisions/CollisionHandlerFactory.h"
#include "core/physics/reactions/ReactionHandlerFactory.h"

#include <H5Cpp.h>
#include <filesystem>
#include <cmath>
#include <vector>

using namespace ICARION;
using namespace ICARION::config;
using namespace ICARION::integrator;
using namespace ICARION::physics;

TEST_CASE("Main integration: Complete simulation pipeline", "[integration][main]") {
    
    // Use real example config (has proper species DB setup)
    std::string config_path = "../examples/ims_basic.json";
    
    if (!std::filesystem::exists(config_path)) {
        SKIP("Example config not found: " << config_path);
    }
    
    SECTION("Load example configuration") {
        FullConfig config = ConfigLoader::load(config_path);
        
        REQUIRE(config.simulation.total_time_s > 0.0);
        REQUIRE(config.simulation.dt_s > 0.0);
        REQUIRE(!config.domains.empty());
        REQUIRE(!config.ions.species.empty());
        REQUIRE(config.species_db.size() > 0);
    }
    
    SECTION("Generate ions from example config") {
        FullConfig config = ConfigLoader::load(config_path);
        
        std::mt19937 rng(config.simulation.rng_seed);
        auto ions = config.generate_ions(rng);
        
        REQUIRE(ions.size() > 0);
        
        // All ions should be properly initialized
        for (const auto& ion : ions) {
            REQUIRE(ion.active == true);
            REQUIRE(ion.born == true);
            REQUIRE(ion.mass_kg > 0.0);
            REQUIRE(ion.ion_charge_C != 0.0);
        }
    }
    
    SECTION("Create integrator from config") {
        FullConfig config = ConfigLoader::load(config_path);
        
        std::shared_ptr<IIntegrationStrategy> strategy;
        
        // Same logic as main.cpp
        if (config.simulation.integrator == "RK4" || config.simulation.integrator == "rk4") {
            strategy = std::make_shared<RK4Strategy>();
        } else if (config.simulation.integrator == "RK45" || config.simulation.integrator == "rk45") {
            strategy = std::make_shared<RK45Strategy>();
        } else if (config.simulation.integrator == "Boris" || config.simulation.integrator == "boris") {
            strategy = std::make_shared<BorisStrategy>();
        }
        
        REQUIRE(strategy != nullptr);
    }
    
    SECTION("Create physics modules from config") {
        FullConfig config = ConfigLoader::load(config_path);
        
        auto force_registry = std::make_shared<ForceRegistry>(config.domains[0]);
        REQUIRE(force_registry != nullptr);
        
        std::shared_ptr<ICollisionHandler> collision_handler = 
            CollisionHandlerFactory::create(config.physics, nullptr, 0.0, false, &config.species_db);
        REQUIRE(collision_handler != nullptr);
        
        std::shared_ptr<IReactionHandler> reaction_handler = 
            ReactionHandlerFactory::create(config.physics);
        REQUIRE(reaction_handler != nullptr);
    }
    
    SECTION("Complete simulation pipeline (short run)") {
        FullConfig config = ConfigLoader::load(config_path);
        
        // Override to very short simulation for test speed
        config.simulation.total_time_s = 1e-7;  // 100 ns
        config.simulation.dt_s = 1e-9;           // 1 ns steps
        config.simulation.write_interval = 10;
        config.output.trajectory_file = "test_main_integration.h5";  // Relative path
        config.output.folder = ".";
        
        // Generate ions
        std::mt19937 rng(config.simulation.rng_seed);
        auto ions = config.generate_ions(rng);
        REQUIRE(ions.size() > 0);
        
        size_t initial_count = ions.size();
        
        // Create physics modules (same as main.cpp)
        // Phase 12: Create ForceRegistry for each domain
        std::vector<std::shared_ptr<ForceRegistry>> force_registries;
        for (const auto& domain : config.domains) {
            force_registries.push_back(std::make_shared<ForceRegistry>(domain));
        }
        
        std::shared_ptr<IIntegrationStrategy> integration_strategy;
        if (config.simulation.integrator == "RK4" || config.simulation.integrator == "rk4") {
            integration_strategy = std::make_shared<RK4Strategy>();
        } else {
            integration_strategy = std::make_shared<RK45Strategy>();
        }
        
        std::shared_ptr<ICollisionHandler> collision_handler = 
            CollisionHandlerFactory::create(config.physics, nullptr, 0.0, false, &config.species_db);
        std::shared_ptr<IReactionHandler> reaction_handler = 
            ReactionHandlerFactory::create(config.physics);
        
        // Create SimulationEngine
        SimulationEngine engine(
            config,
            force_registries,  // Vector of registries (one per domain)
            integration_strategy,
            collision_handler,
            reaction_handler
        );
        
        // Run simulation
        auto ensemble = core::IonEnsemble::from_legacy(ions);
        auto final_ensemble = engine.run(ensemble);
        
        // Verify results
        REQUIRE(final_ensemble.size() == initial_count);

        // One initial post-step snapshot plus one snapshot per configured
        // interval. The former dual trigger produced adjacent duplicates.
        H5::H5File file(config.output.trajectory_file, H5F_ACC_RDONLY);
        auto time_dataset = file.openDataSet("/trajectory/time");
        auto time_space = time_dataset.getSpace();
        hsize_t time_dims[1];
        time_space.getSimpleExtentDims(time_dims);
        REQUIRE(time_dims[0] == 10);

        std::vector<double> times(time_dims[0]);
        time_dataset.read(times.data(), H5::PredType::NATIVE_DOUBLE);
        for (size_t i = 1; i < times.size(); ++i) {
            REQUIRE(times[i] > times[i - 1]);
            REQUIRE_THAT(
                times[i] - times[i - 1],
                Catch::Matchers::WithinAbs(10e-9, 1e-15)
            );
        }

        file.close();
        std::filesystem::remove(config.output.trajectory_file);
    }
}

TEST_CASE("Main integration: Error handling", "[integration][main][error]") {
    
    SECTION("Invalid config file path") {
        REQUIRE_THROWS_AS(
            ConfigLoader::load("/nonexistent/path/to/config.json"),
            std::exception
        );
    }
}

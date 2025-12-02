// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file test_simulation_engine_soa.cpp
 * @brief Unit tests for SimulationEngine SoA integration (Phase 2)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/integrator/SimulationEngine.h"
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/config/types/FullConfig.h"

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::core;
using namespace ICARION::config;
using namespace ICARION::physics;
using Catch::Approx;

/**
 * @brief Create minimal config for testing
 */
FullConfig create_test_config() {
    FullConfig cfg;
    
    // Simulation settings
    cfg.simulation.total_time_s = 1e-6;  // 1 μs
    cfg.simulation.dt_s = 1e-9;          // 1 ns
    cfg.simulation.write_interval = 100;
    cfg.simulation.enable_openmp = false;
    cfg.simulation.rng_seed = 42;
    cfg.simulation.enable_safety_logging = false;
    
    // Single domain (simple cylindrical geometry)
    DomainConfig domain;
    domain.instrument = Instrument::IMS;
    domain.geometry.length_m = 0.1;      // 10 cm
    domain.geometry.radius_m = 0.01;     // 1 cm
    
    domain.environment.temperature_K = 300.0;
    domain.environment.pressure_Pa = 101325.0;  // 1 atm
    domain.environment.gas_species = "N2";
    domain.environment.compute_derived_properties();  // Compute density, etc.
    
    cfg.domains.push_back(domain);
    
    // Output settings
    cfg.output.folder = "/tmp/test_soa";
    cfg.output.trajectory_file = "test_trajectory.h5";
    cfg.output.print_progress = false;
    
    return cfg;
}

/**
 * @brief Create test ions
 */
std::vector<IonState> create_test_ions(size_t count = 10) {
    std::vector<IonState> ions;
    ions.reserve(count);
    
    for (size_t i = 0; i < count; ++i) {
        IonState ion;
        ion.pos = Vec3{0.0, 0.0, 0.001 * i};  // Spread along z
        ion.vel = Vec3{0.0, 0.0, 100.0};      // 100 m/s forward
        ion.mass_kg = 29.0 * 1.66054e-27;     // ~N2H+
        ion.ion_charge_C = 1.602e-19;
        ion.reduced_mobility_cm2_Vs = 2.0;
        ion.CCS_m2 = 1e-19;
        ion.species_id = "N2H+";
        ion.active = true;
        ion.born = true;
        ion.t = 0.0;
        ions.push_back(ion);
    }
    
    return ions;
}

TEST_CASE("SimulationEngine SoA - Basic construction and run", "[soa][integration]") {
    // Create configuration
    auto config = create_test_config();
    
    // Create force registry (minimal - no forces)
    auto force_registry = std::make_shared<ForceRegistry>(config.domains[0]);
    std::vector<std::shared_ptr<ForceRegistry>> force_registries = {force_registry};
    
    // Create integrator
    auto integrator = std::make_shared<RK4Strategy>();
    
    // Create simulation engine
    SimulationEngine engine(config, force_registries, integrator);
    
    SECTION("AoS vs SoA produces same results") {
        // Create separate test ions for each run to avoid OutputManager reuse
        auto ions_aos = create_test_ions(10);
        auto ions_soa_init = create_test_ions(10);  // Fresh copy for SoA test
        
        // Run with AoS (first engine instance)
        auto result_aos = engine.run(ions_aos);
        
        // Create fresh engine for SoA test to avoid OutputManager finalize issue
        SimulationEngine engine_soa(config, force_registries, integrator);
        
        // Convert to SoA and run
        auto ensemble = IonEnsemble::from_legacy(ions_soa_init);
        auto result_soa = engine_soa.run_soa(ensemble);
        
        // Results should be identical (Phase 2 uses same underlying code)
        REQUIRE(result_aos.size() == result_soa.size());
        
        for (size_t i = 0; i < result_aos.size(); ++i) {
            REQUIRE(result_aos[i].pos.x == Approx(result_soa[i].pos.x).margin(1e-12));
            REQUIRE(result_aos[i].pos.y == Approx(result_soa[i].pos.y).margin(1e-12));
            REQUIRE(result_aos[i].pos.z == Approx(result_soa[i].pos.z).margin(1e-12));
            REQUIRE(result_aos[i].vel.x == Approx(result_soa[i].vel.x).margin(1e-12));
            REQUIRE(result_aos[i].vel.y == Approx(result_soa[i].vel.y).margin(1e-12));
            REQUIRE(result_aos[i].vel.z == Approx(result_soa[i].vel.z).margin(1e-12));
            REQUIRE(result_aos[i].active == result_soa[i].active);
        }
    }
}

TEST_CASE("SimulationEngine SoA - Memory footprint", "[soa][memory]") {
    SECTION("IonEnsemble uses less memory than AoS") {
        auto ions = create_test_ions(1000);
        auto ensemble = IonEnsemble::from_legacy(ions);
        
        // AoS: ~180 bytes/ion (after removing domain cache fields)
        size_t aos_footprint = ions.size() * sizeof(IonState);
        
        // SoA: ~120 bytes/ion
        size_t soa_footprint = ensemble.memory_footprint();
        
        // Verify 30-35% reduction (after Phase 4 domain cache cleanup)
        double reduction = 1.0 - static_cast<double>(soa_footprint) / aos_footprint;
        
        INFO("AoS footprint: " << aos_footprint << " bytes (" << aos_footprint/1000 << " bytes/ion)");
        INFO("SoA footprint: " << soa_footprint << " bytes (" << soa_footprint/1000 << " bytes/ion)");
        INFO("Reduction: " << reduction * 100 << "%");
        
        REQUIRE(reduction > 0.25);  // At least 25% reduction (after Phase 4)
        REQUIRE(reduction < 0.35);  // ~30% reduction expected
        
        // Verify per-ion footprint
        size_t per_ion = soa_footprint / ions.size();
        REQUIRE(per_ion >= 100);
        REQUIRE(per_ion <= 140);
    }
}

TEST_CASE("SimulationEngine SoA - Round-trip conversion", "[soa][conversion]") {
    SECTION("AoS -> SoA -> AoS preserves data") {
        auto ions_original = create_test_ions(50);
        
        // Convert to SoA
        auto ensemble = IonEnsemble::from_legacy(ions_original);
        
        // Convert back to AoS
        auto ions_roundtrip = ensemble.to_legacy();
        
        // Verify identical
        REQUIRE(ions_roundtrip.size() == ions_original.size());
        
        for (size_t i = 0; i < ions_original.size(); ++i) {
            REQUIRE(ions_roundtrip[i].pos.x == Approx(ions_original[i].pos.x));
            REQUIRE(ions_roundtrip[i].pos.y == Approx(ions_original[i].pos.y));
            REQUIRE(ions_roundtrip[i].pos.z == Approx(ions_original[i].pos.z));
            REQUIRE(ions_roundtrip[i].vel.x == Approx(ions_original[i].vel.x));
            REQUIRE(ions_roundtrip[i].vel.y == Approx(ions_original[i].vel.y));
            REQUIRE(ions_roundtrip[i].vel.z == Approx(ions_original[i].vel.z));
            REQUIRE(ions_roundtrip[i].mass_kg == Approx(ions_original[i].mass_kg));
            REQUIRE(ions_roundtrip[i].ion_charge_C == Approx(ions_original[i].ion_charge_C));
            REQUIRE(ions_roundtrip[i].species_id == ions_original[i].species_id);
            REQUIRE(ions_roundtrip[i].active == ions_original[i].active);
            REQUIRE(ions_roundtrip[i].born == ions_original[i].born);
        }
    }
}

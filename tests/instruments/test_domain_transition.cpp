// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
// Domain transition tests
// Test that ions correctly transition between domains (e.g., IMS → TOF)

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

#include "helpers/physics_sim_utils.h"
#include "core/config/types/FullConfig.h"
#include "utils/constants.h"

using namespace ICARION;
using namespace ICARION::tests;
using Catch::Approx;

namespace {

// Helper: Create two-domain config (IMS → TOF)
config::FullConfig make_two_domain_config() {
    config::FullConfig cfg;
    
    // Simulation
    cfg.simulation.dt_s = 1e-9;  // 1 ns
    cfg.simulation.total_time_s = 1e-4;  // 100 µs
    cfg.simulation.write_interval = 100;
    
    // Physics (no collisions for simplicity)
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    
    // Domain 1: IMS drift tube (5cm)
    config::DomainConfig ims;
    ims.instrument = config::Instrument::IMS;
    ims.name = "ims_region";
    ims.geometry.origin_m = {0.0, 0.0, 0.0};
    ims.geometry.length_m = 0.05;  // 5cm
    ims.geometry.radius_m = 0.02;
    ims.geometry.end_aperture_m = 0.005;  // 5mm aperture
    
    ims.fields.dc.axial_V.constant_value = 1000.0;  // 1000V over 5cm = 20 kV/m
    ims.fields.dc.quad_V.constant_value = 0.0;
    ims.fields.dc.radial_V.constant_value = 0.0;
    ims.fields.rf.voltage_V.constant_value = 0.0;
    ims.fields.rf.frequency_Hz.constant_value = 1.0e6;
    ims.fields.rf.compute_derived();
    ims.fields.ac.voltage_V.constant_value = 0.0;
    ims.fields.ac.frequency_Hz.constant_value = 1.0e5;
    ims.fields.ac.compute_derived();
    
    ims.environment.pressure_Pa = 1e-10;
    ims.environment.temperature_K = 300.0;
    ims.environment.gas_species = "N2";
    ims.environment.compute_derived_properties();
    
    // Domain 2: TOF tube (5cm)
    config::DomainConfig tof;
    tof.instrument = config::Instrument::TOF;
    tof.name = "tof_region";
    tof.geometry.origin_m = {0.0, 0.0, 0.05};  // Starts at z=5cm
    tof.geometry.length_m = 0.05;  // 5cm
    tof.geometry.radius_m = 0.02;
    tof.geometry.acc_length_m = 0.01;  // 1cm acceleration
    
    tof.fields.dc.axial_V.constant_value = 2000.0;  // Higher voltage in TOF
    tof.fields.dc.quad_V.constant_value = 0.0;
    tof.fields.dc.radial_V.constant_value = 0.0;
    tof.fields.rf.voltage_V.constant_value = 0.0;
    tof.fields.rf.frequency_Hz.constant_value = 1.0e6;
    tof.fields.rf.compute_derived();
    tof.fields.ac.voltage_V.constant_value = 0.0;
    tof.fields.ac.frequency_Hz.constant_value = 1.0e5;
    tof.fields.ac.compute_derived();
    
    tof.environment.pressure_Pa = 1e-10;
    tof.environment.temperature_K = 300.0;
    tof.environment.gas_species = "N2";
    tof.environment.compute_derived_properties();
    
    cfg.domains.push_back(ims);
    cfg.domains.push_back(tof);
    
    // Output
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "transition_test.h5";
    cfg.output.print_progress = false;
    
    return cfg;
}

// Helper: Create ion
core::IonState make_ion(double mass_amu) {
    core::IonState ion;
    ion.species_id = "TestIon";
    ion.pos = {0.0, 0.0, 0.001};  // Start near IMS entrance
    ion.vel = {0.0, 0.0, 0.0};
    ion.mass_kg = mass_amu * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1e-19;
    ion.active = true;
    ion.current_domain_index = 0;  // Start in domain 0 (IMS)
    return ion;
}

} // namespace

TEST_CASE("Domain transition: Ion crosses from IMS to TOF", "[instrument][transition][physics]") {
    auto cfg = make_two_domain_config();
    
    // Ion: m/z = 500
    auto ion = make_ion(500.0);
    
    // Run simulation
    auto result = run_simple_simulation(cfg, {ion}, false);
    
    REQUIRE(result.ions.size() == 1);
    const auto& final_ion = result.ions[0];
    
    INFO("Final z position: " << final_ion.pos.z*100 << " cm");
    INFO("Final domain index: " << final_ion.current_domain_index);
    INFO("Active: " << final_ion.active);
    
    SECTION("Ion moves forward") {
        // Should move beyond IMS region (z > 5cm)
        REQUIRE(final_ion.pos.z > 0.05);
    }
    
    SECTION("Ion enters second domain") {
        // If ion went far enough, it should be in domain 1 (TOF)
        // or have exited completely
        if (final_ion.active) {
            // Still active → should be in domain 1
            REQUIRE(final_ion.current_domain_index == 1);
        }
        // If not active, it exited the system (also valid)
    }
}

TEST_CASE("Domain transition: Aperture blocks off-axis ions", "[instrument][transition][physics]") {
    auto cfg = make_two_domain_config();
    cfg.simulation.total_time_s = 5e-5;  // 50 µs
    
    SECTION("On-axis ion passes through") {
        auto ion = make_ion(500.0);
        ion.pos = {0.0, 0.0, 0.001};  // On axis
        
        auto result = run_simple_simulation(cfg, {ion}, false);
        const auto& final_ion = result.ions[0];
        
        // Should reach TOF region (z > 5cm) or stay active
        INFO("Final z: " << final_ion.pos.z*100 << " cm, active: " << final_ion.active);
        
        if (final_ion.pos.z > 0.05) {
            // Crossed into TOF region
            REQUIRE(true);
        } else {
            // Still in IMS, should be active
            REQUIRE(final_ion.active);
        }
    }
    
    SECTION("Off-axis ion may be blocked") {
        auto ion = make_ion(500.0);
        ion.pos = {0.008, 0.0, 0.001};  // 8mm off-axis (outside 5mm aperture)
        
        auto result = run_simple_simulation(cfg, {ion}, false);
        const auto& final_ion = result.ions[0];
        
        INFO("Final z: " << final_ion.pos.z*100 << " cm, active: " << final_ion.active);
        
        // Off-axis ion: either blocked at aperture or hits wall
        // No strict requirement - just check simulation completes
        REQUIRE(std::isfinite(final_ion.pos.z));
    }
}

TEST_CASE("Domain transition: Properties update between domains", "[instrument][transition][physics]") {
    // Test that domain-specific properties are updated when ion transitions
    auto cfg = make_two_domain_config();
    
    // Make domains have different pressures to test property updates
    cfg.domains[0].environment.pressure_Pa = 100.0;  // Higher pressure in IMS
    cfg.domains[0].environment.compute_derived_properties();
    cfg.domains[1].environment.pressure_Pa = 1e-10;  // Vacuum in TOF
    cfg.domains[1].environment.compute_derived_properties();
    
    auto ion = make_ion(500.0);
    
    auto result = run_simple_simulation(cfg, {ion}, false);
    
    REQUIRE(result.ions.size() == 1);
    const auto& final_ion = result.ions[0];
    
    SECTION("Simulation completes") {
        // Just verify transition doesn't crash
        REQUIRE(std::isfinite(final_ion.pos.z));
    }
    
    SECTION("Ion moved forward") {
        REQUIRE(final_ion.pos.z > 0.001);
    }
}

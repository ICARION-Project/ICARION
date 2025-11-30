// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

// TOF flight time validation tests
// Simple physics validation: ions with different m/z should have different flight times

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

// Helper: Create minimal TOF config
config::FullConfig make_tof_config(double length_m, double acceleration_V) {
    config::FullConfig cfg;
    
    // Simulation
    cfg.simulation.dt_s = 1e-9;  // 1 ns timestep (faster than 0.1ns)
    cfg.simulation.total_time_s = 1e-4;  // 100 µs max (default, overridden in tests)
    cfg.simulation.write_interval = 1000000;
    
    // Physics (no collisions)
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    
    // Domain
    config::DomainConfig dom;
    dom.instrument = config::Instrument::TOF;
    dom.name = "tof_tube";
    dom.geometry.length_m = length_m;
    dom.geometry.radius_m = 0.02;
    dom.geometry.acc_length_m = 0.01;  // 1cm acceleration region
    
    // Acceleration field
    dom.fields.dc.axial_V.constant_value = acceleration_V;
    dom.fields.dc.quad_V.constant_value = 0.0;
    dom.fields.dc.radial_V.constant_value = 0.0;
    dom.fields.rf.voltage_V.constant_value = 0.0;
    dom.fields.rf.frequency_Hz.constant_value = 1.0e6;
    dom.fields.rf.compute_derived();
    dom.fields.ac.voltage_V.constant_value = 0.0;
    dom.fields.ac.frequency_Hz.constant_value = 1.0e5;
    dom.fields.ac.compute_derived();
    
    // Ultra-low pressure (no collisions)
    dom.environment.pressure_Pa = 1e-10;
    dom.environment.temperature_K = 300.0;
    dom.environment.gas_species = "N2";
    dom.environment.compute_derived_properties();
    
    cfg.domains.push_back(dom);
    
    // Output
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "tof_test.h5";
    cfg.output.print_progress = false;
    
    return cfg;
}

// Helper: Create ion with specified m/z
core::IonState make_ion(double mass_amu, int charge) {
    core::IonState ion;
    ion.species_id = "TestIon";
    ion.pos = {0.0, 0.0, 0.0};
    ion.vel = {0.0, 0.0, 0.0};
    ion.mass_kg = mass_amu * AMU_TO_KG;
    ion.ion_charge_C = charge * ELEM_CHARGE_C;
    ion.CCS_m2 = 1e-19;
    ion.active = true;
    return ion;
}

} // namespace

TEST_CASE("TOF: Ion accelerates through drift tube", "[instrument][tof][physics]") {
    // Setup: 10cm tube, 5000V acceleration
    const double flight_length_m = 0.1;
    const double U_accel = 5000.0;
    auto cfg = make_tof_config(flight_length_m, U_accel);
    cfg.simulation.total_time_s = 1e-4;  // 100 μs
    
    // Ion: m/z = 500
    auto ion = make_ion(500.0, 1);
    
    // Run simulation
    auto result = run_simple_simulation(cfg, {ion}, false);
    
    REQUIRE(result.ions.size() == 1);
    const auto& final_ion = result.ions[0];
    
    SECTION("Ion reaches detector") {
        REQUIRE(final_ion.pos.z >= 0.09);  // Near or past 10cm
    }
    
    SECTION("Ion is accelerated") {
        REQUIRE(final_ion.vel.z > 0.0);
        
        // Expected velocity from energy conservation: v = sqrt(2qU/m)
        double v_expected = std::sqrt(2.0 * final_ion.ion_charge_C * U_accel / final_ion.mass_kg);
        REQUIRE(final_ion.vel.z == Approx(v_expected).epsilon(0.05));
    }
}

TEST_CASE("TOF: Heavier ions arrive later", "[instrument][tof][physics]") {
    // Setup: 10cm tube, 5000V
    const double flight_length_m = 0.1;
    const double U_accel = 5000.0;
    auto cfg = make_tof_config(flight_length_m, U_accel);
    cfg.simulation.total_time_s = 1e-4;  // 100 μs
    
    // Two ions with different masses
    auto ion_light = make_ion(100.0, 1);
    auto ion_heavy = make_ion(1000.0, 1);
    
    // Run separate simulations (no space charge)
    auto result_light = run_simple_simulation(cfg, {ion_light}, false);
    auto result_heavy = run_simple_simulation(cfg, {ion_heavy}, false);
    
    REQUIRE(result_light.ions.size() == 1);
    REQUIRE(result_heavy.ions.size() == 1);
    
    const auto& final_light = result_light.ions[0];
    const auto& final_heavy = result_heavy.ions[0];
    
    SECTION("Both reach detector") {
        REQUIRE(final_light.pos.z >= 0.09);
        REQUIRE(final_heavy.pos.z >= 0.09);
    }
    
    SECTION("Light ion has higher velocity") {
        // For same kinetic energy, lighter ions move faster
        REQUIRE(final_light.vel.z > final_heavy.vel.z);
    }
}

TEST_CASE("TOF: Flight time scales with sqrt(m)", "[instrument][tof][physics]") {
    // Setup: 10cm tube, 5000V
    const double flight_length_m = 0.1;
    const double U_accel = 5000.0;
    auto cfg = make_tof_config(flight_length_m, U_accel);
    cfg.simulation.total_time_s = 1e-4;  // 100 μs
    
    // Test three m/z ratios
    std::vector<double> masses = {100.0, 500.0, 1000.0};
    std::vector<double> velocities;
    
    for (double mass : masses) {
        auto ion = make_ion(mass, 1);
        auto result = run_simple_simulation(cfg, {ion}, false);
        
        REQUIRE(result.ions.size() == 1);
        velocities.push_back(result.ions[0].vel.z);
    }
    
    SECTION("Velocity ratio matches mass ratio") {
        // v ∝ 1/sqrt(m) → v1/v2 = sqrt(m2/m1)
        double v_ratio_12 = velocities[0] / velocities[1];
        double m_ratio_12 = std::sqrt(masses[1] / masses[0]);
        
        REQUIRE(v_ratio_12 == Approx(m_ratio_12).epsilon(0.05));
        
        double v_ratio_23 = velocities[1] / velocities[2];
        double m_ratio_23 = std::sqrt(masses[2] / masses[1]);
        
        REQUIRE(v_ratio_23 == Approx(m_ratio_23).epsilon(0.05));
    }
}

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
// LQIT Mathieu stability validation tests
// Simple physics validation: ions should be trapped in stable region, escape in unstable region

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

// Helper: Create LQIT config
config::FullConfig make_lqit_config(double r0_m, double rf_voltage_V, double rf_frequency_Hz) {
    config::FullConfig cfg;
    
    // Simulation
    cfg.simulation.dt_s = 1e-8;  // 10 ns (100 steps per RF cycle at 1MHz)
    cfg.simulation.total_time_s = 200e-6;  // 200 µs (200 RF cycles)
    cfg.simulation.write_interval = 100000;
    
    // Physics (no collisions)
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    
    // Domain
    config::DomainConfig dom;
    dom.instrument = config::Instrument::LQIT;
    dom.name = "lqit_trap";
    dom.geometry.radius_m = r0_m;
    dom.geometry.length_m = 0.1;
    
    // RF trapping field
    dom.fields.rf.voltage_V.constant_value = rf_voltage_V;
    dom.fields.rf.frequency_Hz.constant_value = rf_frequency_Hz;
    dom.fields.rf.compute_derived();
    
    // DC off
    dom.fields.dc.quad_V.constant_value = 0.0;
    dom.fields.dc.axial_V.constant_value = 0.0;
    dom.fields.dc.radial_V.constant_value = 0.0;
    dom.fields.ac.voltage_V.constant_value = 0.0;
    dom.fields.ac.frequency_Hz.constant_value = 1.0e5;
    dom.fields.ac.compute_derived();
    
    // Ultra-low pressure
    dom.environment.pressure_Pa = 1e-10;
    dom.environment.temperature_K = 300.0;
    dom.environment.gas_species = "N2";
    dom.environment.compute_derived_properties();
    
    cfg.domains.push_back(dom);
    
    // Output
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "lqit_test.h5";
    cfg.output.print_progress = false;
    
    return cfg;
}

// Helper: Create ion with displacement
core::IonState make_ion(double mass_amu, double displacement_m) {
    core::IonState ion;
    ion.species_id = "TestIon";
    ion.pos = {displacement_m, 0.0, 0.0};  // Off-axis
    ion.vel = {0.0, 0.0, 0.0};
    ion.mass_kg = mass_amu * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1e-19;
    ion.active = true;
    return ion;
}

// Helper: Calculate Mathieu q parameter
double calculate_q(double mass_kg, double V_rf, double omega, double r0) {
    return (4.0 * ELEM_CHARGE_C * V_rf) / (mass_kg * omega * omega * r0 * r0);
}

// Helper: Radial distance
double radial_distance(const core::IonState& ion) {
    return std::sqrt(ion.pos.x * ion.pos.x + ion.pos.y * ion.pos.y);
}

} // namespace

TEST_CASE("LQIT: Stable trapping (low q parameter)", "[instrument][lqit][physics]") {
    // Setup: r0=5mm, RF 1MHz @ 100V → q ≈ 0.078 for m/z=500 (stable)
    double r0 = 0.005;
    double V_rf = 100.0;
    double f_rf = 1.0e6;
    
    auto cfg = make_lqit_config(r0, V_rf, f_rf);
    
    // Ion: m/z = 500, 1mm off-axis
    auto ion = make_ion(500.0, 0.001);
    
    // Calculate q parameter
    double omega = 2.0 * M_PI * f_rf;
    double q = calculate_q(ion.mass_kg, V_rf, omega, r0);
    
    INFO("Mathieu q = " << q);
    REQUIRE(q < 0.908);  // Must be in stable region
    
    // Run simulation
    auto result = run_simple_simulation(cfg, {ion}, false);
    
    REQUIRE(result.ions.size() == 1);
    const auto& final_ion = result.ions[0];
    
    SECTION("Ion remains confined") {
        double r_final = radial_distance(final_ion);
        REQUIRE(r_final < r0);  // Still inside trap
    }
    
    SECTION("Ion is still active") {
        REQUIRE(final_ion.active);
    }
}

TEST_CASE("LQIT: Unstable region (high q parameter)", "[instrument][lqit][physics]") {
    // Setup: Higher voltage to push into unstable region
    // For q > 0.908: need V > ~1160V (at 1MHz, r0=5mm, m/z=500)
    double r0 = 0.005;
    double V_rf = 1500.0;  // High RF → q ≈ 1.17 (unstable)
    double f_rf = 1.0e6;
    
    auto cfg = make_lqit_config(r0, V_rf, f_rf);
    cfg.simulation.total_time_s = 100e-6;  // 100 µs (enough to see escape)
    
    // Ion: m/z = 500, 1mm off-axis
    auto ion = make_ion(500.0, 0.001);
    
    // Calculate q parameter
    double omega = 2.0 * M_PI * f_rf;
    double q = calculate_q(ion.mass_kg, V_rf, omega, r0);
    
    INFO("Mathieu q = " << q);
    INFO("Expected q ≈ " << 4.0*ELEM_CHARGE_C*V_rf/(500*AMU_TO_KG*omega*omega*r0*r0));
    REQUIRE(q > 0.908);  // Must be in unstable region
    
    // Run simulation
    auto result = run_simple_simulation(cfg, {ion}, false);
    
    REQUIRE(result.ions.size() == 1);
    const auto& final_ion = result.ions[0];
    
    SECTION("Ion escapes or is deactivated") {
        double r_final = radial_distance(final_ion);
        
        // Ion should either escape (r > r0) or be deactivated
        bool escaped = (r_final >= r0);
        bool deactivated = !final_ion.active;
        
        REQUIRE((escaped || deactivated));
    }
}

TEST_CASE("LQIT: Mass-selective stability", "[instrument][lqit][physics]") {
    // Setup: r0=5mm, RF 1MHz @ 250V
    double r0 = 0.005;
    double V_rf = 250.0;
    double f_rf = 1.0e6;
    double omega = 2.0 * M_PI * f_rf;
    
    auto cfg = make_lqit_config(r0, V_rf, f_rf);
    cfg.simulation.total_time_s = 100e-6;  // 100 µs
    
    // Test different masses
    std::vector<double> masses = {200.0, 500.0, 1000.0};
    
    for (double mass : masses) {
        auto ion = make_ion(mass, 0.001);
        
        double q = calculate_q(ion.mass_kg, V_rf, omega, r0);
        bool should_be_stable = (q < 0.908);
        
        INFO("m/z = " << mass << ", q = " << q << ", stable = " << should_be_stable);
        
        auto result = run_simple_simulation(cfg, {ion}, false);
        REQUIRE(result.ions.size() == 1);
        
        const auto& final_ion = result.ions[0];
        double r_final = radial_distance(final_ion);
        
        if (should_be_stable) {
            // Ion should remain confined
            REQUIRE(r_final < r0);
        }
        // Note: We don't test unstable case here as it's already covered above
    }
}

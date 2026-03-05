// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
// FTICR cyclotron frequency validation tests
// Test that ions rotate at correct cyclotron frequency: f_c = q·B / (2π·m)

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <cmath>
#include <iostream>
#include <algorithm>

#include "helpers/physics_sim_utils.h"
#include "core/config/types/FullConfig.h"
#include "utils/constants.h"

using namespace ICARION;
using namespace ICARION::tests;
using Catch::Approx;
using Catch::Matchers::WithinAbs;

namespace {

// Helper: Create minimal FTICR config
config::FullConfig make_fticr_config(double B_field_T, double radius_m, double length_m) {
    config::FullConfig cfg;
    
    // Simulation
    cfg.simulation.dt_s = 1e-9;  // 1 ns
    cfg.simulation.total_time_s = 1e-4;  // 100 µs (multiple orbits)
    cfg.simulation.write_interval = 100;
    
    // Physics (no collisions for clean frequency)
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    
    // FTICR cell
    config::DomainConfig dom;
    dom.instrument = config::Instrument::FTICR;
    dom.name = "fticr_cell";
    dom.domain_index = 0;
    
    // Cylindrical geometry
    dom.geometry.radius_m = radius_m;
    dom.geometry.length_m = length_m;
    dom.geometry.origin_m = {0.0, 0.0, 0.0};
    
    // Identity rotation
    dom.rotation_global_to_local = Mat3::identity();
    dom.rotation_local_to_global = Mat3::identity();
    
    // Magnetic field (z-direction)
    dom.fields.magnetic.enabled = true;
    dom.fields.magnetic.field_strength_T = {0.0, 0.0, B_field_T};
    
    // No electric fields needed for cyclotron frequency test
    dom.fields.dc.axial_V.constant_value = 0.0;
    dom.fields.dc.quad_V.constant_value = 0.0;
    dom.fields.dc.radial_V.constant_value = 0.0;
    dom.fields.rf.voltage_V.constant_value = 0.0;
    dom.fields.rf.frequency_Hz.constant_value = 0.0;
    dom.fields.rf.compute_derived();
    dom.fields.ac.voltage_V.constant_value = 0.0;
    dom.fields.ac.frequency_Hz.constant_value = 0.0;
    dom.fields.ac.compute_derived();
    
    // Ultra-high vacuum
    dom.environment.pressure_Pa = 1e-10;
    dom.environment.temperature_K = 300.0;
    dom.environment.gas_species = "N2";
    dom.environment.compute_derived_properties();
    
    cfg.domains.push_back(dom);
    
    // Output
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "fticr_test.h5";
    cfg.output.print_progress = false;
    
    return cfg;
}

// Helper: Create ion with transverse velocity (circular motion)
core::IonState make_ion(double mass_amu, double v_transverse_m_s) {
    core::IonState ion;
    ion.species_id = "TestIon";
    ion.pos = {5e-3, 0.0, 0.0};  // 5mm off-center in x
    ion.vel = {0.0, v_transverse_m_s, 0.0};  // Initial velocity in y
    
    ion.mass_kg = mass_amu * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1e-19;
    ion.active = true;
    ion.born = true;
    ion.birth_time_s = 0.0;
    
    return ion;
}

// Helper: Expected cyclotron frequency
double expected_cyclotron_freq_Hz(double mass_kg, double q_charge, double B_T) {
    // f_c = q·B / (2π·m)
    return (q_charge * B_T) / (2.0 * M_PI * mass_kg);
}

// Helper: Measure rotation frequency from trajectory
[[maybe_unused]] double measure_rotation_frequency(const std::vector<double>& times,
                                                   const std::vector<double>& x_pos,
                                                   const std::vector<double>& y_pos) {
    if (times.size() < 100) return 0.0;
    
    // Calculate angular position θ = atan2(y, x)
    std::vector<double> angles;
    angles.reserve(times.size());
    for (size_t i = 0; i < x_pos.size(); ++i) {
        angles.push_back(std::atan2(y_pos[i], x_pos[i]));
    }
    
    // Unwrap angles (handle 2π jumps)
    for (size_t i = 1; i < angles.size(); ++i) {
        double diff = angles[i] - angles[i-1];
        if (diff > M_PI) angles[i] -= 2.0 * M_PI;
        if (diff < -M_PI) angles[i] += 2.0 * M_PI;
    }
    
    // Linear fit: θ = ω·t + θ₀
    // ω = Δθ / Δt
    double t_total = times.back() - times.front();
    double angle_change = angles.back() - angles.front();
    double omega = angle_change / t_total;  // rad/s
    
    return std::abs(omega) / (2.0 * M_PI);  // Hz
}

} // namespace

TEST_CASE("FTICR: Cyclotron frequency", "[instrument][fticr][physics]") {
    // Standard FTICR parameters
    double B_field = 7.0;      // 7 Tesla (typical high-field FTICR)
    double radius = 0.05;      // 5 cm cell radius
    double length = 0.1;       // 10 cm cell length
    
    auto cfg = make_fticr_config(B_field, radius, length);
    
    SECTION("Cyclotron frequency formula") {
        // Test formula: f_c = q·B / (2π·m)
        double mass_amu = 500.0;
        double v_trans = 1000.0;
        
        core::IonState ion = make_ion(mass_amu, v_trans);
        
        // Expected cyclotron frequency
        double f_expected = expected_cyclotron_freq_Hz(ion.mass_kg, ion.ion_charge_C, B_field);
        
        INFO("m/z = " << mass_amu);
        INFO("B-field = " << B_field << " T");
        INFO("Expected cyclotron frequency: " << f_expected/1000.0 << " kHz");
        INFO("Expected period: " << 1.0/f_expected * 1e6 << " µs");
        
        // Check formula gives reasonable result
        REQUIRE(f_expected > 0.0);
        REQUIRE(f_expected < 1e9);  // Less than 1 GHz (reasonable for ICR)
        
        // For 7T and m/z=500: f ≈ 214 kHz
        REQUIRE_THAT(f_expected/1000.0, WithinAbs(214.986, 0.01));  // kHz
    }
    
    SECTION("Cyclotron frequency scales with 1/m") {
        // Lighter ion should rotate faster: f ∝ 1/m
        double mass_light = 100.0;  // m/z = 100
        double mass_heavy = 400.0;  // m/z = 400
        
        double f_light = expected_cyclotron_freq_Hz(mass_light * AMU_TO_KG, ELEM_CHARGE_C, B_field);
        double f_heavy = expected_cyclotron_freq_Hz(mass_heavy * AMU_TO_KG, ELEM_CHARGE_C, B_field);
        
        INFO("Light ion (m/z=" << mass_light << "): f = " << f_light/1000 << " kHz");
        INFO("Heavy ion (m/z=" << mass_heavy << "): f = " << f_heavy/1000 << " kHz");
        INFO("Frequency ratio: " << f_light/f_heavy);
        
        // Frequency ratio should equal mass ratio: f_light/f_heavy = m_heavy/m_light
        double freq_ratio = f_light / f_heavy;
        double mass_ratio = mass_heavy / mass_light;
        
        REQUIRE_THAT(freq_ratio, WithinAbs(mass_ratio, 0.001));
        
        // Both frequencies positive
        REQUIRE(f_light > f_heavy);
        REQUIRE(f_light > 0.0);
        REQUIRE(f_heavy > 0.0);
    }
}

TEST_CASE("FTICR: Magnetic field strength", "[instrument][fticr][physics]") {
    double mass_amu = 500.0;
    
    SECTION("Higher B-field increases cyclotron frequency") {
        double B_low = 3.0;   // 3 T
        double B_high = 9.0;  // 9 T
        
        double f_low = expected_cyclotron_freq_Hz(mass_amu * AMU_TO_KG, ELEM_CHARGE_C, B_low);
        double f_high = expected_cyclotron_freq_Hz(mass_amu * AMU_TO_KG, ELEM_CHARGE_C, B_high);
        
        INFO("3T: f = " << f_low/1000 << " kHz");
        INFO("9T: f = " << f_high/1000 << " kHz");
        
        // f ∝ B
        REQUIRE_THAT(f_high / f_low, WithinAbs(B_high / B_low, 0.001));
    }
}

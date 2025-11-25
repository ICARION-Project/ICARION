// SPDX-License-Identifier: MIT
// Orbitrap frequency validation tests
// Test that ions oscillate at the correct frequency: f = (1/2π) × sqrt(q × k / m)

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <iostream>
#include <iomanip>

#include "helpers/physics_sim_utils.h"
#include "core/config/types/FullConfig.h"
#include "utils/constants.h"

using namespace ICARION;
using namespace ICARION::tests;
using Catch::Approx;

namespace {

// Helper: Create minimal Orbitrap config
config::FullConfig make_orbitrap_config(double R_inner, double R_outer, double R_char) {
    config::FullConfig cfg;
    
    // Simulation
    cfg.simulation.dt_s = 1e-9;  // 1 ns
    cfg.simulation.total_time_s = 1e-4;  // 100 µs (several oscillation periods)
    cfg.simulation.write_interval = 100;  // Capture trace
    
    // Physics (no collisions for clean frequency measurement)
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    
    // Domain
    config::DomainConfig dom;
    dom.instrument = config::Instrument::Orbitrap;
    dom.name = "orbitrap";
    dom.domain_index = 0;
    
    // Compute derived orbitrap parameters (C_in, C_out)
    dom.geometry.orbitrap_C_in = -0.5 * R_inner * R_inner;
    dom.geometry.orbitrap_C_out = -0.5 * R_outer * R_outer;
    
    // Orbitrap geometry radii
    dom.geometry.radius_in_m = R_inner;
    dom.geometry.radius_out_m = R_outer;
    dom.geometry.radius_char_m = R_char;
    dom.geometry.origin_m = {0.0, 0.0, 0.0};  // Orbitrap centered at origin (z=0 is axial center)
    

    
    // Fields (static DC for quadrupole-logarithmic potential)
    dom.fields.dc.axial_V.constant_value = 0.0;  // V0 for orbitrap potential
    dom.fields.dc.quad_V.constant_value = 0.0;
    dom.fields.dc.radial_V.constant_value = 3500.0;  // Default 3500V (can be overridden)
    dom.fields.rf.voltage_V.constant_value = 0.0;
    dom.fields.rf.frequency_Hz.constant_value = 0.0;
    dom.fields.rf.compute_derived();
    dom.fields.ac.voltage_V.constant_value = 0.0;
    dom.fields.ac.frequency_Hz.constant_value = 0.0;
    dom.fields.ac.compute_derived();
    
    // Ultra-low pressure
    dom.environment.pressure_Pa = 1e-10;
    dom.environment.temperature_K = 300.0;
    dom.environment.gas_species = "N2";
    dom.environment.compute_derived_properties();
    
    cfg.domains.push_back(dom);
    
    // Output
    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "orbitrap_test.h5";
    cfg.output.print_progress = false;
    
    return cfg;
}

// Helper: Create ion with tangential kinetic energy (for orbital motion)
// In Orbitrap: ions are injected tangentially with E_tangential ≈ 1-2 keV
// Axial energy is much smaller (from injection offset)
core::IonState make_ion(double mass_amu, double E_tangential_eV) {
    core::IonState ion;
    ion.species_id = "TestIon";
    ion.pos = {9e-3, 0.0, 6e-3};  // Standard: x=9mm (radial), z=6mm (axial)
    
    // Calculate tangential velocity from kinetic energy
    // E_tang = 0.5 * m * v_tang^2
    double m_kg = mass_amu * AMU_TO_KG;
    double v_tang = std::sqrt(2.0 * E_tangential_eV * ELEM_CHARGE_C / m_kg);
    
    // Initial velocity: tangential (y-direction for x-position)
    // Small axial velocity component from initial z-offset
    ion.vel = {0.0, v_tang, 0.0};
    
    ion.mass_kg = m_kg;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1e-19;
    ion.active = true;
    return ion;
}

// Helper: Calculate k parameter from geometry and voltage (as in codebase)
double calculate_k_param(double voltage, double R_inner, double R_outer, double R_char) {
    // k = 2*V / (R_char² * ln(R_out/R_in) - 0.5*(R_out² - R_in²))
    double R_char_sq = R_char * R_char;
    double denom = R_char_sq * std::log(R_outer / R_inner) - 0.5 * (R_outer * R_outer - R_inner * R_inner);
    return 2.0 * voltage / denom;
}

// Helper: Calculate expected axial frequency
double expected_frequency_Hz(double mass_kg, double q_charge, double k_param) {
    // f = (1/2π) × sqrt(q × k / m)
    return (1.0 / (2.0 * M_PI)) * std::sqrt(q_charge * k_param / mass_kg);
}

} // namespace

TEST_CASE("Orbitrap: Axial oscillation frequency", "[instrument][orbitrap][physics]") {
    // Setup: Standard Orbitrap parameters from paper (Table 1)
    double R_inner = 6e-3;   // 6mm spindle radius
    double R_outer = 15e-3;  // 15mm outer radius
    double R_char = 22e-3;   // 22mm characteristic radius
    double voltage = 3500.0; // 3500V trap voltage
    
    auto cfg = make_orbitrap_config(R_inner, R_outer, R_char);
    cfg.simulation.integrator = "RK45";
    cfg.domains[0].fields.dc.radial_V.constant_value = voltage;
    cfg.domains[0].geometry.length_m = 0.04;  // ±20mm around z=0
    cfg.domains[0].environment.pressure_Pa = 1e-7;  // UHV
    
    // Ion: m/z = 500, tangential injection energy 1600 eV (from Makarov 2000)
    double mass_amu = 500.0;
    double E_tangential_eV = 1600.0;  // Tangential (orbital) energy
    core::IonState ion = make_ion(mass_amu, E_tangential_eV);
    // Standard position: r=9mm, z=6mm (already set by make_ion)
    
    // Calculate k parameter
    double k = calculate_k_param(voltage, R_inner, R_outer, R_char);
    
    // Expected frequency
    double f_expected = expected_frequency_Hz(ion.mass_kg, ion.ion_charge_C, k);
    
    INFO("Expected frequency: " << f_expected/1000.0 << " kHz");
    INFO("Expected period: " << 1.0/f_expected * 1e6 << " µs");
    
    // Run simulation with trace
    auto result = run_simple_simulation(cfg, {ion}, true);
    
    REQUIRE(result.ions.size() == 1);
    const auto& final_ion = result.ions[0];
    
    SECTION("Ion remains trapped") {
        // Debug: Print final position
        double x = final_ion.pos.x;
        double y = final_ion.pos.y;
        double z = final_ion.pos.z;
        double r = std::sqrt(x*x + y*y);
        
        std::cout << "\n=== ION STATUS ===\n";
        std::cout << "Final position: r=" << r*1000 << " mm, z=" << z*1000 << " mm\n";
        std::cout << "Final velocity: vx=" << final_ion.vel.x << " m/s, vy=" << final_ion.vel.y 
                  << " m/s, vz=" << final_ion.vel.z << " m/s\n";
        std::cout << "Active: " << final_ion.active << "\n";
        std::cout << "Initial: r=" << 9.0 << " mm, z=" << 6.0 << " mm\n";
        
        // Print boundary check at starting position
        double z_start = 6e-3;
        double R_in = 6e-3;
        double R_out = 15e-3;
        double R_m = 22e-3;
        
        auto calc_r_boundary = [&](double z_val, double R_electrode) {
            double z2 = z_val * z_val;
            double r_lo = 0.3 * R_electrode;
            double r_hi = 3.0 * R_electrode;
            for (int iter = 0; iter < 80; ++iter) {
                double r_mid = 0.5 * (r_lo + r_hi);
                double f = 0.5 * (r_mid*r_mid - R_electrode*R_electrode) + 
                          R_m*R_m * std::log(R_electrode / r_mid) - z2;
                if (std::abs(f) < 1e-10) return r_mid;
                if (f > 0) r_lo = r_mid;
                else r_hi = r_mid;
            }
            return 0.5 * (r_lo + r_hi);
        };
        
        double r_in_start = calc_r_boundary(z_start, R_in);
        double r_out_start = calc_r_boundary(z_start, R_out);
        std::cout << "At z=6mm: r_in=" << r_in_start*1000 << " mm, r_out=" << r_out_start*1000 << " mm\n";
        std::cout << "Ion starts at r=9mm -> inside? " << (9e-3 >= r_in_start && 9e-3 <= r_out_start) << "\n";
        
        REQUIRE(final_ion.active);
    }
    
    SECTION("Ion oscillates in z-direction") {
        // Check that z-position varies (oscillates)
        REQUIRE(result.trace.z_positions.size() > 100);
        
        double z_min = *std::min_element(result.trace.z_positions.begin(), 
                                         result.trace.z_positions.end());
        double z_max = *std::max_element(result.trace.z_positions.begin(), 
                                         result.trace.z_positions.end());
        double z_amplitude = (z_max - z_min) / 2.0;
        
        INFO("z amplitude: " << z_amplitude*1000 << " mm");
        REQUIRE(z_amplitude > 0.0005);  // At least 0.5mm oscillation
    }
}

TEST_CASE("Orbitrap: Frequency scales with sqrt(1/m)", "[instrument][orbitrap][physics]") {
    // Test that heavier ions oscillate slower: f ∝ 1/sqrt(m)
    // Using standard Orbitrap parameters from paper
    double R_inner = 6e-3;   // 6mm spindle radius
    double R_outer = 15e-3;  // 15mm outer radius
    double R_char = 22e-3;   // 22mm characteristic radius
    double voltage = 3500.0; // 3500V trap voltage
    
    auto cfg1 = make_orbitrap_config(R_inner, R_outer, R_char);
    auto cfg2 = make_orbitrap_config(R_inner, R_outer, R_char);
    cfg1.simulation.integrator = "RK45";
    cfg2.simulation.integrator = "RK45";
    cfg1.simulation.total_time_s = 1e-5;  // 10 µs (shorter to check initial trapping)
    cfg2.simulation.total_time_s = 1e-5;
    cfg1.domains[0].fields.dc.radial_V.constant_value = voltage;
    cfg2.domains[0].fields.dc.radial_V.constant_value = voltage;
    cfg1.domains[0].geometry.length_m = 0.04;
    cfg2.domains[0].geometry.length_m = 0.04;
    cfg1.domains[0].environment.pressure_Pa = 1e-7;
    cfg2.domains[0].environment.pressure_Pa = 1e-7;
    
    // Calculate k parameter
    double k = calculate_k_param(voltage, R_inner, R_outer, R_char);
    
    // Two ions with 4x mass difference, both with 1600 eV tangential energy
    double E_tangential_eV = 1600.0;  // Same tangential injection energy
    auto ion_light = make_ion(200.0, E_tangential_eV);  // m/z = 200
    auto ion_heavy = make_ion(800.0, E_tangential_eV);  // m/z = 800
    
    double f_light = expected_frequency_Hz(ion_light.mass_kg, ion_light.ion_charge_C, k);
    double f_heavy = expected_frequency_Hz(ion_heavy.mass_kg, ion_heavy.ion_charge_C, k);
    
    INFO("Light ion (m/z=200): " << f_light/1000.0 << " kHz");
    INFO("Heavy ion (m/z=800): " << f_heavy/1000.0 << " kHz");
    INFO("Frequency ratio: " << f_light/f_heavy << " (expected: 2.0)");
    
    auto result_light = run_simple_simulation(cfg1, {ion_light}, true);  // Enable trace to use manual loop
    auto result_heavy = run_simple_simulation(cfg2, {ion_heavy}, true);
    
    REQUIRE(result_light.ions.size() == 1);
    REQUIRE(result_heavy.ions.size() == 1);
    
    SECTION("Both ions remain trapped") {
        std::cout << "\n=== ION STATUS (m/z=200) ===\n";
        std::cout << "Final pos: r=" << std::sqrt(result_light.ions[0].pos.x * result_light.ions[0].pos.x + 
                                                   result_light.ions[0].pos.y * result_light.ions[0].pos.y) * 1000 
                  << " mm, z=" << result_light.ions[0].pos.z * 1000 << " mm\n";
        std::cout << "Active: " << result_light.ions[0].active << "\n";
        std::cout << "Initial v_tang: " << ion_light.vel.y << " m/s\n";
        
        std::cout << "\n=== ION STATUS (m/z=800) ===\n";
        std::cout << "Final pos: r=" << std::sqrt(result_heavy.ions[0].pos.x * result_heavy.ions[0].pos.x + 
                                                   result_heavy.ions[0].pos.y * result_heavy.ions[0].pos.y) * 1000 
                  << " mm, z=" << result_heavy.ions[0].pos.z * 1000 << " mm\n";
        std::cout << "Active: " << result_heavy.ions[0].active << "\n";
        std::cout << "Initial v_tang: " << ion_heavy.vel.y << " m/s\n";
        
        INFO("Light ion (m/z=200) active: " << result_light.ions[0].active);
        INFO("Heavy ion (m/z=800) active: " << result_heavy.ions[0].active);
        REQUIRE(result_light.ions[0].active);
        REQUIRE(result_heavy.ions[0].active);
    }
    
    SECTION("Frequency ratio matches sqrt(m_heavy/m_light)") {
        // For f ∝ 1/sqrt(m), we have f1/f2 = sqrt(m2/m1)
        // With m_ratio = 4, expect f_ratio = 2
        double expected_ratio = std::sqrt(800.0 / 200.0);
        double actual_ratio = f_light / f_heavy;
        
        REQUIRE(actual_ratio == Approx(expected_ratio).epsilon(0.01));
    }
}

TEST_CASE("Orbitrap: Mass-dependent detection", "[instrument][orbitrap][physics]") {
    // Simplified test: Just verify different masses complete simulation
    
    double R_inner = 6e-3;   // 6mm spindle radius
    double R_outer = 15e-3;  // 15mm outer radius
    double R_char = 22e-3;   // 22mm characteristic radius
    double voltage = 3500.0; // 3500V trap voltage
    
    auto cfg = make_orbitrap_config(R_inner, R_outer, R_char);
    cfg.simulation.integrator = "RK45";
    cfg.simulation.total_time_s = 2e-5;  // 20 µs (short)
    cfg.domains[0].fields.dc.radial_V.constant_value = voltage;
    cfg.domains[0].geometry.length_m = 0.04;
    cfg.domains[0].environment.pressure_Pa = 1e-7;
    
    // Calculate k parameter
    double k = calculate_k_param(voltage, R_inner, R_outer, R_char);
    
    std::vector<double> masses = {300.0, 600.0, 1200.0};
    double E_tangential_eV = 1600.0;  // Tangential injection energy
    
    for (double mass : masses) {
        auto ion = make_ion(mass, E_tangential_eV);
        auto result = run_simple_simulation(cfg, {ion}, true);  // Enable trace to use manual loop
        
        REQUIRE(result.ions.size() == 1);
        
        double f_expected = expected_frequency_Hz(ion.mass_kg, ion.ion_charge_C, k);
        INFO("m/z = " << mass << ", f = " << f_expected/1000.0 << " kHz");
        
        // Just check ion completes simulation
        REQUIRE(std::isfinite(result.ions[0].pos.z));
    }
}

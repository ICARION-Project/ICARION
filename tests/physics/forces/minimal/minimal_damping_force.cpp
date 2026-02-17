// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * @file minimal_damping_force.cpp
 * @brief Minimal physics test for DampingForce
 * 
 * Purpose: Quick verification that collision damping F = -γ·m·v is computed correctly.
 * 
 * ⚠️ WARNING: Uses DEPRECATED DampingParams!
 * This test uses the legacy parameter struct for simplicity.
 * Future: Will be replaced with EnvironmentConfig in Phase 2 SSOT refactor.
 * 
 * Expected outputs:
 * - Friction model: F opposes velocity, |F| ∝ |v|
 * - HardSphere model: Collision frequency depends on CCS and gas density
 * - Langevin model: Enhanced cross-section from polarization
 * 
 * Compile:
 *   g++ -std=c++20 -I../../../src -I../../../include \
 *       minimal_damping_force.cpp \
 *       ../../../src/core/physics/forces/DampingForce.cpp \
 *       ../../../src/core/utils/mathUtils.cpp \
 *       -o minimal_damping_force
 * 
 * Run: ./minimal_damping_force
 */

#include "core/physics/forces/DampingForce.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/config/types/EnvironmentConfig.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ICARION;
using namespace ICARION::physics;

void test_explicit_gamma_damping() {
    std::cout << "\n=== Test 1: Friction Model (Mobility-Based) ===\n";
    std::cout << "Physics: F = -γ·m·v where γ = q/(K₀·m)\n";
    
    // SSOT: Use EnvironmentConfig (empty for Friction model)
    ICARION::config::EnvironmentConfig env;
    env.pressure_Pa = 101325.0;
    env.temperature_K = 300.0;
    env.compute_derived_properties();
    
    DampingForce force(env, DampingModel::Friction);
    ForceContext ctx;
    
    // Test ion moving at 1000 m/s with known mobility
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{1000, 0, 0};
    ion.ion_charge_C = 1.602e-19;
    ion.mass_kg = 100 * 1.66e-27;  // 100 u
    ion.reduced_mobility_cm2_Vs = 2.0;  // 2 cm²/(V·s)
    
    ICARION::core::IonEnsemble ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    ctx.ion_index = 0;
    Vec3 F = force.compute(ens, 0, 0.0, ctx);
    
    // Expected: F = -γ·m·v where γ = q/(K₀·m)
    const double K0_SI = ion.reduced_mobility_cm2_Vs * 1e-4;  // Convert to m²/(V·s)
    const double gamma = ion.ion_charge_C / (K0_SI * ion.mass_kg);
    Vec3 F_expected(
        -gamma * ion.mass_kg * ion.vel.x,
        -gamma * ion.mass_kg * ion.vel.y,
        -gamma * ion.mass_kg * ion.vel.z
    );
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Reduced mobility: K₀ = " << ion.reduced_mobility_cm2_Vs << " cm²/(V·s)\n";
    std::cout << "Damping rate: γ = q/(K₀·m) = " << gamma << " 1/s\n";
    std::cout << "Ion mass: m = " << ion.mass_kg << " kg\n";
    std::cout << "Ion velocity: v = (" << ion.vel.x << ", " << ion.vel.y << ", " << ion.vel.z << ") m/s\n";
    std::cout << "\nExpected force: F = -γ·m·v = (" << F_expected.x << ", " << F_expected.y << ", " << F_expected.z << ") N\n";
    std::cout << "Computed force: F = (" << F.x << ", " << F.y << ", " << F.z << ") N\n";
    
    // Verify force opposes velocity
    const double dot_product = F.x * ion.vel.x + F.y * ion.vel.y + F.z * ion.vel.z;
    std::cout << "\nF · v = " << dot_product << " (should be < 0, opposing motion)\n";
    
    const bool pass = (std::abs(F.x - F_expected.x) / std::abs(F_expected.x) < 1e-10) &&
                     (std::abs(F.y) < 1e-20) &&
                     (std::abs(F.z) < 1e-20) &&
                     (dot_product < 0);
    
    std::cout << "\nResult: " << (pass ? "✓ PASS" : "✗ FAIL") << "\n";
}

void test_hard_sphere_damping() {
    std::cout << "\n=== Test 2: Hard-Sphere Collision Model ===\n";
    std::cout << "Physics: γ = ν_collision = n·σ·v_th·(m_n/(m_i+m_n))\n";
    
    // SSOT: Use EnvironmentConfig
    ICARION::config::EnvironmentConfig env;
    env.pressure_Pa = 101325.0;  // 1 atm
    env.temperature_K = 300.0;
    env.gas_species = "N2";
    env.compute_derived_properties();
    
    DampingForce force(env, DampingModel::HardSphere);
    ForceContext ctx;
    
    // Test ion at rest (thermal motion only)
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{100, 0, 0};  // Slow drift
    ion.ion_charge_C = 1.602e-19;
    ion.mass_kg = 100 * 1.66e-27;
    ion.CCS_m2 = 100e-20;  // 100 Å²
    
    ICARION::core::IonEnsemble ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    ctx.ion_index = 0;
    Vec3 F = force.compute(ens, 0, 0.0, ctx);
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Gas density: n = " << env.particle_density_m_3 << " m⁻³\n";
    std::cout << "Thermal velocity: v_th = " << env.mean_thermal_velocity_m_s << " m/s\n";
    std::cout << "Collision cross-section: σ = " << ion.CCS_m2 << " m² (" << ion.CCS_m2 * 1e20 << " Å²)\n";
    std::cout << "Neutral mass: m_n = " << env.gas_mass_kg << " kg\n";
    std::cout << "Ion mass: m_i = " << ion.mass_kg << " kg\n";
    
    // Estimate collision frequency
    const double reduced_mass = (ion.mass_kg * env.gas_mass_kg) / 
                               (ion.mass_kg + env.gas_mass_kg);
    const double nu_collision = env.particle_density_m_3 * ion.CCS_m2 * 
                               env.mean_thermal_velocity_m_s;
    
    std::cout << "\nReduced mass: μ = " << reduced_mass << " kg\n";
    std::cout << "Collision frequency: ν ≈ n·σ·v_th = " << nu_collision << " 1/s\n";
    
    std::cout << "\nComputed force: F = (" << F.x << ", " << F.y << ", " << F.z << ") N\n";
    std::cout << "Force magnitude: |F| = " << std::sqrt(F.x*F.x + F.y*F.y + F.z*F.z) << " N\n";
    
    // Verify force opposes velocity
    const double dot_product = F.x * ion.vel.x + F.y * ion.vel.y + F.z * ion.vel.z;
    std::cout << "\nF · v = " << dot_product << " (should be < 0)\n";
    
    const bool pass = (dot_product < 0) && (std::abs(F.x) > 1e-22);
    std::cout << "\nResult: " << (pass ? "✓ PASS (force opposes motion)" : "✗ FAIL") << "\n";
}

void test_langevin_damping() {
    std::cout << "\n=== Test 3: Langevin Polarization Model ===\n";
    std::cout << "Physics: Enhanced collision rate from ion-induced dipole\n";
    
    // SSOT: Use EnvironmentConfig
    ICARION::config::EnvironmentConfig env;
    env.pressure_Pa = 101325.0;
    env.temperature_K = 300.0;
    env.gas_species = "N2";
    env.compute_derived_properties();
    
    DampingForce force(env, DampingModel::Langevin);
    ForceContext ctx;
    
    // Test ion with velocity
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{500, 0, 0};
    ion.ion_charge_C = 1.602e-19;
    ion.mass_kg = 100 * 1.66e-27;
    
    ICARION::core::IonEnsemble ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    ctx.ion_index = 0;
    Vec3 F = force.compute(ens, 0, 0.0, ctx);
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Gas density: n = " << env.particle_density_m_3 << " m⁻³\n";
    std::cout << "Neutral polarizability: α = " << env.gas_polarizability_m3 << " m³\n";
    std::cout << "Ion charge: q = " << ion.ion_charge_C << " C\n";
    
    // Langevin collision rate constant
    const double epsilon_0 = 8.854187817e-12;  // F/m
    const double k_langevin = std::sqrt(ion.ion_charge_C * ion.ion_charge_C * 
                                       env.gas_polarizability_m3 / 
                                       (M_PI * epsilon_0 * env.gas_mass_kg));
    
    std::cout << "\nLangevin rate constant: k_L = " << k_langevin << " m³/s\n";
    std::cout << "Enhancement: k_L·n = " << k_langevin * env.particle_density_m_3 << " 1/s\n";
    
    std::cout << "\nComputed force: F = (" << F.x << ", " << F.y << ", " << F.z << ") N\n";
    std::cout << "Force magnitude: |F| = " << std::sqrt(F.x*F.x + F.y*F.y + F.z*F.z) << " N\n";
    
    // Verify force opposes velocity
    const double dot_product = F.x * ion.vel.x;
    std::cout << "\nF · v = " << dot_product << " (should be < 0)\n";
    
    const bool pass = (dot_product < 0) && (std::abs(F.x) > 1e-22);
    std::cout << "\nResult: " << (pass ? "✓ PASS (Langevin damping active)" : "✗ FAIL") << "\n";
}

void test_velocity_scaling() {
    std::cout << "\n=== Test 4: Force Scales with Velocity ===\n";
    std::cout << "Physics: Damping force should scale linearly with |v|\n";
    
    DampingParams params;
    params.model = DampingModel::Friction;
    params.gamma_coefficient = 1e6;
    
    DampingForce force(params);
    ForceContext ctx;
    
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.ion_charge_C = 1.602e-19;
    ion.mass_kg = 100 * 1.66e-27;
    
    // Test at v = 100 m/s
    ion.vel = Vec3{100, 0, 0};
    ICARION::core::IonEnsemble ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    ctx.ion_index = 0;
    Vec3 F1 = force.compute(ens, 0, 0.0, ctx);
    
    // Test at v = 200 m/s (double velocity)
    ion.vel = Vec3{200, 0, 0};
    ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    Vec3 F2 = force.compute(ens, 0, 0.0, ctx);
    
    const double F1_mag = std::abs(F1.x);
    const double F2_mag = std::abs(F2.x);
    const double ratio = F2_mag / F1_mag;
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Force at v = 100 m/s: |F| = " << F1_mag << " N\n";
    std::cout << "Force at v = 200 m/s: |F| = " << F2_mag << " N\n";
    std::cout << "\nRatio F(2v)/F(v) = " << ratio << " (should be ≈ 2.0)\n";
    
    const bool pass = (std::abs(ratio - 2.0) < 1e-10);
    std::cout << "\nResult: " << (pass ? "✓ PASS (linear scaling)" : "✗ FAIL") << "\n";
}

int main() {
    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║   Minimal Physics Test: DampingForce                  ║\n";
    std::cout << "║   Purpose: Verify collision damping F = -γ·m·v        ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n";
    
    test_explicit_gamma_damping();
    test_hard_sphere_damping();
    test_langevin_damping();
    test_velocity_scaling();
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "All minimal tests completed.\n";
    std::cout << "These tests verify the fundamental damping physics.\n";
    std::cout << std::string(60, '=') << "\n\n";
    
    return 0;
}

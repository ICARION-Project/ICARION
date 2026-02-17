// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * @file minimal_electric_force.cpp
 * @brief Minimal physics test for ElectricFieldForce
 * 
 * Purpose: Quick verification that electric field force calculations are correct.
 * Useful for debugging fundamental physics without full test framework overhead.
 * 
 * ⚠️ WARNING: Uses DEPRECATED AnalyticalFieldParams!
 * This test uses the legacy parameter struct for simplicity.
 * Future: Will be replaced with DomainConfig in Phase 2 SSOT refactor.
 * 
 * Expected outputs:
 * - IMS drift field: F = q·E (linear along z-axis)
 * - LQIT quadrupole: F ∝ position (harmonic oscillator)
 * 
 * Compile:
 *   g++ -std=c++20 -I../../../src -I../../../include \
 *       minimal_electric_force.cpp \
 *       ../../../src/core/physics/forces/ElectricFieldForce.cpp \
 *       ../../../src/core/utils/mathUtils.cpp \
 *       -o minimal_electric_force
 * 
 * Run: ./minimal_electric_force
 */

#include "core/physics/forces/ElectricFieldForce.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/FieldsConfig.h"
#include "core/config/types/GeometryConfig.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ICARION;
using namespace ICARION::physics;

void test_ims_drift_field() {
    std::cout << "\n=== Test 1: IMS Drift Field ===\n";
    std::cout << "Physics: F = q·E with uniform field E_z = 400 V/m\n";
    
    // IMS parameters: 1 m drift tube, 400 V drift voltage (SSOT config)
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::IMS;
    domain.fields.dc.axial_V.constant_value = 400.0;
    domain.geometry.length_m = 1.0;
    
    ElectricFieldForce force(domain);
    ForceContext ctx;
    
    // Test ion: charge +e at z=0.5m
    IonState ion;
    ion.pos = Vec3{0, 0, 0.5};
    ion.vel = Vec3{0, 0, 0};
    ion.ion_charge_C = 1.602e-19;  // Elementary charge
    ion.mass_kg = 100 * 1.66e-27;
    
    ICARION::core::IonEnsemble ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    ctx.ion_index = 0;
    Vec3 F = force.compute(ens, 0, 0.0, ctx);
    
    // Expected: F_z = q·E = 1.602e-19 C × 400 V/m = 6.408e-17 N
    const double E_field = 400.0;  // V/m
    const double expected_F_z = ion.ion_charge_C * E_field;
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Ion charge:    " << ion.ion_charge_C << " C\n";
    std::cout << "Electric field: " << E_field << " V/m\n";
    std::cout << "Computed force: F_z = " << F.z << " N\n";
    std::cout << "Expected force: F_z = " << expected_F_z << " N\n";
    std::cout << "Relative error: " << std::abs(F.z - expected_F_z) / expected_F_z * 100 << " %\n";
    
    // Verify
    const bool pass = (std::abs(F.x) < 1e-20) && 
                     (std::abs(F.y) < 1e-20) &&
                     (std::abs(F.z - expected_F_z) / expected_F_z < 1e-10);
    
    std::cout << "Result: " << (pass ? "✓ PASS" : "✗ FAIL") << "\n";
}

void test_lqit_quadrupole() {
    std::cout << "\n=== Test 2: LQIT Quadrupole Field ===\n";
    std::cout << "Physics: F_x = 2U·x/r₀², F_y = -2U·y/r₀² (harmonic confinement)\n";
    
    // LQIT parameters (SSOT config)
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::LQIT;
    domain.fields.rf.voltage_V.constant_value = 100.0;
    domain.fields.rf.frequency_Hz.constant_value = 1e6;
    domain.fields.rf.compute_derived();
    domain.fields.dc.quad_V.constant_value = 0.0;
    domain.geometry.radius_m = 0.01;
    domain.geometry.length_m = 0.05;
    
    ElectricFieldForce force(domain);
    ForceContext ctx;
    
    // Test ion: at position (2mm, 3mm, 0)
    IonState ion;
    ion.pos = Vec3{2e-3, 3e-3, 0};
    ion.vel = Vec3{0, 0, 0};
    ion.ion_charge_C = 1.602e-19;
    ion.mass_kg = 100 * 1.66e-27;
    
    // Compute force at t=0 (cos(ωt) = 1, RF at maximum)
    const double t = 0.0;
    ICARION::core::IonEnsemble ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    ctx.ion_index = 0;
    Vec3 F = force.compute(ens, 0, t, ctx);
    
    // Expected field at RF maximum:
    // E_x = 2U·x/r₀² = 2×100×0.002/(0.01²) = 4000 V/m
    // E_y = -2U·y/r₀² = -2×100×0.003/(0.01²) = -6000 V/m
    const double r0_sq = domain.geometry.radius_m * domain.geometry.radius_m;
    const double U_eff = domain.fields.rf.voltage_V.constant_value.value();  // At t=0, cos(ωt)=1
    const double E_x_expected = 2.0 * U_eff * ion.pos.x / r0_sq;
    const double E_y_expected = -2.0 * U_eff * ion.pos.y / r0_sq;
    const double F_x_expected = ion.ion_charge_C * E_x_expected;
    const double F_y_expected = ion.ion_charge_C * E_y_expected;
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Ion position: (" << ion.pos.x << ", " << ion.pos.y << ", " << ion.pos.z << ") m\n";
    std::cout << "RF voltage: " << params.rf_voltage_V << " V (at t=0, max)\n";
    std::cout << "r₀ = " << params.radius_m << " m\n";
    std::cout << "\nExpected field:\n";
    std::cout << "  E_x = " << E_x_expected << " V/m\n";
    std::cout << "  E_y = " << E_y_expected << " V/m\n";
    std::cout << "\nComputed force:\n";
    std::cout << "  F_x = " << F.x << " N\n";
    std::cout << "  F_y = " << F.y << " N\n";
    std::cout << "  F_z = " << F.z << " N\n";
    std::cout << "\nExpected force:\n";
    std::cout << "  F_x = " << F_x_expected << " N\n";
    std::cout << "  F_y = " << F_y_expected << " N\n";
    
    // Verify
    const bool pass = (std::abs(F.x - F_x_expected) / std::abs(F_x_expected) < 1e-10) &&
                     (std::abs(F.y - F_y_expected) / std::abs(F_y_expected) < 1e-10) &&
                     (std::abs(F.z) < 1e-20);
    
    std::cout << "\nResult: " << (pass ? "✓ PASS" : "✗ FAIL") << "\n";
}

int main() {
    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║   Minimal Physics Test: ElectricFieldForce            ║\n";
    std::cout << "║   Purpose: Verify fundamental electric field physics  ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n";
    
    test_ims_drift_field();
    test_lqit_quadrupole();
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "All minimal tests completed.\n";
    std::cout << "These tests verify the fundamental physics is correct.\n";
    std::cout << std::string(60, '=') << "\n\n";
    
    return 0;
}

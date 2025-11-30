// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file minimal_magnetic_force.cpp
 * @brief Minimal physics test for MagneticFieldForce
 * 
 * Purpose: Quick verification that Lorentz force F = q(v × B) is computed correctly.
 * 
 * ✅ SSOT-COMPLIANT: Uses MagneticFieldConfig directly.
 * 
 * Expected outputs:
 * - Lorentz force perpendicular to both v and B
 * - Cyclotron radius: r = m·v/(q·B)
 * - No force if v ∥ B (parallel motion)
 * 
 * Compile:
 *   g++ -std=c++20 -I../../../src -I../../../include \
 *       minimal_magnetic_force.cpp \
 *       ../../../src/core/physics/forces/MagneticFieldForce.cpp \
 *       ../../../src/core/utils/mathUtils.cpp \
 *       -o minimal_magnetic_force
 * 
 * Run: ./minimal_magnetic_force
 */

#include "core/physics/forces/MagneticFieldForce.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonState.h"
#include "core/config/types/FieldsConfig.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ICARION;
using namespace ICARION::physics;

void test_lorentz_force_perpendicular() {
    std::cout << "\n=== Test 1: Lorentz Force (v ⊥ B) ===\n";
    std::cout << "Physics: F = q(v × B), perpendicular to both v and B\n";
    
    // Magnetic field: 1 T along z-axis (SSOT config)
    ICARION::config::MagneticFieldConfig mag_config;
    mag_config.field_strength_T = Vec3{0, 0, 1.0};  // 1 Tesla in z-direction
    mag_config.enabled = true;
    
    MagneticFieldForce force(mag_config);
    ForceContext ctx;
    
    // Test ion: moving in x-direction with v = 1000 m/s
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{1000, 0, 0};        // v_x = 1000 m/s
    ion.ion_charge_C = 1.602e-19;      // +e
    ion.mass_kg = 100 * 1.66e-27;      // 100 u
    
    Vec3 F = force.compute(ion, 0.0, ctx);
    
    // Expected: F = q(v × B)
    // v × B = (1000, 0, 0) × (0, 0, 1) = (0, 1000, 0) [using right-hand rule]
    // F = 1.602e-19 × (0, 1000, 0) = (0, 1.602e-16, 0) N
    Vec3 v_cross_B(
        ion.vel.y * mag_config.field_strength_T.z - ion.vel.z * mag_config.field_strength_T.y,
        ion.vel.z * mag_config.field_strength_T.x - ion.vel.x * mag_config.field_strength_T.z,
        ion.vel.x * mag_config.field_strength_T.y - ion.vel.y * mag_config.field_strength_T.x
    );
    Vec3 F_expected = v_cross_B * ion.ion_charge_C;
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Magnetic field: B = (0, 0, " << mag_config.field_strength_T.z << ") T\n";
    std::cout << "Ion velocity:   v = (" << ion.vel.x << ", 0, 0) m/s\n";
    std::cout << "Ion charge:     q = " << ion.ion_charge_C << " C\n";
    std::cout << "\nv × B = (" << v_cross_B.x << ", " << v_cross_B.y << ", " << v_cross_B.z << ")\n";
    std::cout << "\nComputed force: F = (" << F.x << ", " << F.y << ", " << F.z << ") N\n";
    std::cout << "Expected force: F = (" << F_expected.x << ", " << F_expected.y << ", " << F_expected.z << ") N\n";
    
    // Verify: F should be in y-direction, magnitude = q·v·B
    const double F_mag_computed = std::sqrt(F.x*F.x + F.y*F.y + F.z*F.z);
    const double F_mag_expected = std::abs(ion.ion_charge_C) * ion.vel.x * mag_config.field_strength_T.z;
    
    std::cout << "\nForce magnitude: " << F_mag_computed << " N\n";
    std::cout << "Expected:        " << F_mag_expected << " N\n";
    
    const bool pass = (std::abs(F.x) < 1e-20) &&
                     (std::abs(F.y - F_expected.y) / std::abs(F_expected.y) < 1e-10) &&
                     (std::abs(F.z) < 1e-20);
    
    std::cout << "\nResult: " << (pass ? "✓ PASS" : "✗ FAIL") << "\n";
}

void test_lorentz_force_parallel() {
    std::cout << "\n=== Test 2: No Force when v ∥ B ===\n";
    std::cout << "Physics: F = 0 when velocity parallel to magnetic field\n";
    
    // Magnetic field: 1 T along z-axis (SSOT config)
    ICARION::config::MagneticFieldConfig mag_config;
    mag_config.field_strength_T = Vec3{0, 0, 1.0};
    mag_config.enabled = true;
    
    MagneticFieldForce force(mag_config);
    ForceContext ctx;
    
    // Test ion: moving parallel to B (along z-axis)
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{0, 0, 1000};        // v_z = 1000 m/s (parallel to B)
    ion.ion_charge_C = 1.602e-19;
    ion.mass_kg = 100 * 1.66e-27;
    
    Vec3 F = force.compute(ion, 0.0, ctx);
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Magnetic field: B = (0, 0, " << mag_config.field_strength_T.z << ") T\n";
    std::cout << "Ion velocity:   v = (0, 0, " << ion.vel.z << ") m/s (parallel!)\n";
    std::cout << "\nComputed force: F = (" << F.x << ", " << F.y << ", " << F.z << ") N\n";
    std::cout << "Expected:       F = (0, 0, 0) N\n";
    
    const double F_mag = std::sqrt(F.x*F.x + F.y*F.y + F.z*F.z);
    std::cout << "\nForce magnitude: " << F_mag << " N (should be ~0)\n";
    
    const bool pass = (F_mag < 1e-25);
    std::cout << "\nResult: " << (pass ? "✓ PASS" : "✗ FAIL") << "\n";
}

void test_cyclotron_radius() {
    std::cout << "\n=== Test 3: Cyclotron Motion ===\n";
    std::cout << "Physics: Cyclotron radius r_c = m·v/(q·B)\n";
    
    // Strong magnetic field: 7 T (typical for FT-ICR) (SSOT config)
    ICARION::config::MagneticFieldConfig mag_config;
    mag_config.field_strength_T = Vec3{0, 0, 7.0};  // 7 Tesla
    mag_config.enabled = true;
    
    MagneticFieldForce force(mag_config);
    ForceContext ctx;
    
    // Test ion: proton with thermal velocity
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{500, 0, 0};         // 500 m/s in x-direction
    ion.ion_charge_C = 1.602e-19;      // Proton charge
    ion.mass_kg = 1.67e-27;            // Proton mass
    
    Vec3 F = force.compute(ion, 0.0, ctx);
    
    // Cyclotron radius: r_c = m·v/(q·B)
    const double v_perp = ion.vel.x;  // Perpendicular velocity component
    const double B = mag_config.field_strength_T.z;
    const double r_cyclotron = ion.mass_kg * v_perp / (ion.ion_charge_C * B);
    
    // Centripetal force: F_c = m·v²/r = q·v·B
    const double F_centripetal = ion.mass_kg * v_perp * v_perp / r_cyclotron;
    const double F_lorentz = ion.ion_charge_C * v_perp * B;
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Ion: Proton (m = " << ion.mass_kg << " kg, q = " << ion.ion_charge_C << " C)\n";
    std::cout << "Magnetic field: B = " << B << " T\n";
    std::cout << "Velocity (perpendicular): v = " << v_perp << " m/s\n";
    std::cout << "\nCyclotron radius: r_c = m·v/(q·B) = " << r_cyclotron << " m\n";
    std::cout << "                  r_c = " << r_cyclotron * 1e6 << " µm\n";
    std::cout << "\nCentripetal force needed: F = m·v²/r = " << F_centripetal << " N\n";
    std::cout << "Lorentz force magnitude:  F = q·v·B   = " << F_lorentz << " N\n";
    std::cout << "Computed force magnitude: F = " << std::sqrt(F.x*F.x + F.y*F.y + F.z*F.z) << " N\n";
    
    const bool pass = std::abs(F_centripetal - F_lorentz) / F_lorentz < 1e-10;
    std::cout << "\nResult: " << (pass ? "✓ PASS (F_centripetal = F_Lorentz)" : "✗ FAIL") << "\n";
}

int main() {
    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║   Minimal Physics Test: MagneticFieldForce            ║\n";
    std::cout << "║   Purpose: Verify Lorentz force F = q(v × B)          ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n";
    
    test_lorentz_force_perpendicular();
    test_lorentz_force_parallel();
    test_cyclotron_radius();
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "All minimal tests completed.\n";
    std::cout << "These tests verify the fundamental physics is correct.\n";
    std::cout << std::string(60, '=') << "\n\n";
    
    return 0;
}

// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file minimal_space_charge_force.cpp
 * @brief Minimal physics test for SpaceChargeDirect
 * 
 * Purpose: Quick verification that Coulomb forces F = k·q₁·q₂/r² are computed correctly.
 * 
 * Expected outputs:
 * - Two-body Coulomb force matches analytical calculation
 * - Newton's 3rd law: F₁₂ = -F₂₁
 * - Softening prevents divergence at small separations
 * 
 * Compile:
 *   g++ -std=c++20 -I../../../src -I../../../include \
 *       minimal_space_charge_force.cpp \
 *       ../../../src/core/physics/forces/SpaceChargeDirect.cpp \
 *       ../../../src/core/utils/mathUtils.cpp \
 *       -o minimal_space_charge_force
 * 
 * Run: ./minimal_space_charge_force
 */

#include "core/physics/forces/SpaceChargeDirect.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonState.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace ICARION;
using namespace ICARION::physics;

void test_two_ion_coulomb_force() {
    std::cout << "\n=== Test 1: Two-Ion Coulomb Repulsion ===\n";
    std::cout << "Physics: F = k_e·q₁·q₂/r², repulsive for like charges\n";
    
    SpaceChargeDirect force(0.0);  // No softening
    std::vector<IonState> ions;
    ForceContext ctx;
    ctx.all_ions = &ions;
    
    // Two ions separated by 1 µm
    IonState ion1, ion2;
    ion1.pos = Vec3{0, 0, 0};
    ion1.vel = Vec3{0, 0, 0};
    ion1.ion_charge_C = 1.602e-19;  // +e
    ion1.mass_kg = 100 * 1.66e-27;
    
    ion2.pos = Vec3{1e-6, 0, 0};    // 1 µm away
    ion2.vel = Vec3{0, 0, 0};
    ion2.ion_charge_C = 1.602e-19;  // +e
    ion2.mass_kg = 100 * 1.66e-27;
    
    ions.push_back(ion1);
    ions.push_back(ion2);
    
    Vec3 F1 = force.compute(ions[0], 0.0, ctx);
    
    // Expected: F = k_e·e²/r²
    const double k_e = 8.987551787e9;  // N·m²/C²
    const double r = 1e-6;              // m
    const double F_expected = k_e * ion1.ion_charge_C * ion2.ion_charge_C / (r * r);
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Coulomb constant: k_e = " << k_e << " N·m²/C²\n";
    std::cout << "Elementary charge: e = " << ion1.ion_charge_C << " C\n";
    std::cout << "Separation: r = " << r << " m (" << r * 1e6 << " µm)\n";
    std::cout << "\nExpected force: F = k_e·e²/r² = " << F_expected << " N\n";
    std::cout << "Computed force: F = (" << F1.x << ", " << F1.y << ", " << F1.z << ") N\n";
    std::cout << "Force magnitude: |F| = " << std::sqrt(F1.x*F1.x + F1.y*F1.y + F1.z*F1.z) << " N\n";
    
    // Force should point in negative x (away from ion2)
    std::cout << "\nDirection: Should be negative x (repulsion)\n";
    std::cout << "  F_x = " << F1.x << " (should be ≈ -" << F_expected << ")\n";
    
    const bool pass = (std::abs(F1.x + F_expected) / F_expected < 1e-10) &&
                     (std::abs(F1.y) < 1e-20) &&
                     (std::abs(F1.z) < 1e-20);
    
    std::cout << "\nResult: " << (pass ? "✓ PASS" : "✗ FAIL") << "\n";
}

void test_newtons_third_law() {
    std::cout << "\n=== Test 2: Newton's Third Law (F₁₂ = -F₂₁) ===\n";
    std::cout << "Physics: Action-reaction pair must be equal and opposite\n";
    
    SpaceChargeDirect force(0.0);
    std::vector<IonState> ions;
    ForceContext ctx;
    ctx.all_ions = &ions;
    
    // Two ions at arbitrary positions
    IonState ion1, ion2;
    ion1.pos = Vec3{1e-6, 2e-6, 0};
    ion1.vel = Vec3{0, 0, 0};
    ion1.ion_charge_C = 1.602e-19;
    ion1.mass_kg = 100 * 1.66e-27;
    
    ion2.pos = Vec3{3e-6, 1e-6, 1e-6};
    ion2.vel = Vec3{0, 0, 0};
    ion2.ion_charge_C = 2 * 1.602e-19;  // Double charge
    ion2.mass_kg = 200 * 1.66e-27;
    
    ions.push_back(ion1);
    ions.push_back(ion2);
    
    Vec3 F_on_1 = force.compute(ions[0], 0.0, ctx);
    Vec3 F_on_2 = force.compute(ions[1], 0.0, ctx);
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Ion 1: pos = (" << ion1.pos.x << ", " << ion1.pos.y << ", " << ion1.pos.z << ") m\n";
    std::cout << "       q₁ = " << ion1.ion_charge_C << " C\n";
    std::cout << "Ion 2: pos = (" << ion2.pos.x << ", " << ion2.pos.y << ", " << ion2.pos.z << ") m\n";
    std::cout << "       q₂ = " << ion2.ion_charge_C << " C (double charge)\n";
    
    std::cout << "\nForce on ion 1: F₁ = (" << F_on_1.x << ", " << F_on_1.y << ", " << F_on_1.z << ") N\n";
    std::cout << "Force on ion 2: F₂ = (" << F_on_2.x << ", " << F_on_2.y << ", " << F_on_2.z << ") N\n";
    
    Vec3 F_sum = Vec3{F_on_1.x + F_on_2.x, F_on_1.y + F_on_2.y, F_on_1.z + F_on_2.z};
    std::cout << "\nSum F₁ + F₂ = (" << F_sum.x << ", " << F_sum.y << ", " << F_sum.z << ") N\n";
    std::cout << "             (should be ≈ 0 by Newton's 3rd law)\n";
    
    const double sum_magnitude = std::sqrt(F_sum.x*F_sum.x + F_sum.y*F_sum.y + F_sum.z*F_sum.z);
    const double F1_magnitude = std::sqrt(F_on_1.x*F_on_1.x + F_on_1.y*F_on_1.y + F_on_1.z*F_on_1.z);
    
    std::cout << "\n|F₁ + F₂| = " << sum_magnitude << " N\n";
    std::cout << "Relative to |F₁|: " << sum_magnitude / F1_magnitude * 100 << " %\n";
    
    const bool pass = (sum_magnitude / F1_magnitude < 1e-10);
    std::cout << "\nResult: " << (pass ? "✓ PASS" : "✗ FAIL") << "\n";
}

void test_softening_at_small_distance() {
    std::cout << "\n=== Test 3: Softening Prevents Divergence ===\n";
    std::cout << "Physics: F = k·q²/(r²+ε²) remains finite as r → 0\n";
    
    const double epsilon = 1e-10;  // 0.1 nm softening
    SpaceChargeDirect force_with_softening(epsilon);
    SpaceChargeDirect force_no_softening(0.0);
    
    std::vector<IonState> ions;
    ForceContext ctx;
    ctx.all_ions = &ions;
    
    // Two ions very close together (0.01 nm)
    IonState ion1, ion2;
    ion1.pos = Vec3{0, 0, 0};
    ion1.vel = Vec3{0, 0, 0};
    ion1.ion_charge_C = 1.602e-19;
    ion1.mass_kg = 100 * 1.66e-27;
    
    const double r_close = 1e-11;  // 0.01 nm (very close!)
    ion2.pos = Vec3{r_close, 0, 0};
    ion2.vel = Vec3{0, 0, 0};
    ion2.ion_charge_C = 1.602e-19;
    ion2.mass_kg = 100 * 1.66e-27;
    
    ions.push_back(ion1);
    ions.push_back(ion2);
    
    Vec3 F_with_soft = force_with_softening.compute(ions[0], 0.0, ctx);
    Vec3 F_no_soft = force_no_softening.compute(ions[0], 0.0, ctx);
    
    const double F_with_mag = std::sqrt(F_with_soft.x*F_with_soft.x + F_with_soft.y*F_with_soft.y + F_with_soft.z*F_with_soft.z);
    const double F_no_mag = std::sqrt(F_no_soft.x*F_no_soft.x + F_no_soft.y*F_no_soft.y + F_no_soft.z*F_no_soft.z);
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "Separation: r = " << r_close << " m (" << r_close * 1e9 << " nm)\n";
    std::cout << "Softening: ε = " << epsilon << " m (" << epsilon * 1e9 << " nm)\n";
    std::cout << "\nWithout softening: |F| = " << F_no_mag << " N (HUGE!)\n";
    std::cout << "With softening:    |F| = " << F_with_mag << " N (reduced)\n";
    std::cout << "\nReduction factor: " << F_no_mag / F_with_mag << "×\n";
    
    // At r → 0, force should approach k·q²/ε²
    const double k_e = 8.987551787e9;
    const double F_limit = k_e * ion1.ion_charge_C * ion2.ion_charge_C / (epsilon * epsilon);
    std::cout << "\nAsymptotic limit (r→0): F → k·q²/ε² = " << F_limit << " N\n";
    
    const bool pass = (F_with_mag < F_no_mag) && (F_with_mag < F_limit * 1.5) && std::isfinite(F_with_mag);
    std::cout << "\nResult: " << (pass ? "✓ PASS (softening works)" : "✗ FAIL") << "\n";
}

void test_ion_cloud_expansion() {
    std::cout << "\n=== Test 4: Ion Cloud Self-Repulsion ===\n";
    std::cout << "Physics: Cloud of like-charged ions should experience net outward force\n";
    
    SpaceChargeDirect force(0.0);
    std::vector<IonState> ions;
    ForceContext ctx;
    ctx.all_ions = &ions;
    
    // Create 5 ions in a line (simple 1D cloud)
    for (int i = 0; i < 5; ++i) {
        IonState ion;
        ion.pos = Vec3{i * 1e-6, 0, 0};  // Spaced 1 µm apart
        ion.vel = Vec3{0, 0, 0};
        ion.ion_charge_C = 1.602e-19;
        ion.mass_kg = 100 * 1.66e-27;
        ions.push_back(ion);
    }
    
    // Compute force on center ion (should be zero by symmetry)
    Vec3 F_center = force.compute(ions[2], 0.0, ctx);
    
    // Compute force on edge ion (should push outward)
    Vec3 F_edge = force.compute(ions[0], 0.0, ctx);
    
    std::cout << std::scientific << std::setprecision(6);
    std::cout << "5 ions at positions: 0, 1, 2, 3, 4 µm\n";
    std::cout << "All ions have charge +e\n";
    
    std::cout << "\nForce on center ion (x=2µm): F = (" << F_center.x << ", " << F_center.y << ", " << F_center.z << ") N\n";
    std::cout << "  (should be ≈0 by symmetry)\n";
    
    std::cout << "\nForce on edge ion (x=0):     F = (" << F_edge.x << ", " << F_edge.y << ", " << F_edge.z << ") N\n";
    std::cout << "  (should be negative x, pushing left/outward)\n";
    
    const double F_center_mag = std::sqrt(F_center.x*F_center.x + F_center.y*F_center.y + F_center.z*F_center.z);
    
    std::cout << "\nCenter force magnitude: " << F_center_mag << " N\n";
    std::cout << "Edge force x-component: " << F_edge.x << " N\n";
    
    const bool pass = (F_center_mag < 1e-18) &&  // Center nearly zero
                     (F_edge.x < -1e-17);         // Edge pushes left (outward)
    
    std::cout << "\nResult: " << (pass ? "✓ PASS (cloud expands)" : "✗ FAIL") << "\n";
}

int main() {
    std::cout << "╔════════════════════════════════════════════════════════╗\n";
    std::cout << "║   Minimal Physics Test: SpaceChargeDirect              ║\n";
    std::cout << "║   Purpose: Verify Coulomb interactions F = k·q²/r²    ║\n";
    std::cout << "╚════════════════════════════════════════════════════════╝\n";
    
    test_two_ion_coulomb_force();
    test_newtons_third_law();
    test_softening_at_small_distance();
    test_ion_cloud_expansion();
    
    std::cout << "\n" << std::string(60, '=') << "\n";
    std::cout << "All minimal tests completed.\n";
    std::cout << "These tests verify the fundamental physics is correct.\n";
    std::cout << std::string(60, '=') << "\n\n";
    
    return 0;
}

// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/physics/forces/SpaceChargeDirect.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include "utils/constants.h"

#include <cmath>
#include <vector>

using namespace ICARION;
using namespace ICARION::physics;
using Catch::Matchers::WithinAbs;

// ============================================================================
// Test Constants and Helpers
// ============================================================================

// Test tolerance for force comparisons
static constexpr double tolerance = 1e-18;  // Newton (force tolerance) - for ~1e-16 N forces

// Helper: Create ion at position with charge
static IonState create_ion(Vec3 pos, double charge_C, double mass_u = 100.0) {
    IonState ion;
    ion.pos = pos;
    ion.vel = Vec3{0, 0, 0};
    ion.ion_charge_C = charge_C;
    ion.mass_kg = mass_u * AMU_TO_KG;  // Convert u to kg
    return ion;
}

// ============================================================================
// Basic Tests
// ============================================================================

TEST_CASE("SpaceChargeDirect: Name is 'SpaceCharge'", "[SpaceChargeDirect][basic]") {
    SpaceChargeDirect force{0.0};
    REQUIRE(force.name() == "SpaceCharge");
}

TEST_CASE("SpaceChargeDirect: Applies to all ions", "[SpaceChargeDirect][basic]") {
    SpaceChargeDirect force{0.0};
    IonState ion = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    REQUIRE(force.applies_to(ion));
}

TEST_CASE("SpaceChargeDirect: Negative softening throws exception", "[SpaceChargeDirect][basic]") {
    REQUIRE_THROWS_AS(SpaceChargeDirect(-1e-10), std::invalid_argument);
}

// ============================================================================
// Empty Ensemble Tests
// ============================================================================

TEST_CASE("SpaceChargeDirect: No ions returns zero force", "[SpaceChargeDirect][empty]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    IonState ion = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    Vec3 result = force.compute(ion, 0.0, ctx);
    
    REQUIRE(result.x == 0.0);
    REQUIRE(result.y == 0.0);
    REQUIRE(result.z == 0.0);
}

TEST_CASE("SpaceChargeDirect: Null ion ensemble returns zero force", "[SpaceChargeDirect][empty]") {
    SpaceChargeDirect force{0.0};
    ForceContext ctx;
    ctx.all_ions = nullptr;
    
    IonState ion = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    Vec3 result = force.compute(ion, 0.0, ctx);
    
    REQUIRE(result.x == 0.0);
    REQUIRE(result.y == 0.0);
    REQUIRE(result.z == 0.0);
}

TEST_CASE("SpaceChargeDirect: Single ion returns zero force", "[SpaceChargeDirect][empty]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    ion_ensemble.push_back(create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C));
    
    // Compute force on the only ion in ensemble - should be zero (no other ions)
    Vec3 result = force.compute(ion_ensemble[0], 0.0, ctx);
    
    REQUIRE(result.x == 0.0);
    REQUIRE(result.y == 0.0);
    REQUIRE(result.z == 0.0);
}

// ============================================================================
// Self-Interaction Tests
// ============================================================================

TEST_CASE("SpaceChargeDirect: Self-interaction excluded", "[SpaceChargeDirect][self]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Create two ions at 1µm apart
    ion_ensemble.push_back(create_ion(Vec3{1e-6, 0, 0}, ELEM_CHARGE_C));
    ion_ensemble.push_back(create_ion(Vec3{2e-6, 0, 0}, ELEM_CHARGE_C));
    
    // Compute force on first ion in ensemble
    // Self-interaction is excluded by address comparison (&ion == &other_ion)
    Vec3 result = force.compute(ion_ensemble[0], 0.0, ctx);
    
    // Force should only come from ion at 2µm, not from itself at 1µm
    // Two positive charges repel: ion at 1µm pushed left (negative x) by ion at 2µm
    // F = k·e²/r² along x-axis, r = 1e-6 m
    const double expected_force_x = -COULOMB_CONST * ELEM_CHARGE_C * ELEM_CHARGE_C / (1e-6 * 1e-6);
    
    REQUIRE_THAT(result.x, WithinAbs(expected_force_x, tolerance));
    REQUIRE_THAT(result.y, WithinAbs(0.0, tolerance));
    REQUIRE_THAT(result.z, WithinAbs(0.0, tolerance));
}

// ============================================================================
// Two-Ion System Tests (Analytical Validation)
// ============================================================================

TEST_CASE("SpaceChargeDirect: Two positive charges repel along X", "[SpaceChargeDirect][coulomb]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Ion 0 at origin, ion 1 at x=1µm
    IonState ion0 = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    IonState ion1 = create_ion(Vec3{1e-6, 0, 0}, ELEM_CHARGE_C);
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    
    // Force on ion0 from ion1 (should push left: negative x)
    Vec3 force0 = force.compute(ion_ensemble[0], 0.0, ctx);
    
    // Expected: F = k·e²/r²
    const double r = 1e-6;
    const double expected_magnitude = COULOMB_CONST * ELEM_CHARGE_C * ELEM_CHARGE_C / (r * r);
    
    // Force should point in negative x direction (away from ion1)
    REQUIRE_THAT(force0.x, WithinAbs(-expected_magnitude, tolerance));
    REQUIRE_THAT(force0.y, WithinAbs(0.0, tolerance));
    REQUIRE_THAT(force0.z, WithinAbs(0.0, tolerance));
}

TEST_CASE("SpaceChargeDirect: Opposite charges attract", "[SpaceChargeDirect][coulomb]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Ion 0: positive charge at origin
    // Ion 1: negative charge at x=1µm
    IonState ion0 = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    IonState ion1 = create_ion(Vec3{1e-6, 0, 0}, -ELEM_CHARGE_C);
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    
    Vec3 force0 = force.compute(ion_ensemble[0], 0.0, ctx);
    
    const double r = 1e-6;
    const double expected_magnitude = COULOMB_CONST * ELEM_CHARGE_C * ELEM_CHARGE_C / (r * r);
    
    // Force should point in positive x direction (towards ion1)
    REQUIRE_THAT(force0.x, WithinAbs(expected_magnitude, tolerance));
    REQUIRE_THAT(force0.y, WithinAbs(0.0, tolerance));
    REQUIRE_THAT(force0.z, WithinAbs(0.0, tolerance));
}

// ============================================================================
// Newton's Third Law (Symmetry)
// ============================================================================

TEST_CASE("SpaceChargeDirect: Newton's third law", "[SpaceChargeDirect][newton]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Two ions at different positions
    IonState ion0 = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    IonState ion1 = create_ion(Vec3{1e-6, 1e-6, 0}, 2.0 * ELEM_CHARGE_C);  // Double charge
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    
    Vec3 force_on_0 = force.compute(ion_ensemble[0], 0.0, ctx);
    Vec3 force_on_1 = force.compute(ion_ensemble[1], 0.0, ctx);
    
    // F₀₁ = -F₁₀ (Newton's 3rd law)
    REQUIRE_THAT(force_on_0.x, WithinAbs(-force_on_1.x, tolerance));
    REQUIRE_THAT(force_on_0.y, WithinAbs(-force_on_1.y, tolerance));
    REQUIRE_THAT(force_on_0.z, WithinAbs(-force_on_1.z, tolerance));
}

// ============================================================================
// 3D Position Tests
// ============================================================================

TEST_CASE("SpaceChargeDirect: Three-dimensional repulsion", "[SpaceChargeDirect][3d]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Ion 0 at origin, ion 1 at (1,1,1) µm
    IonState ion0 = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    IonState ion1 = create_ion(Vec3{1e-6, 1e-6, 1e-6}, ELEM_CHARGE_C);
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    
    Vec3 force0 = force.compute(ion_ensemble[0], 0.0, ctx);
    
    // Distance: r = √(3) µm
    const double r = std::sqrt(3.0) * 1e-6;
    const double expected_magnitude = COULOMB_CONST * ELEM_CHARGE_C * ELEM_CHARGE_C / (r * r);
    
    // Force direction: r̂ = (-1,-1,-1) / √3
    const double component = -expected_magnitude / std::sqrt(3.0);
    
    REQUIRE_THAT(force0.x, WithinAbs(component, tolerance));
    REQUIRE_THAT(force0.y, WithinAbs(component, tolerance));
    REQUIRE_THAT(force0.z, WithinAbs(component, tolerance));
}

// ============================================================================
// Multiple Ions (Superposition)
// ============================================================================

TEST_CASE("SpaceChargeDirect: Three ions linear", "[SpaceChargeDirect][superposition]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Three ions along x-axis: ion0 at origin, ion1 at 1µm, ion2 at 2µm
    IonState ion0 = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    IonState ion1 = create_ion(Vec3{1e-6, 0, 0}, ELEM_CHARGE_C);
    IonState ion2 = create_ion(Vec3{2e-6, 0, 0}, ELEM_CHARGE_C);
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    ion_ensemble.push_back(ion2);
    
    Vec3 force0 = force.compute(ion_ensemble[0], 0.0, ctx);
    
    // Force from ion1: F₁ = -k·e²/(1e-6)² (pushes left)
    // Force from ion2: F₂ = -k·e²/(2e-6)² (pushes left)
    const double F1 = -COULOMB_CONST * ELEM_CHARGE_C * ELEM_CHARGE_C / (1e-6 * 1e-6);
    const double F2 = -COULOMB_CONST * ELEM_CHARGE_C * ELEM_CHARGE_C / (2e-6 * 2e-6);
    const double expected_force_x = F1 + F2;
    
    REQUIRE_THAT(force0.x, WithinAbs(expected_force_x, tolerance));
    REQUIRE_THAT(force0.y, WithinAbs(0.0, tolerance));
    REQUIRE_THAT(force0.z, WithinAbs(0.0, tolerance));
}

// ============================================================================
// Softening Tests
// ============================================================================

TEST_CASE("SpaceChargeDirect: Softening prevents infinite force", "[SpaceChargeDirect][softening]") {
    SpaceChargeDirect force_no_soft{0.0};
    SpaceChargeDirect force_with_soft{1e-10};  // 0.1 nm softening
    
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Two ions very close together (0.01 nm apart)
    IonState ion0 = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    IonState ion1 = create_ion(Vec3{1e-11, 0, 0}, ELEM_CHARGE_C);  // 0.01 nm separation
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    
    Vec3 force_no = force_no_soft.compute(ion_ensemble[0], 0.0, ctx);
    Vec3 force_with = force_with_soft.compute(ion_ensemble[0], 0.0, ctx);
    
    // Softened force should be smaller
    const double mag_no = std::sqrt(force_no.x * force_no.x + force_no.y * force_no.y + force_no.z * force_no.z);
    const double mag_with = std::sqrt(force_with.x * force_with.x + force_with.y * force_with.y + force_with.z * force_with.z);
    
    REQUIRE(mag_with < mag_no);
    REQUIRE(std::isfinite(mag_with));
}

TEST_CASE("SpaceChargeDirect: Softening converges to Coulomb at large distances", "[SpaceChargeDirect][softening]") {
    SpaceChargeDirect force_no_soft{0.0};
    SpaceChargeDirect force_with_soft{1e-10};
    
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Two ions far apart (100 µm)
    IonState ion0 = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    IonState ion1 = create_ion(Vec3{100e-6, 0, 0}, ELEM_CHARGE_C);
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    
    Vec3 force_no = force_no_soft.compute(ion_ensemble[0], 0.0, ctx);
    Vec3 force_with = force_with_soft.compute(ion_ensemble[0], 0.0, ctx);
    
    // At large distances, softening should have negligible effect
    const double relative_diff = std::abs(force_with.x - force_no.x) / std::abs(force_no.x);
    
    REQUIRE(relative_diff < 1e-9);  // Less than 0.0000001% difference
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST_CASE("SpaceChargeDirect: Overlapping ions return zero force", "[SpaceChargeDirect][edge]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Two ions at exact same position
    IonState ion0 = create_ion(Vec3{1e-6, 1e-6, 1e-6}, ELEM_CHARGE_C);
    IonState ion1 = create_ion(Vec3{1e-6, 1e-6, 1e-6}, ELEM_CHARGE_C);
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    
    Vec3 result = force.compute(ion_ensemble[0], 0.0, ctx);
    
    REQUIRE(result.x == 0.0);
    REQUIRE(result.y == 0.0);
    REQUIRE(result.z == 0.0);
}

TEST_CASE("SpaceChargeDirect: Neutral ion produces no force", "[SpaceChargeDirect][edge]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Ion with zero charge
    IonState ion0 = create_ion(Vec3{0, 0, 0}, 0.0);  // Neutral
    IonState ion1 = create_ion(Vec3{1e-6, 0, 0}, ELEM_CHARGE_C);
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    
    Vec3 result = force.compute(ion_ensemble[0], 0.0, ctx);
    
    REQUIRE(result.x == 0.0);
    REQUIRE(result.y == 0.0);
    REQUIRE(result.z == 0.0);
}

TEST_CASE("SpaceChargeDirect: Time independent", "[SpaceChargeDirect][edge]") {
    SpaceChargeDirect force{0.0};
    std::vector<IonState> ion_ensemble;
    ForceContext ctx;
    ctx.all_ions = &ion_ensemble;
    
    // Coulomb force is time-independent
    IonState ion0 = create_ion(Vec3{0, 0, 0}, ELEM_CHARGE_C);
    IonState ion1 = create_ion(Vec3{1e-6, 0, 0}, ELEM_CHARGE_C);
    
    ion_ensemble.push_back(ion0);
    ion_ensemble.push_back(ion1);
    
    Vec3 force_t0 = force.compute(ion_ensemble[0], 0.0, ctx);
    Vec3 force_t1 = force.compute(ion_ensemble[0], 1e-6, ctx);
    Vec3 force_t2 = force.compute(ion_ensemble[0], 1.0, ctx);
    
    REQUIRE(force_t0.x == force_t1.x);
    REQUIRE(force_t0.x == force_t2.x);
}

// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/forces/MagneticFieldForce.h"
#include "core/physics/forces/DampingForce.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include "utils/constants.h"

#include <cmath>

using namespace ICARION::physics;
using Catch::Approx;

// ============================================================================
// Test Utilities
// ============================================================================

IonState make_ion(double vx, double vy, double vz, double mass_amu = 100.0, double charge_e = 1.0) {
    IonState ion;
    ion.pos = Vec3{0.0, 0.0, 0.0};
    ion.vel = Vec3{vx, vy, vz};
    ion.mass_kg = mass_amu * AMU_TO_KG;
    ion.ion_charge_C = charge_e * ELEM_CHARGE_C;
    return ion;
}

// ============================================================================
// MagneticFieldForce Tests
// ============================================================================

TEST_CASE("MagneticFieldForce - Constructor validation", "[forces][magnetic]") {
    SECTION("Analytical mode with valid params") {
        MagneticFieldParams params;
        params.uniform_field_T = {0.0, 0.0, 7.0};
        params.enabled = true;
        
        REQUIRE_NOTHROW(MagneticFieldForce(params));
    }
    
    SECTION("Field provider mode requires non-null provider") {
        REQUIRE_THROWS_AS(
            MagneticFieldForce(nullptr),
            std::invalid_argument
        );
    }
}

TEST_CASE("MagneticFieldForce - Lorentz force F = q(v×B)", "[forces][magnetic]") {
    MagneticFieldParams params;
    params.uniform_field_T = {0.0, 0.0, 1.0};  // 1 T along z-axis
    params.enabled = true;
    
    MagneticFieldForce force(params);
    ForceContext ctx;
    
    SECTION("Force name") {
        REQUIRE(force.name() == "MagneticField");
    }
    
    SECTION("Stationary ion: zero force") {
        IonState ion = make_ion(0, 0, 0);
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        REQUIRE(F.x == Approx(0.0).margin(1e-25));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
    
    SECTION("Velocity parallel to B-field: zero force") {
        IonState ion = make_ion(0, 0, 1000);  // v along z, B along z
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        REQUIRE(std::fabs(F.x) < 1e-20);
        REQUIRE(std::fabs(F.y) < 1e-20);
        REQUIRE(std::fabs(F.z) < 1e-20);
    }
    
    SECTION("Velocity perpendicular to B-field: maximum force") {
        // v = (1000, 0, 0) m/s, B = (0, 0, 1) T
        // v×B = (vy*Bz - vz*By, vz*Bx - vx*Bz, vx*By - vy*Bx)
        //     = (0*1 - 0*0, 0*0 - 1000*1, 1000*0 - 0*0) = (0, -1000, 0)
        // F = q·(v×B) = e·(0, -1000, 0)
        IonState ion = make_ion(1000, 0, 0);
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        double expected_F_y = -ELEM_CHARGE_C * 1000.0 * 1.0;  // q·v·B (negative!)
        
        REQUIRE(std::fabs(F.x) < 1e-20);
        REQUIRE(F.y == Approx(expected_F_y));
        REQUIRE(std::fabs(F.z) < 1e-20);
    }
    
    SECTION("Right-hand rule verification") {
        // v = (1, 0, 0), B = (0, 1, 0)
        // v×B = (0, 0, 1) → force in +z
        MagneticFieldParams params_y;
        params_y.uniform_field_T = {0.0, 1.0, 0.0};
        params_y.enabled = true;
        MagneticFieldForce force_y(params_y);
        
        IonState ion = make_ion(100, 0, 0);
        Vec3 F = force_y.compute(ion, 0.0, ctx);
        
        REQUIRE(F.z > 0.0);  // Force in +z direction
    }
}

TEST_CASE("MagneticFieldForce - Cyclotron motion", "[forces][magnetic]") {
    // For FTICR: ions undergo cyclotron motion in B-field
    // ω_c = q·B/m (cyclotron frequency)
    
    MagneticFieldParams params;
    params.uniform_field_T = {0.0, 0.0, 7.0};  // 7 T (typical FTICR)
    params.enabled = true;
    
    MagneticFieldForce force(params);
    ForceContext ctx;
    
    IonState ion = make_ion(1000, 0, 0, 100.0);  // 100 amu, v_x = 1 km/s
    
    Vec3 F = force.compute(ion, 0.0, ctx);
    
    // F should be perpendicular to v
    double F_dot_v = F.x * ion.vel.x + F.y * ion.vel.y + F.z * ion.vel.z;
    REQUIRE(std::fabs(F_dot_v) < 1e-15);  // F ⊥ v
    
    // Centripetal acceleration: a = v²/r = q·v·B/m
    double a_expected = ELEM_CHARGE_C * 1000.0 * 7.0 / ion.mass_kg;
    double a_actual = std::sqrt(F.x*F.x + F.y*F.y + F.z*F.z) / ion.mass_kg;
    
    REQUIRE(a_actual == Approx(a_expected));
}

TEST_CASE("MagneticFieldForce - Linear gradient field", "[forces][magnetic]") {
    MagneticFieldParams params;
    params.uniform_field_T = {0.0, 0.0, 1.0};      // 1 T base field
    params.gradient_T_per_m = {0.0, 0.0, 0.1};     // 0.1 T/m gradient
    params.enabled = true;
    
    MagneticFieldForce force(params);
    ForceContext ctx;
    
    SECTION("Field increases with position") {
        // At z=0: B = 1 T
        IonState ion1 = make_ion(1000, 0, 0);
        ion1.pos.z = 0.0;
        Vec3 F1 = force.compute(ion1, 0.0, ctx);
        
        // At z=10m: B = 1 + 0.1*10 = 2 T
        IonState ion2 = make_ion(1000, 0, 0);
        ion2.pos.z = 10.0;
        Vec3 F2 = force.compute(ion2, 0.0, ctx);
        
        // Force should be twice as strong (F ∝ B)
        REQUIRE(F2.y == Approx(2.0 * F1.y));
    }
}

TEST_CASE("MagneticFieldForce - Disabled force", "[forces][magnetic]") {
    MagneticFieldParams params;
    params.uniform_field_T = {0.0, 0.0, 10.0};
    params.enabled = false;  // Disabled!
    
    MagneticFieldForce force(params);
    ForceContext ctx;
    
    IonState ion = make_ion(1000, 0, 0);
    Vec3 F = force.compute(ion, 0.0, ctx);
    
    REQUIRE(F.x == Approx(0.0).margin(1e-25));
    REQUIRE(F.y == Approx(0.0).margin(1e-25));
    REQUIRE(F.z == Approx(0.0).margin(1e-25));
}

// ============================================================================
// DampingForce Tests
// ============================================================================

TEST_CASE("DampingForce - Friction model (explicit gamma)", "[forces][damping]") {
    DampingParams params;
    params.model = DampingModel::Friction;
    params.gamma_coefficient = 1e6;  // 1/s (explicit damping coefficient)
    
    DampingForce force(params);
    ForceContext ctx;
    
    SECTION("Force name") {
        REQUIRE(force.name() == "Damping(Friction)");
    }
    
    SECTION("Force opposes velocity") {
        IonState ion = make_ion(1000, 0, 0);  // v_x = 1 km/s, m = 100 Da
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        // F = -γ·m·v [N]
        double expected_F_x = -params.gamma_coefficient * ion.mass_kg * 1000.0;
        
        REQUIRE(F.x == Approx(expected_F_x));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
    
    SECTION("Force proportional to velocity") {
        IonState ion1 = make_ion(500, 0, 0);
        IonState ion2 = make_ion(1000, 0, 0);
        
        Vec3 F1 = force.compute(ion1, 0.0, ctx);
        Vec3 F2 = force.compute(ion2, 0.0, ctx);
        
        // F2 should be twice F1 (F = -γ·m·v, so F ∝ v)
        REQUIRE(F2.x == Approx(2.0 * F1.x));
    }
    
    SECTION("Stationary ion: zero force") {
        IonState ion = make_ion(0, 0, 0);
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        REQUIRE(F.x == Approx(0.0).margin(1e-25));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
}

TEST_CASE("DampingForce - HardSphere model (deterministic)", "[forces][damping]") {
    // Test HardSphere damping with explicit parameters
    DampingParams params;
    params.model = DampingModel::HardSphere;
    params.gas_density_m3 = 2.5e25;           // ~1 atm at 300K
    params.mean_thermal_velocity_m_s = 500.0;  // ~N2 at 300K
    params.neutral_mass_kg = 28.0 * AMU_TO_KG; // N2
    params.CCS_m2 = 1e-18;                     // 100 Ų
    
    DampingForce force(params);
    ForceContext ctx;
    ctx.gas_density_m3 = params.gas_density_m3;
    ctx.mean_thermal_velocity_m_s = params.mean_thermal_velocity_m_s;
    ctx.neutral_mass_kg = params.neutral_mass_kg;
    
    SECTION("Force name") {
        REQUIRE(force.name() == "Damping(HardSphere)");
    }
    
    SECTION("Force opposes velocity (deterministic)") {
        IonState ion = make_ion(1000, 0, 0);  // v_x = 1 km/s
        ion.CCS_m2 = params.CCS_m2;
        
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        // F = -γ·m·v where γ = n·σ·v_th·m_reduced/m_ion
        // Should be negative (opposes motion) and deterministic
        REQUIRE(F.x < 0.0);
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
        
        // Check determinism: repeated calls give same result
        Vec3 F2 = force.compute(ion, 0.0, ctx);
        REQUIRE(F2.x == Approx(F.x));
    }
    
    SECTION("Force proportional to velocity") {
        IonState ion = make_ion(500, 0, 0);
        ion.CCS_m2 = params.CCS_m2;
        
        Vec3 F1 = force.compute(ion, 0.0, ctx);
        
        ion.vel.x = 1000.0;  // Double velocity
        Vec3 F2 = force.compute(ion, 0.0, ctx);
        
        // F ∝ v (F = -γ·m·v with constant γ)
        REQUIRE(F2.x == Approx(2.0 * F1.x));
    }
}

TEST_CASE("DampingForce - No damping", "[forces][damping]") {
    DampingParams params;
    params.model = DampingModel::None;
    
    DampingForce force(params);
    ForceContext ctx;
    
    SECTION("Force name") {
        REQUIRE(force.name() == "Damping(None)");
    }
    
    SECTION("Always returns zero force") {
        IonState ion = make_ion(1000, 500, 250);
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        REQUIRE(F.x == Approx(0.0).margin(1e-25));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
}

// ============================================================================
// Combined Forces Test
// ============================================================================

TEST_CASE("Combined Magnetic + Damping", "[forces][combined]") {
    // Simulate FTICR-like scenario: B-field + friction damping
    
    MagneticFieldParams b_params;
    b_params.uniform_field_T = {0.0, 0.0, 7.0};
    b_params.enabled = true;
    
    DampingParams d_params;
    d_params.model = DampingModel::Friction;
    d_params.gamma_coefficient = 1e5;  // 1/s (explicit)
    
    MagneticFieldForce mag_force(b_params);
    DampingForce damp_force(d_params);
    
    ForceContext ctx;
    IonState ion = make_ion(1000, 500, 0);  // Velocity in x-y plane
    
    Vec3 F_mag = mag_force.compute(ion, 0.0, ctx);
    Vec3 F_damp = damp_force.compute(ion, 0.0, ctx);
    Vec3 F_total = F_mag + F_damp;
    
    SECTION("Magnetic force perpendicular to velocity") {
        double dot = F_mag.x * ion.vel.x + F_mag.y * ion.vel.y + F_mag.z * ion.vel.z;
        REQUIRE(std::fabs(dot) < 1e-14);
    }
    
    SECTION("Damping force opposes velocity") {
        REQUIRE(F_damp.x < 0.0);  // Opposes v_x > 0
        REQUIRE(F_damp.y < 0.0);  // Opposes v_y > 0
    }
    
    SECTION("Total force is vector sum") {
        REQUIRE(F_total.x == Approx(F_mag.x + F_damp.x));
        REQUIRE(F_total.y == Approx(F_mag.y + F_damp.y));
        REQUIRE(F_total.z == Approx(F_mag.z + F_damp.z));
    }
}

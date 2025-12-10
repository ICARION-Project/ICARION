// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/forces/MagneticFieldForce.h"
#include "core/physics/forces/DampingForce.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/types/Vec3.h"
#include "core/config/types/FieldsConfig.h"
#include "core/config/types/EnvironmentConfig.h"
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

Vec3 compute_force(const MagneticFieldForce& force, const IonState& ion, ForceContext ctx, double t = 0.0) {
    ICARION::core::IonEnsemble ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    ctx.ion_index = 0;
    return force.compute(ens, 0, t, ctx);
}

// ============================================================================
// MagneticFieldForce Tests
// ============================================================================

TEST_CASE("MagneticFieldForce - Constructor validation", "[forces][magnetic]") {
    SECTION("Analytical mode with valid config") {
        ICARION::config::MagneticFieldConfig config;
        config.field_strength_T = {0.0, 0.0, 7.0};
        config.enabled = true;
        
        REQUIRE_NOTHROW(MagneticFieldForce(config));
    }
    
    SECTION("Field provider mode requires non-null provider") {
        REQUIRE_THROWS_AS(
            MagneticFieldForce(nullptr),
            std::invalid_argument
        );
    }
}

TEST_CASE("MagneticFieldForce - Lorentz force F = q(v×B)", "[forces][magnetic]") {
    ICARION::config::MagneticFieldConfig config;
    config.field_strength_T = {0.0, 0.0, 1.0};  // 1 T along z-axis
    config.enabled = true;
    
    MagneticFieldForce force(config);
    ForceContext ctx;
    
    SECTION("Force name") {
        REQUIRE(force.name() == "MagneticField");
    }
    
    SECTION("Stationary ion: zero force") {
        IonState ion = make_ion(0, 0, 0);
        Vec3 F = compute_force(force, ion, ctx);
        
        REQUIRE(F.x == Approx(0.0).margin(1e-25));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
    
    SECTION("Velocity parallel to B-field: zero force") {
        IonState ion = make_ion(0, 0, 1000);  // v along z, B along z
        Vec3 F = compute_force(force, ion, ctx);
        
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
        Vec3 F = compute_force(force, ion, ctx);
        
        double expected_F_y = -ELEM_CHARGE_C * 1000.0 * 1.0;  // q·v·B (negative!)
        
        REQUIRE(std::fabs(F.x) < 1e-20);
        REQUIRE(F.y == Approx(expected_F_y));
        REQUIRE(std::fabs(F.z) < 1e-20);
    }
    
    SECTION("Right-hand rule verification") {
        // v = (1, 0, 0), B = (0, 1, 0)
        // v×B = (0, 0, 1) → force in +z
        ICARION::config::MagneticFieldConfig config_y;
        config_y.field_strength_T = {0.0, 1.0, 0.0};
        config_y.enabled = true;
        MagneticFieldForce force_y(config_y);
        
        IonState ion = make_ion(100, 0, 0);
        Vec3 F = compute_force(force_y, ion, ctx);
        
        REQUIRE(F.z > 0.0);  // Force in +z direction
    }
}

TEST_CASE("MagneticFieldForce - Cyclotron motion", "[forces][magnetic]") {
    // For FTICR: ions undergo cyclotron motion in B-field
    // ω_c = q·B/m (cyclotron frequency)
    
    ICARION::config::MagneticFieldConfig config;
    config.field_strength_T = {0.0, 0.0, 7.0};  // 7 T (typical FTICR)
    config.enabled = true;
    
    MagneticFieldForce force(config);
    ForceContext ctx;
    
    IonState ion = make_ion(1000, 0, 0, 100.0);  // 100 amu, v_x = 1 km/s
    
    Vec3 F = compute_force(force, ion, ctx);
    
    // F should be perpendicular to v
    double F_dot_v = F.x * ion.vel.x + F.y * ion.vel.y + F.z * ion.vel.z;
    REQUIRE(std::fabs(F_dot_v) < 1e-15);  // F ⊥ v
    
    // Centripetal acceleration: a = v²/r = q·v·B/m
    double a_expected = ELEM_CHARGE_C * 1000.0 * 7.0 / ion.mass_kg;
    double a_actual = std::sqrt(F.x*F.x + F.y*F.y + F.z*F.z) / ion.mass_kg;
    
    REQUIRE(a_actual == Approx(a_expected));
}

TEST_CASE("MagneticFieldForce - Linear gradient field", "[forces][magnetic]") {
    ICARION::config::MagneticFieldConfig config;
    config.field_strength_T = {0.0, 0.0, 1.0};      // 1 T base field
    config.field_gradient_T_m = {0.0, 0.0, 0.1};     // 0.1 T/m gradient
    config.enabled = true;
    
    MagneticFieldForce force(config);
    ForceContext ctx;
    
    SECTION("Field increases with position") {
        // At z=0: B = 1 T
        IonState ion1 = make_ion(1000, 0, 0);
        ion1.pos.z = 0.0;
        Vec3 F1 = compute_force(force, ion1, ctx);
        
        // At z=10m: B = 1 + 0.1*10 = 2 T
        IonState ion2 = make_ion(1000, 0, 0);
        ion2.pos.z = 10.0;
        Vec3 F2 = compute_force(force, ion2, ctx);
        
        // Force should be twice as strong (F ∝ B)
        REQUIRE(F2.y == Approx(2.0 * F1.y));
    }
}

TEST_CASE("MagneticFieldForce - Disabled force", "[forces][magnetic]") {
    ICARION::config::MagneticFieldConfig config;
    config.field_strength_T = {0.0, 0.0, 10.0};
    config.enabled = false;  // Disabled!
    
    MagneticFieldForce force(config);
    ForceContext ctx;
    
    IonState ion = make_ion(1000, 0, 0);
    Vec3 F = compute_force(force, ion, ctx);
    
    REQUIRE(F.x == Approx(0.0).margin(1e-25));
    REQUIRE(F.y == Approx(0.0).margin(1e-25));
    REQUIRE(F.z == Approx(0.0).margin(1e-25));
}

// ============================================================================
// DampingForce Tests
// ============================================================================

TEST_CASE("DampingForce - Friction model (mobility-based)", "[forces][damping]") {
    ICARION::config::EnvironmentConfig env;
    env.pressure_Pa = 101325.0;
    env.temperature_K = 300.0;
    env.compute_derived_properties();
    
    DampingForce force(env, DampingModel::Friction);
    ForceContext ctx;
    
    SECTION("Force name") {
        REQUIRE(force.name() == "Damping(Friction)");
    }
    
    SECTION("Force opposes velocity") {
        IonState ion = make_ion(1000, 0, 0);  // v_x = 1 km/s, m = 100 Da
        ion.reduced_mobility_cm2_Vs = 2.0;  // Need mobility for Friction model
        Vec3 F = compute_force(force, ion, ctx);
        
        // F = -γ·m·v [N] where γ = q/(K·m)
        // K = K₀ · (n/n₀) where n = gas_density, n₀ = LOSCHMIDT_CONSTANT
        double K0_m2_Vs = ion.reduced_mobility_cm2_Vs * 1e-4;
        double K_m2_Vs = K0_m2_Vs * LOSCHMIDT_CONSTANT / env.particle_density_m_3;
        double gamma = ion.ion_charge_C / (K_m2_Vs * ion.mass_kg);
        double expected_F_x = -gamma * ion.mass_kg * 1000.0;
        
        REQUIRE(F.x == Approx(expected_F_x));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
    
    SECTION("Force proportional to velocity") {
        IonState ion1 = make_ion(500, 0, 0);
        IonState ion2 = make_ion(1000, 0, 0);
        
        Vec3 F1 = compute_force(force, ion1, ctx);
        Vec3 F2 = compute_force(force, ion2, ctx);
        
        // F2 should be twice F1 (F = -γ·m·v, so F ∝ v)
        REQUIRE(F2.x == Approx(2.0 * F1.x));
    }
    
    SECTION("Stationary ion: zero force") {
        IonState ion = make_ion(0, 0, 0);
        Vec3 F = compute_force(force, ion, ctx);
        
        REQUIRE(F.x == Approx(0.0).margin(1e-25));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
}

TEST_CASE("DampingForce - HardSphere model (deterministic)", "[forces][damping]") {
    // Test HardSphere damping with EnvironmentConfig
    ICARION::config::EnvironmentConfig env;
    env.pressure_Pa = 101325.0;  // 1 atm
    env.temperature_K = 300.0;
    env.gas_species = "N2";
    env.compute_derived_properties();
    
    DampingForce force(env, DampingModel::HardSphere);
    ForceContext ctx;
    
    SECTION("Force name") {
        REQUIRE(force.name() == "Damping(HardSphere)");
    }
    
    SECTION("Force opposes velocity (deterministic)") {
        IonState ion = make_ion(1000, 0, 0);  // v_x = 1 km/s
        ion.CCS_m2 = 1e-18;  // 100 Ų

        Vec3 F = compute_force(force, ion, ctx);
        
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
        ion.CCS_m2 = 1e-18;  // 100 Ų
        
        Vec3 F1 = compute_force(force, ion, ctx);
        
        ion.vel.x = 1000.0;  // Double velocity
        Vec3 F2 = compute_force(force, ion, ctx);
        
        // F ∝ v (F = -γ·m·v with constant γ)
        REQUIRE(F2.x == Approx(2.0 * F1.x));
    }
}

TEST_CASE("DampingForce - No damping", "[forces][damping]") {
    ICARION::config::EnvironmentConfig env;
    env.compute_derived_properties();
    
    DampingForce force(env, DampingModel::None);
    ForceContext ctx;
    
    SECTION("Force name") {
        REQUIRE(force.name() == "Damping(None)");
    }
    
    SECTION("Always returns zero force") {
        IonState ion = make_ion(1000, 500, 250);
        Vec3 F = compute_force(force, ion, ctx);
        
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
    
    ICARION::config::MagneticFieldConfig mag_config;
    mag_config.field_strength_T = {0.0, 0.0, 7.0};
    mag_config.enabled = true;
    
    ICARION::config::EnvironmentConfig env;
    env.pressure_Pa = 101325.0;
    env.temperature_K = 300.0;
    env.compute_derived_properties();
    
    MagneticFieldForce mag_force(mag_config);
    DampingForce damp_force(env, DampingModel::Friction);
    
    ForceContext ctx;
    IonState ion = make_ion(1000, 500, 0);  // Velocity in x-y plane
    ion.reduced_mobility_cm2_Vs = 2.0;  // Need mobility
    
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

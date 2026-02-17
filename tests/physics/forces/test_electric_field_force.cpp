// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/forces/ElectricFieldForce.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include "fieldsolver/utils/IFieldProvider.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/FieldsConfig.h"
#include "core/config/types/GeometryConfig.h"
#include "utils/constants.h"

#include <cmath>

using namespace ICARION::physics;
using namespace ICARION::instrument;
using Catch::Approx;

// ============================================================================
// Mock Field Provider for Testing
// ============================================================================

class MockFieldProvider : public IFieldProvider {
public:
    MockFieldProvider(const Vec3& uniform_field) : uniform_field_(uniform_field) {}
    
    Vec3 get_E(const Vec3& pos) const override {
        (void)pos;  // Position-independent uniform field
        return uniform_field_;
    }
    
private:
    Vec3 uniform_field_;
};

// ============================================================================
// Test Utilities
// ============================================================================

// Create a test ion at position (x, y, z) with charge q
IonState make_test_ion(double x, double y, double z, double charge_C = ELEM_CHARGE_C) {
    IonState ion;
    ion.pos = Vec3{x, y, z};
    ion.mass_kg = 100.0 * AMU_TO_KG;  // 100 amu
    ion.ion_charge_C = charge_C;
    return ion;
}

Vec3 compute_force(ElectricFieldForce& force, const IonState& ion, ForceContext ctx, double t = 0.0) {
    ICARION::core::IonEnsemble ens = ICARION::core::IonEnsemble::from_legacy({ion});
    ctx.ion_ensemble = &ens;
    ctx.ion_index = 0;
    return force.compute(ens, 0, t, ctx);
}

// ============================================================================
// Basic Functionality Tests
// ============================================================================

TEST_CASE("ElectricFieldForce - Constructor validation", "[forces][electric]") {
    SECTION("Analytical mode with valid config") {
        ICARION::config::DomainConfig domain;
        domain.instrument = ICARION::config::Instrument::IMS;
        domain.geometry.length_m = 0.1;
        domain.fields.dc.axial_V.constant_value = 1000.0;
        
        REQUIRE_NOTHROW(ElectricFieldForce(domain));
    }
    
    SECTION("Field provider mode requires non-null provider") {
        REQUIRE_THROWS_AS(
            ElectricFieldForce(nullptr),
            std::invalid_argument
        );
    }
    
    SECTION("Valid field provider constructor") {
        auto provider = std::make_shared<MockFieldProvider>(Vec3{100, 0, 0});
        REQUIRE_NOTHROW(ElectricFieldForce(provider));
    }
}

TEST_CASE("ElectricFieldForce - Field provider mode", "[forces][electric]") {
    // Uniform field: E = (1000, 0, 0) V/m
    auto provider = std::make_shared<MockFieldProvider>(Vec3{1000.0, 0.0, 0.0});
    ElectricFieldForce force(provider);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(FieldProvider)");
    }
    
    SECTION("Force proportional to charge") {
        IonState ion1 = make_test_ion(0, 0, 0, ELEM_CHARGE_C);
        IonState ion2 = make_test_ion(0, 0, 0, 2.0 * ELEM_CHARGE_C);
        
        ForceContext ctx;
        Vec3 F1 = compute_force(force, ion1, ctx);
        Vec3 F2 = compute_force(force, ion2, ctx);
        
        // F = q*E, so F2 should be twice F1
        REQUIRE(F2.x == Approx(2.0 * F1.x));
    }
    
    SECTION("Force calculation: F = q*E") {
        IonState ion = make_test_ion(0, 0, 0, ELEM_CHARGE_C);
        
        ForceContext ctx;
        Vec3 F = compute_force(force, ion, ctx);
        
        // F = q * E = 1.602e-19 * 1000 = 1.602e-16 N
        REQUIRE(F.x == Approx(ELEM_CHARGE_C * 1000.0));
        REQUIRE(F.y == Approx(0.0));
        REQUIRE(F.z == Approx(0.0));
    }
    
    SECTION("Context field provider overrides constructor provider") {
        // Constructor provider: E = (1000, 0, 0)
        // Context provider: E = (0, 2000, 0)
        auto ctx_provider = std::make_shared<MockFieldProvider>(Vec3{0.0, 2000.0, 0.0});
        ForceContext ctx;
        ctx.field_provider = ctx_provider.get();
        
        IonState ion = make_test_ion(0, 0, 0, ELEM_CHARGE_C);
        Vec3 F = compute_force(force, ion, ctx);
        
        // Should use context provider (y-direction field)
        REQUIRE(F.x == Approx(0.0).margin(1e-20));
        REQUIRE(F.y == Approx(ELEM_CHARGE_C * 2000.0));
        REQUIRE(F.z == Approx(0.0).margin(1e-20));
    }
}

// ============================================================================
// IMS: Ion Mobility Spectrometry
// ============================================================================

TEST_CASE("ElectricFieldForce - IMS (Ion Mobility Spectrometry)", "[forces][electric][ims]") {
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::IMS;
    domain.geometry.length_m = 0.1;  // 10 cm drift tube
    // Initialize all field values
    domain.fields.dc.axial_V.constant_value = 1000.0;  // 1000 V -> 10 kV/m
    domain.fields.dc.quad_V.constant_value = 0.0;  // Initialize
    domain.fields.rf.voltage_V.constant_value = 0.0;  // Initialize
    domain.fields.rf.frequency_Hz.constant_value = 0.0;  // Initialize
    domain.fields.rf.compute_derived();
    
    ElectricFieldForce force(domain);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(IMS)");
    }
    
    SECTION("Uniform axial field") {
        // E_z = V/L = 1000/0.1 = 10000 V/m
        IonState ion1 = make_test_ion(0.00, 0.00, 0.00);
        IonState ion2 = make_test_ion(0.01, 0.02, 0.05);
        
        ForceContext ctx;
        Vec3 F1 = compute_force(force, ion1, ctx);
        Vec3 F2 = compute_force(force, ion2, ctx);
        
        // Field should be position-independent
        REQUIRE(F1.x == Approx(F2.x).margin(1e-20));
        REQUIRE(F1.y == Approx(F2.y).margin(1e-20));
        REQUIRE(F1.z == Approx(F2.z));
    }
    
    SECTION("Correct field magnitude") {
        IonState ion = make_test_ion(0, 0, 0);
        
        ForceContext ctx;
        Vec3 F = compute_force(force, ion, ctx);
        
        // F_z = q * E_z = q * (V/L)
        double expected_F_z = ELEM_CHARGE_C * (1000.0 / 0.1);
        REQUIRE(F.z == Approx(expected_F_z));
        REQUIRE(F.x == Approx(0.0).margin(1e-20));
        REQUIRE(F.y == Approx(0.0).margin(1e-20));
    }
}

// ============================================================================
// TOF: Time-of-Flight
// ============================================================================

TEST_CASE("ElectricFieldForce - TOF (Time-of-Flight)", "[forces][electric][tof]") {
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::TOF;
    domain.geometry.length_m = 1.0;  // 1 m flight tube
    domain.geometry.acc_length_m = 0.05;  // 5 cm acceleration region
    domain.fields.dc.axial_V.constant_value = 20000.0;  // 20 kV acceleration
    
    ElectricFieldForce force(domain);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(TOF)");
    }
    
    SECTION("Uniform extraction field in acceleration region") {
        IonState ion = make_test_ion(0, 0, 0.025);  // z = 2.5 cm (within acc region)
        
        ForceContext ctx;
        Vec3 F = compute_force(force, ion, ctx);
        
        // E_z = V_axial / acc_length = 20000 / 0.05 = 400000 V/m
        double expected_F_z = ELEM_CHARGE_C * 400000.0;
        REQUIRE(F.z == Approx(expected_F_z));
    }
    
    SECTION("Field-free drift region") {
        IonState ion = make_test_ion(0, 0, 0.5);  // z = 0.5 m (beyond acc region)
        
        ForceContext ctx;
        Vec3 F = compute_force(force, ion, ctx);
        
        // Should be zero in drift region
        REQUIRE(F.x == Approx(0.0).margin(1e-20));
        REQUIRE(F.y == Approx(0.0).margin(1e-20));
        REQUIRE(F.z == Approx(0.0).margin(1e-20));
    }
}

// ============================================================================
// LQIT: Linear Quadrupole Ion Trap
// ============================================================================

TEST_CASE("ElectricFieldForce - LQIT (Linear Quadrupole Ion Trap)", "[forces][electric][lqit]") {
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::LQIT;
    domain.geometry.radius_m = 0.005;  // 5 mm trap radius
    domain.geometry.length_m = 0.05;  // 5 cm length
    // Initialize all field values (required by ValueOrWaveform)
    domain.fields.rf.voltage_V.constant_value = 1000.0;  // 1 kV RF
    domain.fields.rf.frequency_Hz.constant_value = 1e6;  // 1 MHz
    domain.fields.rf.compute_derived();
    domain.fields.dc.quad_V.constant_value = 0.0;  // No DC offset
    domain.fields.dc.axial_V.constant_value = 0.0;  // Initialize
    domain.fields.ac.voltage_V.constant_value = 0.0;  // Initialize
    domain.fields.ac.frequency_Hz.constant_value = 0.0;  // Initialize
    domain.fields.ac.compute_derived();
    
    ElectricFieldForce force(domain);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(LQIT)");
    }
    
    SECTION("RF field at t=0 (cos(0)=1)") {
        IonState ion = make_test_ion(0.001, 0.002, 0.0);  // x=1mm, y=2mm
        
        ForceContext ctx;
        Vec3 F = compute_force(force, ion, ctx);
        
        // At t=0: U_eff = U_RF * cos(0) = 1000 V
        // E_x = 2*x*U_eff/r0^2, E_y = -2*y*U_eff/r0^2
        double r0_sq = domain.geometry.radius_m * domain.geometry.radius_m;
        double E_x = -2.0 * ion.pos.x * 1000.0 / r0_sq;
        double E_y = 2.0 * ion.pos.y * 1000.0 / r0_sq;
        
        REQUIRE(F.x == Approx(ELEM_CHARGE_C * E_x));
        REQUIRE(F.y == Approx(ELEM_CHARGE_C * E_y));
        REQUIRE(F.z == Approx(0.0).margin(1e-20));
    }
    
    SECTION("RF field time-dependence") {
        IonState ion = make_test_ion(0.001, 0.0, 0.0);
        ForceContext ctx;
        
        // At t=0: cos(ωt) = 1
        Vec3 F0 = compute_force(force, ion, ctx);
        
        // At t = T/4: cos(ωt) = 0 (quarter period)
        double T = 1.0 / domain.fields.rf.frequency_Hz.constant_value.value();
        Vec3 F_quarter = compute_force(force, ion, ctx, T / 4.0);
        
        // At t = T/2: cos(ωt) = -1 (half period)
        Vec3 F_half = compute_force(force, ion, ctx, T / 2.0);
        
        // Force should oscillate
        REQUIRE(std::fabs(F0.x) > 1e-20);  // Non-zero at t=0
        REQUIRE(std::fabs(F_quarter.x) < 1e-18);  // Zero at quarter period
        REQUIRE(F_half.x == Approx(-F0.x));  // Opposite at half period
    }
    
    SECTION("LQIT with DC offset") {
        ICARION::config::DomainConfig domain_dc = domain;
        domain_dc.fields.dc.quad_V.constant_value = 100.0;  // 100 V DC offset
        
        ElectricFieldForce force_dc(domain_dc);
        
        IonState ion = make_test_ion(0.001, 0.0, 0.0);
        ForceContext ctx;
        
        // At t=0: U_eff = DC + RF*cos(0) = 100 + 1000 = 1100 V
        Vec3 F = compute_force(force_dc, ion, ctx);
        
        double r0_sq = domain_dc.geometry.radius_m * domain_dc.geometry.radius_m;
        double E_x = -2.0 * ion.pos.x * 1100.0 / r0_sq;
        
        REQUIRE(F.x == Approx(ELEM_CHARGE_C * E_x));
    }
    
    SECTION("LQIT with DC endcap field") {
        ICARION::config::DomainConfig domain_axial = domain;
        domain_axial.fields.dc.axial_V.constant_value = 500.0;  // 500 V endcap potential
        domain_axial.geometry.length_m = 0.1;  // 10 cm trap
        
        ElectricFieldForce force_axial(domain_axial);
        ForceContext ctx;
        
        // Test center region: should be field-free
        IonState ion_center = make_test_ion(0.0, 0.0, 0.05);  // z = L/2 (center)
        Vec3 F_center = compute_force(force_axial, ion_center, ctx);
        REQUIRE(std::fabs(F_center.z) < 1e-20);  // No field in center
        
        // Test left endcap (z < 0.1*L): should push right (E_z > 0)
        IonState ion_left = make_test_ion(0.0, 0.0, 0.005);  // z = 5 mm (in left 10%)
        Vec3 F_left = compute_force(force_axial, ion_left, ctx);
        REQUIRE(F_left.z > 0.0);  // Pushes away from left endcap
        
        // Test right endcap (z > 0.9*L): should push left (E_z < 0)
        IonState ion_right = make_test_ion(0.0, 0.0, 0.095);  // z = 95 mm (in right 10%)
        Vec3 F_right = compute_force(force_axial, ion_right, ctx);
        REQUIRE(F_right.z < 0.0);  // Pushes away from right endcap
    }
    
    SECTION("LQIT with AC field (fixed x-direction for v1.0.0)") {
        ICARION::config::DomainConfig domain_ac;
        domain_ac.instrument = ICARION::config::Instrument::LQIT;
        domain_ac.geometry.radius_m = 0.005;
        domain_ac.geometry.length_m = 0.1;
        // Initialize all field values
        domain_ac.fields.rf.voltage_V.constant_value = 0.0;  // No RF for clean AC test
        domain_ac.fields.rf.frequency_Hz.constant_value = 0.0;  // Initialize
        domain_ac.fields.rf.compute_derived();
        domain_ac.fields.dc.quad_V.constant_value = 0.0;  // Initialize
        domain_ac.fields.dc.axial_V.constant_value = 0.0;  // Initialize
        domain_ac.fields.ac.voltage_V.constant_value = 100.0;  // 100 V AC
        domain_ac.fields.ac.frequency_Hz.constant_value = 5e5;  // 500 kHz
        domain_ac.fields.ac.compute_derived();
        
        ElectricFieldForce force_ac(domain_ac);
        
        IonState ion = make_test_ion(0.0, 0.0, 0.05);  // At trap center to isolate AC field
        ForceContext ctx;
        
        // At t=0: cos(0) = 1, AC field magnitude: E_x = -V_ac / r0
        Vec3 F0 = compute_force(force_ac, ion, ctx);
        double expected_E_ac = -domain_ac.fields.ac.voltage_V.constant_value.value() / domain_ac.geometry.radius_m;
        REQUIRE(F0.x == Approx(ELEM_CHARGE_C * expected_E_ac));
        
        // At t=T/4: cos(ωt) = 0, AC field should be zero
        double T = 1.0 / domain_ac.fields.ac.frequency_Hz.constant_value.value();
        Vec3 F_quarter = compute_force(force_ac, ion, ctx, T / 4.0);
        REQUIRE(std::fabs(F_quarter.x) <= 1e-18);
        
        // At t=T/2: cos(ωt) = -1, AC field should be reversed
        Vec3 F_half = compute_force(force_ac, ion, ctx, T / 2.0);
        REQUIRE(F_half.x == Approx(-F0.x));
        
        // Verify oscillation in x-direction only (no y or z component)
        REQUIRE(std::fabs(F0.y) < 1e-20);
        REQUIRE(std::fabs(F0.z) < 1e-20);
    }
}

// ============================================================================
// FTICR: Fourier Transform Ion Cyclotron Resonance
// ============================================================================

TEST_CASE("ElectricFieldForce - FTICR", "[forces][electric][fticr]") {
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::FTICR;
    domain.fields.dc.radial_V.constant_value = 10.0;  // 10 V trapping voltage
    domain.geometry.radius_m = 0.05;  // 5 cm radius
    domain.geometry.length_m = 0.1;  // 10 cm cell length
    
    ElectricFieldForce force(domain);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(FTICR)");
    }
    
    SECTION("Radial field increases with radius") {
        IonState ion1 = make_test_ion(0.001, 0.0, 0.05);  // r = 1 mm
        IonState ion2 = make_test_ion(0.002, 0.0, 0.05);  // r = 2 mm
        
        ForceContext ctx;
        Vec3 F1 = compute_force(force, ion1, ctx);
        Vec3 F2 = compute_force(force, ion2, ctx);
        
        // E_x = factor * x, so E should double when x doubles
        REQUIRE(F2.x == Approx(2.0 * F1.x));
    }
    
    SECTION("Axial field restoring (harmonic)") {
        // Field should be zero at center, increase linearly with z
        IonState ion_center = make_test_ion(0.0, 0.0, 0.05);  // At z_center
        IonState ion_offset = make_test_ion(0.0, 0.0, 0.06);  // 1 cm off-center
        
        ForceContext ctx;
        Vec3 F_center = compute_force(force, ion_center, ctx);
        Vec3 F_offset = compute_force(force, ion_offset, ctx);
        
        // At center: z - 0.5*L = 0.05 - 0.05 = 0, so F_z should be ~0
        REQUIRE(std::fabs(F_center.z) < 1e-18);
        
        // Off-center: should have restoring force
        REQUIRE(std::fabs(F_offset.z) > 1e-20);
    }
}

// ============================================================================
// Orbitrap
// ============================================================================

TEST_CASE("ElectricFieldForce - Orbitrap", "[forces][electric][orbitrap]") {
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::Orbitrap;
    domain.fields.dc.radial_V.constant_value = 3500.0;  // 3.5 kV
    domain.geometry.radius_char_m = 0.01;  // 1 cm characteristic radius
    domain.geometry.radius_in_m = 0.008;   // 8 mm inner
    domain.geometry.radius_out_m = 0.012;  // 12 mm outer
    domain.geometry.length_m = 0.05;  // 5 cm trap length
    domain.geometry.origin_m = {0.0, 0.0, 0.025};  // Trap center at z=0 (local)
    
    ElectricFieldForce force(domain);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(Orbitrap)");
    }
    
    SECTION("Axial harmonic confinement") {
        // At trap center (z = L/2), axial force should be zero
        IonState ion_center = make_test_ion(0.005, 0.0, 0.025);
        
        ForceContext ctx;
        Vec3 F_center = compute_force(force, ion_center, ctx);
        
        // z_center = 0.025 - 0.025 = 0
        REQUIRE(std::fabs(F_center.z) < 1e-18);
    }
    
    SECTION("Axial restoring force") {
        // Off-axis: should have linear restoring force
        IonState ion1 = make_test_ion(0.005, 0.0, 0.030);  // 0.5 cm above center
        IonState ion2 = make_test_ion(0.005, 0.0, 0.020);  // 0.5 cm below center
        
        ForceContext ctx;
        Vec3 F1 = compute_force(force, ion1, ctx);
        Vec3 F2 = compute_force(force, ion2, ctx);
        
        // Forces should be opposite (restoring)
        REQUIRE(F1.z == Approx(-F2.z));
        
        // Force proportional to displacement: F_z = -k * z_center
        // Compute k from config (same as ElectricFieldForce::compute_orbitrap_field)
        double r_char_sq = domain.geometry.radius_char_m * domain.geometry.radius_char_m;
        double k = 2.0 * domain.fields.dc.radial_V.constant_value.value() / 
                   (r_char_sq * std::log(domain.geometry.radius_out_m / domain.geometry.radius_in_m)
                    - 0.5 * (domain.geometry.radius_out_m * domain.geometry.radius_out_m 
                             - domain.geometry.radius_in_m * domain.geometry.radius_in_m));
        double z_center1 = 0.030 - 0.025;
        double expected_F_z1 = -ELEM_CHARGE_C * k * z_center1;
        REQUIRE(F1.z == Approx(expected_F_z1));
    }
    
    SECTION("Radial field depends on radius") {
        // r_char = 0.01 m, so test at r < r_char and r > r_char
        IonState ion1 = make_test_ion(0.005, 0.0, 0.025);  // r = 5 mm < r_char
        IonState ion2 = make_test_ion(0.015, 0.0, 0.025);  // r = 15 mm > r_char
        
        ForceContext ctx;
        Vec3 F1 = compute_force(force, ion1, ctx);
        Vec3 F2 = compute_force(force, ion2, ctx);
        
        // Radial field is not simply proportional (hyperlogarithmic)
        // Just check that both have radial components
        REQUIRE(std::fabs(F1.x) > 1e-20);
        REQUIRE(std::fabs(F2.x) > 1e-20);
    }
}

// ============================================================================
// QuadrupoleRF (generic, includes SLIM)
// ============================================================================

TEST_CASE("ElectricFieldForce - QuadrupoleRF", "[forces][electric][quadrupole]") {
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::QuadrupoleRF;
    domain.geometry.radius_m = 0.003;  // 3 mm
    // Initialize all field values
    domain.fields.rf.voltage_V.constant_value = 500.0;
    domain.fields.rf.frequency_Hz.constant_value = 5e5;  // 500 kHz
    domain.fields.rf.compute_derived();
    domain.fields.dc.quad_V.constant_value = 0.0;  // Initialize
    domain.fields.dc.axial_V.constant_value = 0.0;  // Initialize
    domain.fields.ac.voltage_V.constant_value = 0.0;  // Initialize
    domain.fields.ac.frequency_Hz.constant_value = 0.0;  // Initialize
    domain.fields.ac.compute_derived();
    
    ElectricFieldForce force(domain);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(QuadrupoleRF)");
    }
    
    SECTION("Behaves like LQIT") {
        // QuadrupoleRF uses same formula as LQIT
        IonState ion = make_test_ion(0.001, 0.001, 0.0);
        
        ForceContext ctx;
        Vec3 F = compute_force(force, ion, ctx);
        
        // Should have radial components
        REQUIRE(std::fabs(F.x) > 1e-20);
        REQUIRE(std::fabs(F.y) > 1e-20);
    }
}

// ============================================================================
// Edge Cases and Validation
// ============================================================================

TEST_CASE("ElectricFieldForce - Edge cases", "[forces][electric][edge]") {
    SECTION("Zero field produces zero force") {
        auto provider = std::make_shared<MockFieldProvider>(Vec3{0, 0, 0});
        ElectricFieldForce force(provider);
        
        IonState ion = make_test_ion(1, 2, 3);
        ForceContext ctx;
        Vec3 F = compute_force(force, ion, ctx);
        
        REQUIRE(F.x == Approx(0.0).margin(1e-25));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
    
    SECTION("Negative charge reverses force direction") {
        auto provider = std::make_shared<MockFieldProvider>(Vec3{1000, 0, 0});
        ElectricFieldForce force(provider);
        
        IonState ion_pos = make_test_ion(0, 0, 0, +ELEM_CHARGE_C);
        IonState ion_neg = make_test_ion(0, 0, 0, -ELEM_CHARGE_C);
        
        ForceContext ctx;
        Vec3 F_pos = compute_force(force, ion_pos, ctx);
        Vec3 F_neg = compute_force(force, ion_neg, ctx);
        
        REQUIRE(F_neg.x == Approx(-F_pos.x));
    }
    
    SECTION("NoFixedInstrument returns zero field") {
        ICARION::config::DomainConfig domain;
        domain.instrument = ICARION::config::Instrument::NoFixedInstrument;
        
        ElectricFieldForce force(domain);
        
        IonState ion = make_test_ion(1, 2, 3);
        ForceContext ctx;
        Vec3 F = compute_force(force, ion, ctx);
        
        REQUIRE(F.x == Approx(0.0).margin(1e-25));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
}

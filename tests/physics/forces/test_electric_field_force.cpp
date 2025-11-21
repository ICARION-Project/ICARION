// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/forces/ElectricFieldForce.h"
#include "core/physics/forces/ForceContext.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include "fieldsolver/utils/IFieldProvider.h"
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

// ============================================================================
// Basic Functionality Tests
// ============================================================================

TEST_CASE("ElectricFieldForce - Constructor validation", "[forces][electric]") {
    SECTION("Analytical mode requires valid instrument type") {
        AnalyticalFieldParams params;
        params.instrument_type = InstrumentType::UnknownInstrument;
        
        REQUIRE_THROWS_AS(
            ElectricFieldForce(params),
            std::invalid_argument
        );
    }
    
    SECTION("Field provider mode requires non-null provider") {
        REQUIRE_THROWS_AS(
            ElectricFieldForce(nullptr),
            std::invalid_argument
        );
    }
    
    SECTION("Valid analytical constructor") {
        AnalyticalFieldParams params;
        params.instrument_type = InstrumentType::IMS;
        params.length_m = 0.1;
        params.dc_axial_voltage_V = 1000.0;
        
        REQUIRE_NOTHROW(ElectricFieldForce(params));
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
        Vec3 F1 = force.compute(ion1, 0.0, ctx);
        Vec3 F2 = force.compute(ion2, 0.0, ctx);
        
        // F = q*E, so F2 should be twice F1
        REQUIRE(F2.x == Approx(2.0 * F1.x));
    }
    
    SECTION("Force calculation: F = q*E") {
        IonState ion = make_test_ion(0, 0, 0, ELEM_CHARGE_C);
        
        ForceContext ctx;
        Vec3 F = force.compute(ion, 0.0, ctx);
        
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
        Vec3 F = force.compute(ion, 0.0, ctx);
        
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
    AnalyticalFieldParams params;
    params.instrument_type = InstrumentType::IMS;
    params.length_m = 0.1;  // 10 cm drift tube
    params.dc_axial_voltage_V = 1000.0;  // 1000 V -> 10 kV/m
    
    ElectricFieldForce force(params);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(IMS)");
    }
    
    SECTION("Uniform axial field") {
        // E_z = V/L = 1000/0.1 = 10000 V/m
        IonState ion1 = make_test_ion(0.00, 0.00, 0.00);
        IonState ion2 = make_test_ion(0.01, 0.02, 0.05);
        
        ForceContext ctx;
        Vec3 F1 = force.compute(ion1, 0.0, ctx);
        Vec3 F2 = force.compute(ion2, 0.0, ctx);
        
        // Field should be position-independent
        REQUIRE(F1.x == Approx(F2.x).margin(1e-20));
        REQUIRE(F1.y == Approx(F2.y).margin(1e-20));
        REQUIRE(F1.z == Approx(F2.z));
    }
    
    SECTION("Correct field magnitude") {
        IonState ion = make_test_ion(0, 0, 0);
        
        ForceContext ctx;
        Vec3 F = force.compute(ion, 0.0, ctx);
        
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
    AnalyticalFieldParams params;
    params.instrument_type = InstrumentType::TOF;
    params.length_m = 1.0;  // 1 m flight tube
    params.dc_axial_voltage_V = 5000.0;  // 5 kV extraction
    
    ElectricFieldForce force(params);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(TOF)");
    }
    
    SECTION("Uniform extraction field") {
        IonState ion = make_test_ion(0, 0, 0.5);
        
        ForceContext ctx;
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        // E_z = 5000/1.0 = 5000 V/m
        double expected_F_z = ELEM_CHARGE_C * 5000.0;
        REQUIRE(F.z == Approx(expected_F_z));
    }
}

// ============================================================================
// LQIT: Linear Quadrupole Ion Trap
// ============================================================================

TEST_CASE("ElectricFieldForce - LQIT (Linear Quadrupole Ion Trap)", "[forces][electric][lqit]") {
    AnalyticalFieldParams params;
    params.instrument_type = InstrumentType::LQIT;
    params.radius_m = 0.005;  // 5 mm trap radius
    params.rf_voltage_V = 1000.0;  // 1 kV RF
    params.rf_frequency_Hz = 1e6;  // 1 MHz
    params.dc_quad_voltage_V = 0.0;  // No DC offset
    
    ElectricFieldForce force(params);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(LQIT)");
    }
    
    SECTION("RF field at t=0 (cos(0)=1)") {
        IonState ion = make_test_ion(0.001, 0.002, 0.0);  // x=1mm, y=2mm
        
        ForceContext ctx;
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        // At t=0: U_eff = U_RF * cos(0) = 1000 V
        // E_x = 2*x*U_eff/r0^2, E_y = -2*y*U_eff/r0^2
        double r0_sq = params.radius_m * params.radius_m;
        double E_x = 2.0 * ion.pos.x * 1000.0 / r0_sq;
        double E_y = -2.0 * ion.pos.y * 1000.0 / r0_sq;
        
        REQUIRE(F.x == Approx(ELEM_CHARGE_C * E_x));
        REQUIRE(F.y == Approx(ELEM_CHARGE_C * E_y));
        REQUIRE(F.z == Approx(0.0).margin(1e-20));
    }
    
    SECTION("RF field time-dependence") {
        IonState ion = make_test_ion(0.001, 0.0, 0.0);
        ForceContext ctx;
        
        // At t=0: cos(ωt) = 1
        Vec3 F0 = force.compute(ion, 0.0, ctx);
        
        // At t = T/4: cos(ωt) = 0 (quarter period)
        double T = 1.0 / params.rf_frequency_Hz;
        Vec3 F_quarter = force.compute(ion, T / 4.0, ctx);
        
        // At t = T/2: cos(ωt) = -1 (half period)
        Vec3 F_half = force.compute(ion, T / 2.0, ctx);
        
        // Force should oscillate
        REQUIRE(std::fabs(F0.x) > 1e-20);  // Non-zero at t=0
        REQUIRE(std::fabs(F_quarter.x) < 1e-18);  // Zero at quarter period
        REQUIRE(F_half.x == Approx(-F0.x));  // Opposite at half period
    }
    
    SECTION("LQIT with DC offset") {
        AnalyticalFieldParams params_dc = params;
        params_dc.dc_quad_voltage_V = 100.0;  // 100 V DC offset
        
        ElectricFieldForce force_dc(params_dc);
        
        IonState ion = make_test_ion(0.001, 0.0, 0.0);
        ForceContext ctx;
        
        // At t=0: U_eff = DC + RF*cos(0) = 100 + 1000 = 1100 V
        Vec3 F = force_dc.compute(ion, 0.0, ctx);
        
        double r0_sq = params.radius_m * params.radius_m;
        double E_x = 2.0 * ion.pos.x * 1100.0 / r0_sq;
        
        REQUIRE(F.x == Approx(ELEM_CHARGE_C * E_x));
    }
    
    SECTION("LQIT with DC endcap field") {
        AnalyticalFieldParams params_axial = params;
        params_axial.dc_axial_voltage_V = 500.0;  // 500 V endcap potential
        params_axial.length_m = 0.1;  // 10 cm trap
        
        ElectricFieldForce force_axial(params_axial);
        ForceContext ctx;
        
        // Test center region: should be field-free
        IonState ion_center = make_test_ion(0.0, 0.0, 0.05);  // z = L/2 (center)
        Vec3 F_center = force_axial.compute(ion_center, 0.0, ctx);
        REQUIRE(std::fabs(F_center.z) < 1e-20);  // No field in center
        
        // Test left endcap (z < 0.1*L): should push right (E_z > 0)
        IonState ion_left = make_test_ion(0.0, 0.0, 0.005);  // z = 5 mm (in left 10%)
        Vec3 F_left = force_axial.compute(ion_left, 0.0, ctx);
        REQUIRE(F_left.z > 0.0);  // Pushes away from left endcap
        
        // Test right endcap (z > 0.9*L): should push left (E_z < 0)
        IonState ion_right = make_test_ion(0.0, 0.0, 0.095);  // z = 95 mm (in right 10%)
        Vec3 F_right = force_axial.compute(ion_right, 0.0, ctx);
        REQUIRE(F_right.z < 0.0);  // Pushes away from right endcap
    }
    
    SECTION("LQIT with AC field (fixed x-direction for v1.0)") {
        AnalyticalFieldParams params_ac;
        params_ac.instrument_type = InstrumentType::LQIT;
        params_ac.radius_m = 0.005;
        params_ac.length_m = 0.1;
        params_ac.rf_voltage_V = 0.0;  // No RF for clean AC test
        params_ac.ac_voltage_V = 100.0;  // 100 V AC
        params_ac.ac_frequency_Hz = 5e5;  // 500 kHz
        
        ElectricFieldForce force_ac(params_ac);
        
        IonState ion = make_test_ion(0.0, 0.0, 0.05);  // x=0 to isolate AC field
        ForceContext ctx;
        
        // At t=0: cos(0) = 1, AC field should be maximum in +x direction
        Vec3 F0 = force_ac.compute(ion, 0.0, ctx);
        double expected_E_ac = params_ac.ac_voltage_V / params_ac.radius_m;
        REQUIRE(F0.x == Approx(ELEM_CHARGE_C * expected_E_ac));
        
        // At t=T/4: cos(ωt) = 0, AC field should be zero
        double T = 1.0 / params_ac.ac_frequency_Hz;
        Vec3 F_quarter = force_ac.compute(ion, T / 4.0, ctx);
        REQUIRE(std::fabs(F_quarter.x) < 1e-18);
        
        // At t=T/2: cos(ωt) = -1, AC field should be reversed
        Vec3 F_half = force_ac.compute(ion, T / 2.0, ctx);
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
    AnalyticalFieldParams params;
    params.instrument_type = InstrumentType::FTICR;
    params.fticr_voltage_V = 10.0;  // 10 V trapping voltage
    params.fticr_char_length_m = 0.05;  // 5 cm characteristic length
    params.length_m = 0.1;  // 10 cm cell length
    
    ElectricFieldForce force(params);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(FTICR)");
    }
    
    SECTION("Radial field increases with radius") {
        IonState ion1 = make_test_ion(0.001, 0.0, 0.05);  // r = 1 mm
        IonState ion2 = make_test_ion(0.002, 0.0, 0.05);  // r = 2 mm
        
        ForceContext ctx;
        Vec3 F1 = force.compute(ion1, 0.0, ctx);
        Vec3 F2 = force.compute(ion2, 0.0, ctx);
        
        // E_x = factor * x, so E should double when x doubles
        REQUIRE(F2.x == Approx(2.0 * F1.x));
    }
    
    SECTION("Axial field restoring (harmonic)") {
        // Field should be zero at center, increase linearly with z
        IonState ion_center = make_test_ion(0.0, 0.0, 0.05);  // At z_center
        IonState ion_offset = make_test_ion(0.0, 0.0, 0.06);  // 1 cm off-center
        
        ForceContext ctx;
        Vec3 F_center = force.compute(ion_center, 0.0, ctx);
        Vec3 F_offset = force.compute(ion_offset, 0.0, ctx);
        
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
    AnalyticalFieldParams params;
    params.instrument_type = InstrumentType::Orbitrap;
    params.orbitrap_k = 10000.0;  // 10 kV/m² curvature
    params.orbitrap_r_char = 0.01;  // 1 cm characteristic radius
    params.length_m = 0.05;  // 5 cm trap length
    
    ElectricFieldForce force(params);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(Orbitrap)");
    }
    
    SECTION("Axial harmonic confinement") {
        // At trap center (z = L/2), axial force should be zero
        IonState ion_center = make_test_ion(0.005, 0.0, 0.025);
        
        ForceContext ctx;
        Vec3 F_center = force.compute(ion_center, 0.0, ctx);
        
        // z_center = 0.025 - 0.025 = 0
        REQUIRE(std::fabs(F_center.z) < 1e-18);
    }
    
    SECTION("Axial restoring force") {
        // Off-axis: should have linear restoring force
        IonState ion1 = make_test_ion(0.005, 0.0, 0.030);  // 0.5 cm above center
        IonState ion2 = make_test_ion(0.005, 0.0, 0.020);  // 0.5 cm below center
        
        ForceContext ctx;
        Vec3 F1 = force.compute(ion1, 0.0, ctx);
        Vec3 F2 = force.compute(ion2, 0.0, ctx);
        
        // Forces should be opposite (restoring)
        REQUIRE(F1.z == Approx(-F2.z));
        
        // Force proportional to displacement: F_z = -k * z_center
        double z_center1 = 0.030 - 0.025;
        double expected_F_z1 = -ELEM_CHARGE_C * params.orbitrap_k * z_center1;
        REQUIRE(F1.z == Approx(expected_F_z1));
    }
    
    SECTION("Radial field depends on radius") {
        // r_char = 0.01 m, so test at r < r_char and r > r_char
        IonState ion1 = make_test_ion(0.005, 0.0, 0.025);  // r = 5 mm < r_char
        IonState ion2 = make_test_ion(0.015, 0.0, 0.025);  // r = 15 mm > r_char
        
        ForceContext ctx;
        Vec3 F1 = force.compute(ion1, 0.0, ctx);
        Vec3 F2 = force.compute(ion2, 0.0, ctx);
        
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
    AnalyticalFieldParams params;
    params.instrument_type = InstrumentType::QuadrupoleRF;
    params.radius_m = 0.003;  // 3 mm
    params.rf_voltage_V = 500.0;
    params.rf_frequency_Hz = 5e5;  // 500 kHz
    
    ElectricFieldForce force(params);
    
    SECTION("Force name") {
        REQUIRE(force.name() == "ElectricField(QuadrupoleRF)");
    }
    
    SECTION("Behaves like LQIT") {
        // QuadrupoleRF uses same formula as LQIT
        IonState ion = make_test_ion(0.001, 0.001, 0.0);
        
        ForceContext ctx;
        Vec3 F = force.compute(ion, 0.0, ctx);
        
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
        Vec3 F = force.compute(ion, 0.0, ctx);
        
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
        Vec3 F_pos = force.compute(ion_pos, 0.0, ctx);
        Vec3 F_neg = force.compute(ion_neg, 0.0, ctx);
        
        REQUIRE(F_neg.x == Approx(-F_pos.x));
    }
    
    SECTION("NoFixedInstrument returns zero field") {
        AnalyticalFieldParams params;
        params.instrument_type = InstrumentType::NoFixedInstrument;
        
        ElectricFieldForce force(params);
        
        IonState ion = make_test_ion(1, 2, 3);
        ForceContext ctx;
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        REQUIRE(F.x == Approx(0.0).margin(1e-25));
        REQUIRE(F.y == Approx(0.0).margin(1e-25));
        REQUIRE(F.z == Approx(0.0).margin(1e-25));
    }
}

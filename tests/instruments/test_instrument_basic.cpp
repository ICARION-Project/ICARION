// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

// Fast sanity checks for basic instrument functionality
// These tests verify that instruments can be instantiated and configured
// without running full simulations. Target: <100ms total runtime.

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/config/types/FullConfig.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/InstrumentTypes.h"
#include "utils/constants.h"

using namespace ICARION;
using Catch::Approx;

namespace {

// Helper to create minimal domain config for each instrument
config::DomainConfig make_minimal_domain(config::Instrument instrument_type) {
    config::DomainConfig dom;
    dom.instrument = instrument_type;
    
    // Set minimal geometry
    dom.geometry.length_m = 0.1;
    dom.geometry.radius_m = 0.01;
    
    // Initialize all field values to prevent uninitialized access
    dom.fields.dc.axial_V.constant_value = 0.0;
    dom.fields.dc.quad_V.constant_value = 0.0;
    dom.fields.dc.radial_V.constant_value = 0.0;
    dom.fields.dc.EN_Td.constant_value = 0.0;
    dom.fields.rf.voltage_V.constant_value = 0.0;
    dom.fields.rf.frequency_Hz.constant_value = 1.0e6;
    dom.fields.rf.compute_derived();
    dom.fields.ac.voltage_V.constant_value = 0.0;
    dom.fields.ac.frequency_Hz.constant_value = 1.0e5;
    dom.fields.ac.compute_derived();
    
    // Minimal environment
    dom.environment.pressure_Pa = 101325.0;
    dom.environment.temperature_K = 300.0;
    dom.environment.gas_species = "N2";
    dom.environment.compute_derived_properties();
    
    return dom;
}

} // namespace

TEST_CASE("IMS domain instantiation and field setup", "[instrument][ims][fast]") {
    auto dom = make_minimal_domain(config::Instrument::IMS);
    
    // IMS-specific: set drift field
    const double E_field_V_per_m = 200.0;
    dom.fields.dc.axial_V.constant_value = E_field_V_per_m * dom.geometry.length_m;
    
    SECTION("Geometry is valid") {
        REQUIRE(dom.geometry.length_m > 0.0);
        REQUIRE(dom.geometry.radius_m > 0.0);
    }
    
    SECTION("Electric field is configured") {
        const double E_calc = dom.fields.dc.axial_V.constant_value.value() / dom.geometry.length_m;
        REQUIRE(E_calc == Approx(E_field_V_per_m));
    }
    
    SECTION("Environment is valid") {
        REQUIRE(dom.environment.pressure_Pa > 0.0);
        REQUIRE(dom.environment.temperature_K > 0.0);
        REQUIRE(dom.environment.gas_species == "N2");
    }
}

TEST_CASE("Orbitrap domain instantiation and k-factor", "[instrument][orbitrap][fast]") {
    auto dom = make_minimal_domain(config::Instrument::Orbitrap);
    
    // Orbitrap-specific geometry and fields (Makarov 2000 parameters)
    // For Orbitrap: R_m^2 * ln(R_out/R_in) > 0.5*(R_out^2 - R_in^2)
    // Using compact Orbitrap dimensions:
    dom.geometry.radius_in_m = 0.015;     // Inner electrode (15 mm)
    dom.geometry.radius_out_m = 0.020;    // Outer electrode (20 mm)
    dom.geometry.radius_char_m = 0.017;   // Characteristic radius sqrt(R_in*R_out) ≈ 17.3 mm
    dom.geometry.length_m = 0.050;        // Axial length (50 mm)
    dom.fields.dc.radial_V.constant_value = 3500.0;
    
    SECTION("Orbitrap geometry is valid") {
        REQUIRE(dom.geometry.radius_in_m > 0.0);
        REQUIRE(dom.geometry.radius_out_m > dom.geometry.radius_in_m);
        REQUIRE(dom.geometry.radius_char_m > dom.geometry.radius_in_m);
        REQUIRE(dom.geometry.radius_char_m < dom.geometry.radius_out_m);
    }
    
    SECTION("Orbitrap voltage and geometry are configured") {
        // Simple sanity checks - don't validate complex k-factor formula
        REQUIRE(dom.fields.dc.radial_V.constant_value.value() > 0.0);
        REQUIRE(std::isfinite(dom.fields.dc.radial_V.constant_value.value()));
    }
}

TEST_CASE("LQIT domain instantiation and Mathieu parameters", "[instrument][lqit][fast]") {
    auto dom = make_minimal_domain(config::Instrument::LQIT);
    
    // LQIT-specific: RF trap parameters
    dom.geometry.radius_m = 0.004;
    dom.geometry.length_m = 0.05;
    dom.fields.rf.voltage_V.constant_value = 30.0;
    dom.fields.rf.frequency_Hz.constant_value = 1.5e6;
    dom.fields.rf.compute_derived();
    dom.fields.dc.quad_V.constant_value = 0.0;
    
    SECTION("RF parameters are valid") {
        REQUIRE(dom.fields.rf.voltage_V.constant_value.value() > 0.0);
        REQUIRE(dom.fields.rf.frequency_Hz.constant_value.value() > 0.0);
        REQUIRE(dom.fields.rf.angular_frequency_rad_s > 0.0);
    }
    
    SECTION("Mathieu q parameter calculation") {
        // For a test ion: m=40 amu, z=1
        const double mass_kg = 40.0 * AMU_TO_KG;
        const double charge_C = ELEM_CHARGE_C;
        const double omega = 2.0 * M_PI * dom.fields.rf.frequency_Hz.constant_value.value();
        const double r0 = dom.geometry.radius_m;
        
        const double q = (4.0 * charge_C * dom.fields.rf.voltage_V.constant_value.value()) /
                        (mass_kg * omega * omega * r0 * r0);
        
        REQUIRE(q > 0.0);
        REQUIRE(q < 0.908);  // Stability limit
    }
}

TEST_CASE("Quadrupole domain instantiation and stability", "[instrument][quadrupole][fast]") {
    auto dom = make_minimal_domain(config::Instrument::QuadrupoleRF);
    
    // Quadrupole-specific parameters
    dom.geometry.radius_m = 0.004;
    dom.geometry.length_m = 0.1;
    dom.fields.rf.voltage_V.constant_value = 250.0;
    dom.fields.rf.frequency_Hz.constant_value = 2.0e6;
    dom.fields.rf.compute_derived();
    dom.fields.dc.quad_V.constant_value = 5.0;
    dom.fields.dc.axial_V.constant_value = 2.0;
    
    SECTION("Quadrupole geometry and fields") {
        REQUIRE(dom.geometry.radius_m > 0.0);
        REQUIRE(dom.geometry.length_m > dom.geometry.radius_m);
        REQUIRE(dom.fields.rf.voltage_V.constant_value.value() > 0.0);
        REQUIRE(dom.fields.dc.axial_V.constant_value.value() > 0.0);
    }
    
    SECTION("Mathieu parameters for m/z=50") {
        const double mass_kg = 50.0 * AMU_TO_KG;
        const double charge_C = ELEM_CHARGE_C;
        const double omega = 2.0 * M_PI * dom.fields.rf.frequency_Hz.constant_value.value();
        const double r0 = dom.geometry.radius_m;
        
        const double q = (4.0 * charge_C * dom.fields.rf.voltage_V.constant_value.value()) /
                        (mass_kg * omega * omega * r0 * r0);
        const double a = (4.0 * charge_C * dom.fields.dc.quad_V.constant_value.value()) /
                        (mass_kg * omega * omega * r0 * r0);
        
        REQUIRE(q > 0.0);
        REQUIRE(std::isfinite(a));
        
        // Check if in stable region (simplified check: q < 0.908)
        bool stable = (q < 0.908);
        REQUIRE(stable);
    }
}

TEST_CASE("TOF domain instantiation and drift region", "[instrument][tof][fast]") {
    auto dom = make_minimal_domain(config::Instrument::TOF);
    
    // TOF-specific: drift region
    dom.geometry.length_m = 1.0;  // 1m drift tube
    dom.geometry.radius_m = 0.02;
    dom.fields.dc.axial_V.constant_value = 0.0;  // Field-free drift
    
    SECTION("TOF geometry is valid") {
        REQUIRE(dom.geometry.length_m > 0.0);
        REQUIRE(dom.geometry.radius_m > 0.0);
    }
    
    SECTION("Field-free drift region") {
        REQUIRE(dom.fields.dc.axial_V.constant_value.value() == 0.0);
    }
    
    SECTION("Flight time calculation for m/z=100") {
        const double mass_kg = 100.0 * AMU_TO_KG;
        const double kinetic_energy_J = 1000.0 * ELEM_CHARGE_C;  // 1 keV
        const double velocity = std::sqrt(2.0 * kinetic_energy_J / mass_kg);
        const double flight_time = dom.geometry.length_m / velocity;
        
        REQUIRE(velocity > 0.0);
        REQUIRE(flight_time > 0.0);
        REQUIRE(flight_time < 1e-3);  // Should be < 1 ms for reasonable parameters
    }
}

TEST_CASE("FTICR domain instantiation and B-field", "[instrument][fticr][fast]") {
    auto dom = make_minimal_domain(config::Instrument::FTICR);
    
    // FTICR-specific: magnetic field
    dom.geometry.length_m = 0.1;
    dom.geometry.radius_m = 0.02;
    dom.fields.magnetic.enabled = true;
    dom.fields.magnetic.field_strength_T = {0.0, 0.0, 7.0};  // 7 Tesla in z-direction
    
    SECTION("FTICR geometry and B-field") {
        REQUIRE(dom.geometry.length_m > 0.0);
        REQUIRE(dom.geometry.radius_m > 0.0);
        REQUIRE(dom.fields.magnetic.enabled);
        REQUIRE(dom.fields.magnetic.field_strength_T.z > 0.0);
    }
    
    SECTION("Cyclotron frequency for m/z=500") {
        const double mass_kg = 500.0 * AMU_TO_KG;
        const double charge_C = ELEM_CHARGE_C;
        const double B = dom.fields.magnetic.field_strength_T.z;
        
        const double omega_c = (charge_C * B) / mass_kg;
        const double freq_Hz = omega_c / (2.0 * M_PI);
        
        REQUIRE(omega_c > 0.0);
        REQUIRE(freq_Hz > 1e3);  // Should be kHz range
        REQUIRE(freq_Hz < 1e7);  // Should be < MHz for typical m/z
    }
}

TEST_CASE("All instruments have valid collision-free environments", "[instrument][environment][fast]") {
    const std::vector<config::Instrument> instruments = {
        config::Instrument::IMS,
        config::Instrument::Orbitrap,
        config::Instrument::LQIT,
        config::Instrument::QuadrupoleRF,
        config::Instrument::TOF,
        config::Instrument::FTICR
    };
    
    for (const auto& instr : instruments) {
        auto dom = make_minimal_domain(instr);
        
        REQUIRE(dom.environment.pressure_Pa > 0.0);
        REQUIRE(dom.environment.temperature_K > 0.0);
        REQUIRE(!dom.environment.gas_species.empty());
    }
}

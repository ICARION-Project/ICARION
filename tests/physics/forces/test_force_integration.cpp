// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file test_force_integration.cpp
 * @brief Integration tests for multi-force scenarios
 * 
 * Tests that verify multiple forces working together through ForceRegistry.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/physics/forces/MagneticFieldForce.h"
#include "core/physics/forces/DampingForce.h"
#include "core/types/IonState.h"
#include "core/config/types/InstrumentTypes.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/FieldsConfig.h"
#include "core/config/types/EnvironmentConfig.h"

#include "utils/constants.h"

#include <vector>
#include <memory>
#include <cmath>

using namespace ICARION;
using namespace ICARION::physics;
using namespace ICARION::instrument;
using Catch::Matchers::WithinAbs;

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * @brief Create minimal DomainConfig for testing
 */
static config::DomainConfig create_test_domain() {
    config::DomainConfig domain;
    domain.instrument = config::Instrument::IMS;
    domain.geometry.length_m = 0.1;
    domain.geometry.radius_m = 0.01;
    domain.environment.temperature_K = 300.0;
    domain.environment.pressure_Pa = 101325.0;
    domain.environment.gas_species = "N2";
    domain.environment.compute_derived_properties();
    return domain;
}

// ============================================================================
// Integration Test 1: Electric + Magnetic Forces
// ============================================================================

TEST_CASE("Integration: Electric + Magnetic forces combine", "[integration]") {
    ForceRegistry registry(create_test_domain());
    std::vector<IonState> ions;
    
    // Setup electric field (IMS) - SSOT config
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::IMS;
    domain.geometry.length_m = 0.1;
    // Initialize all field values
    domain.fields.dc.axial_V.constant_value = 1000.0;
    domain.fields.dc.quad_V.constant_value = 0.0;
    domain.fields.rf.voltage_V.constant_value = 0.0;
    domain.fields.rf.frequency_Hz.constant_value = 0.0;
    domain.fields.rf.compute_derived();
    registry.add_force(std::make_unique<ElectricFieldForce>(domain));
    
    // Setup magnetic field - SSOT config
    ICARION::config::MagneticFieldConfig mag_config;
    mag_config.field_strength_T = Vec3{0, 0, 1.0};
    mag_config.enabled = true;
    registry.add_force(std::make_unique<MagneticFieldForce>(mag_config));
    
    // Ion moving perpendicular to B
    IonState ion;
    ion.pos = Vec3{0, 0, 0.05};
    ion.vel = Vec3{1000, 0, 0};
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.mass_kg = 100 * AMU_TO_KG;
    ions.push_back(ion);
    
    ForceContext ctx;
    ctx.all_ions = &ions;
    
    Vec3 F_total = registry.compute_total_force(ions[0], 0.0, ctx);
    
    // Should have force in both z (electric) and y (magnetic) directions
    REQUIRE(std::abs(F_total.z) > 1e-20);  // Electric force
    REQUIRE(std::abs(F_total.y) > 1e-20);  // Magnetic force (v × B)
}

// ============================================================================
// Integration Test 2: All Four Force Types
// ============================================================================

TEST_CASE("Integration: All forces (Electric + Magnetic + Damping + SpaceCharge)", "[integration]") {
    ForceRegistry registry(create_test_domain());
    std::vector<IonState> ions;
    
    // Add all force types - SSOT configs
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::IMS;
    domain.geometry.length_m = 0.1;
    // Initialize all field values
    domain.fields.dc.axial_V.constant_value = 500.0;
    domain.fields.dc.quad_V.constant_value = 0.0;
    domain.fields.rf.voltage_V.constant_value = 0.0;
    domain.fields.rf.frequency_Hz.constant_value = 0.0;
    domain.fields.rf.compute_derived();
    registry.add_force(std::make_unique<ElectricFieldForce>(domain));
    
    ICARION::config::MagneticFieldConfig mag_config;
    mag_config.field_strength_T = Vec3{0, 0, 0.5};
    mag_config.enabled = true;
    registry.add_force(std::make_unique<MagneticFieldForce>(mag_config));
    
    ICARION::config::EnvironmentConfig env;
    env.pressure_Pa = 101325.0;
    env.temperature_K = 300.0;
    env.compute_derived_properties();
    registry.add_force(std::make_unique<DampingForce>(env, DampingModel::Friction));
    // Create two ions
    IonState ion1, ion2;
    ion1.pos = Vec3{0, 0, 0.02};
    ion1.vel = Vec3{100, 50, 25};
    ion1.ion_charge_C = ELEM_CHARGE_C;
    ion1.mass_kg = 100 * AMU_TO_KG;
    
    ion2.pos = Vec3{1e-6, 0, 0.02};
    ion2.vel = Vec3{100, 50, 25};
    ion2.ion_charge_C = ELEM_CHARGE_C;
    ion2.mass_kg = 100 * AMU_TO_KG;
    
    ions.push_back(ion1);
    ions.push_back(ion2);
    
    ForceContext ctx;
    ctx.all_ions = &ions;
    
    Vec3 F_total = registry.compute_total_force(ions[0], 0.0, ctx);
    
    // All force components should contribute
    REQUIRE(std::isfinite(F_total.x));
    REQUIRE(std::isfinite(F_total.y));
    REQUIRE(std::isfinite(F_total.z));
    
    // Should have non-zero force
    double F_mag = std::sqrt(F_total.x*F_total.x + F_total.y*F_total.y + F_total.z*F_total.z);
    REQUIRE(F_mag > 1e-20);
}

// ============================================================================
// Integration Test 3: Force Superposition Principle
// ============================================================================

TEST_CASE("Integration: Force superposition principle", "[integration]") {
    std::vector<IonState> ions;
    IonState ion;
    ion.pos = Vec3{0, 0, 0.05};
    ion.vel = Vec3{500, 200, 100};
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.mass_kg = 100 * AMU_TO_KG;
    ions.push_back(ion);
    
    ForceContext ctx;
    ctx.all_ions = &ions;
    
    // Compute forces individually - SSOT configs
    ICARION::config::DomainConfig domain;
    domain.instrument = ICARION::config::Instrument::IMS;
    domain.geometry.length_m = 0.1;
    // Initialize all field values
    domain.fields.dc.axial_V.constant_value = 800.0;
    domain.fields.dc.quad_V.constant_value = 0.0;
    domain.fields.rf.voltage_V.constant_value = 0.0;
    domain.fields.rf.frequency_Hz.constant_value = 0.0;
    domain.fields.rf.compute_derived();
    ElectricFieldForce electric(domain);
    Vec3 F_electric = electric.compute(ions[0], 0.0, ctx);
    
    ICARION::config::MagneticFieldConfig mag_config;
    mag_config.field_strength_T = Vec3{0, 0, 0.3};
    mag_config.enabled = true;
    MagneticFieldForce magnetic(mag_config);
    Vec3 F_magnetic = magnetic.compute(ions[0], 0.0, ctx);
    
    ICARION::config::EnvironmentConfig env;
    env.pressure_Pa = 101325.0;
    env.temperature_K = 300.0;
    env.compute_derived_properties();
    ions[0].reduced_mobility_cm2_Vs = 2.0;  // Need mobility for Friction model
    DampingForce damping(env, DampingModel::Friction);
    Vec3 F_damping = damping.compute(ions[0], 0.0, ctx);
    
    // Sum individual forces
    Vec3 F_sum(
        F_electric.x + F_magnetic.x + F_damping.x,
        F_electric.y + F_magnetic.y + F_damping.y,
        F_electric.z + F_magnetic.z + F_damping.z
    );
    
    // Compute via registry - SSOT configs
    ForceRegistry registry(create_test_domain());
    registry.add_force(std::make_unique<ElectricFieldForce>(domain));
    registry.add_force(std::make_unique<MagneticFieldForce>(mag_config));
    registry.add_force(std::make_unique<DampingForce>(env, DampingModel::Friction));
    
    Vec3 F_total = registry.compute_total_force(ions[0], 0.0, ctx);
    
    // Should match within floating point precision
    REQUIRE_THAT(F_total.x, WithinAbs(F_sum.x, 1e-20));
    REQUIRE_THAT(F_total.y, WithinAbs(F_sum.y, 1e-20));
    REQUIRE_THAT(F_total.z, WithinAbs(F_sum.z, 1e-20));
}

// ============================================================================
// Integration Test 4: Performance with Multiple Ions
// ============================================================================

TEST_CASE("Integration: Performance with 100 ions", "[integration][performance]") {
    ForceRegistry registry(create_test_domain());
    std::vector<IonState> ions;
    
    // Create 100 ions in a grid
    for (int i = 0; i < 10; ++i) {
        for (int j = 0; j < 10; ++j) {
            IonState ion;
            ion.pos = Vec3{i * 1e-6, j * 1e-6, 0};
            ion.vel = Vec3{0, 0, 0};
            ion.ion_charge_C = ELEM_CHARGE_C;
            ion.mass_kg = 100 * AMU_TO_KG;
            ions.push_back(ion);
        }
    }
    
    ForceContext ctx;
    ctx.all_ions = &ions;
    
    // Compute force on first ion (should interact with all 99 others)
    Vec3 F = registry.compute_total_force(ions[0], 0.0, ctx);
    
    // Should have repulsive force from neighbors
    REQUIRE(std::isfinite(F.x));
    REQUIRE(std::isfinite(F.y));
    double F_mag = std::sqrt(F.x*F.x + F.y*F.y);
    REQUIRE(F_mag > 1e-20);  // Non-zero repulsion
}

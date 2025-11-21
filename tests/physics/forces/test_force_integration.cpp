// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

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
#include "core/physics/forces/SpaceChargeForce.h"
#include "core/types/IonState.h"
#include "instrument/InstrumentTypes.h"

#include "utils/constants.h"

#include <vector>
#include <memory>
#include <cmath>

using namespace ICARION;
using namespace ICARION::physics;
using namespace ICARION::instrument;
using Catch::Matchers::WithinAbs;

// ============================================================================
// Integration Test 1: Electric + Magnetic Forces
// ============================================================================

TEST_CASE("Integration: Electric + Magnetic forces combine", "[integration]") {
    ForceRegistry registry;
    std::vector<IonState> ions;
    
    // Setup electric field (IMS)
    AnalyticalFieldParams e_params;
    e_params.instrument_type = InstrumentType::IMS;
    e_params.length_m = 0.1;
    e_params.dc_axial_voltage_V = 1000.0;
    registry.add_force(std::make_unique<ElectricFieldForce>(e_params));
    
    // Setup magnetic field
    MagneticFieldParams m_params;
    m_params.uniform_field_T = Vec3{0, 0, 1.0};
    m_params.enabled = true;
    registry.add_force(std::make_unique<MagneticFieldForce>(m_params));
    
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
    ForceRegistry registry;
    std::vector<IonState> ions;
    
    // Add all force types
    AnalyticalFieldParams e_params;
    e_params.instrument_type = InstrumentType::IMS;
    e_params.length_m = 0.1;
    e_params.dc_axial_voltage_V = 500.0;
    registry.add_force(std::make_unique<ElectricFieldForce>(e_params));
    
    MagneticFieldParams m_params;
    m_params.uniform_field_T = Vec3{0, 0, 0.5};
    m_params.enabled = true;
    registry.add_force(std::make_unique<MagneticFieldForce>(m_params));
    
    DampingParams d_params;
    d_params.model = DampingModel::Friction;
    d_params.gamma_coefficient = 1e5;
    registry.add_force(std::make_unique<DampingForce>(d_params));
    
    registry.add_force(std::make_unique<SpaceChargeForce>(1e-10));
    
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
    
    // Compute forces individually
    AnalyticalFieldParams e_params;
    e_params.instrument_type = InstrumentType::IMS;
    e_params.length_m = 0.1;
    e_params.dc_axial_voltage_V = 800.0;
    ElectricFieldForce electric(e_params);
    Vec3 F_electric = electric.compute(ions[0], 0.0, ctx);
    
    MagneticFieldParams m_params;
    m_params.uniform_field_T = Vec3{0, 0, 0.3};
    m_params.enabled = true;
    MagneticFieldForce magnetic(m_params);
    Vec3 F_magnetic = magnetic.compute(ions[0], 0.0, ctx);
    
    DampingParams d_params;
    d_params.model = DampingModel::Friction;
    d_params.gamma_coefficient = 2e5;
    DampingForce damping(d_params);
    Vec3 F_damping = damping.compute(ions[0], 0.0, ctx);
    
    // Sum individual forces
    Vec3 F_sum(
        F_electric.x + F_magnetic.x + F_damping.x,
        F_electric.y + F_magnetic.y + F_damping.y,
        F_electric.z + F_magnetic.z + F_damping.z
    );
    
    // Compute via registry
    ForceRegistry registry;
    registry.add_force(std::make_unique<ElectricFieldForce>(e_params));
    registry.add_force(std::make_unique<MagneticFieldForce>(m_params));
    registry.add_force(std::make_unique<DampingForce>(d_params));
    
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
    ForceRegistry registry;
    std::vector<IonState> ions;
    
    // Add space charge (O(N²) complexity)
    registry.add_force(std::make_unique<SpaceChargeForce>(1e-10));
    
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

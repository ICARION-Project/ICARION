// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file test_space_charge_integration.cpp
 * @brief Integration tests for complete space charge pipeline (CIC + Poisson + Field)
 * 
 * Tests the full workflow:
 * 1. Ion distribution → CIC charge deposition → Grid charge density
 * 2. Grid charge density → Poisson solver → Potential field
 * 3. Potential field → E-field gradient → Forces on ions
 * 
 * Verifies:
 * - Direct (exact) vs Grid (approximate) methods agree within tolerance
 * - CIC smoothness improves Poisson accuracy vs NGP
 * - Two-ion Coulomb repulsion follows analytical expectations
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/physics/spacecharge/spaceChargeSolver.h"
#include "core/physics/forces/SpaceChargeDirect.h"
#include "core/physics/forces/SpaceChargeGrid.h"
#include "core/types/IonState.h"
#include "core/physics/forces/ForceContext.h"
#include "utils/constants.h"

#include "core/physics/spacecharge/SpaceChargeGridModel.h"
#include "core/physics/spacecharge/SpaceChargeDirectModel.h"
#include "core/config/types/CylindricalGeometry.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonEnsemble.h"
#include "core/utils/mathUtils.h"

#include <vector>
#include <cmath>
#include <limits>
#include <memory>

using namespace ICARION;
using namespace ICARION::physics;
using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

// ====================================================================================
// TEST CASE 0: Geometry-driven grid parity using new models
// ====================================================================================

TEST_CASE("SpaceCharge Integration: geometry-aware grid vs direct model",
          "[spacecharge][integration][geometry]") {
    config::DomainConfig dom;
    dom.name = "rotated_cylinder";
    dom.geometry.length_m = 0.06;
    dom.geometry.radius_m = 0.008;
    dom.geometry.end_aperture_m = 0.008;
    dom.geometry.origin_m = Vec3{0.012, -0.018, 0.01};

    const double theta = M_PI / 4.0; // 45 degrees
    Mat3 R;
    R.m[0][0] = std::cos(theta); R.m[0][1] = -std::sin(theta); R.m[0][2] = 0.0;
    R.m[1][0] = std::sin(theta); R.m[1][1] =  std::cos(theta); R.m[1][2] = 0.0;
    R.m[2][0] = 0.0;             R.m[2][1] = 0.0;              R.m[2][2] = 1.0;
    dom.rotation_global_to_local = R;
    dom.rotation_local_to_global = transpose(R);

    config::CylindricalGeometry geometry(dom);

    std::vector<IonState> ions;
    ions.reserve(32);
    for (int ix = 0; ix < 2; ++ix) {
        for (int iy = 0; iy < 2; ++iy) {
            for (int iz = 0; iz < 4; ++iz) {
                Vec3 local_pos{
                    -0.003 + ix * 0.0025,
                    -0.002 + iy * 0.0025,
                    0.002 + iz * 0.01};
                IonState ion;
                ion.pos = geometry.local_to_global_pos(local_pos);
                ion.vel = Vec3{0.0, 0.0, 0.0};
                ion.ion_charge_C = ELEM_CHARGE_C;
                ion.mass_kg = AMU_TO_KG;
                ion.active = true;
                ion.born = true;
                ions.push_back(ion);
            }
        }
    }

    core::IonEnsemble ensemble = core::IonEnsemble::from_legacy(ions);

    physics::SpaceChargeDirectModel direct_model(1e-10);
    direct_model.update_fields(ensemble, 0.0);

    auto geometry_model = std::make_unique<config::CylindricalGeometry>(dom);
    physics::SpaceChargeGridModel grid_model(dom, std::move(geometry_model), 0.0005, 64);
    grid_model.update_fields(ensemble, 0.0);

    constexpr double SIGNIFICANCE = 1e-3;
    double sum_rel_sq = 0.0;
    size_t rel_count = 0;
    double max_rel_significant = 0.0;
    double max_abs_small = 0.0;

    for (size_t i = 0; i < ensemble.size(); ++i) {
        Vec3 E_direct = direct_model.sample_electric_field(i);
        Vec3 E_grid = grid_model.sample_electric_field(i);
        Vec3 diff = E_grid - E_direct;

        double mag_direct = norm(E_direct);
        double mag_diff = norm(diff);

        if (mag_direct > SIGNIFICANCE) {
            double rel = mag_diff / std::max(1e-12, mag_direct);
            max_rel_significant = std::max(max_rel_significant, rel);
            sum_rel_sq += rel * rel;
            rel_count++;
        } else {
            max_abs_small = std::max(max_abs_small, mag_diff);
        }
    }
    double rms_rel = rel_count > 0 ? std::sqrt(sum_rel_sq / rel_count) : 0.0;

    INFO("Geometry-aware grid parity rms_rel=" << rms_rel
         << " max_rel(sig)=" << max_rel_significant
         << " max_abs_small=" << max_abs_small
         << " rel_count=" << rel_count);

    REQUIRE(max_rel_significant < 0.7);
    REQUIRE(max_abs_small < 2e-3);
    REQUIRE(rms_rel < 0.35);
}

// ====================================================================================
// TEST CASE 1: Two-ion Coulomb repulsion (Direct vs Grid comparison)
// ====================================================================================

TEST_CASE("SpaceCharge Integration: Two-ion Coulomb repulsion (Direct vs Grid)", "[spacecharge][integration]") {
    SECTION("Two protons 1mm apart") {
        // Setup: Two protons separated by 1mm along x-axis
        std::vector<IonState> ions(2);
        
        ions[0].pos = Vec3{-0.0005, 0.0, 0.0};  // -0.5mm
        ions[1].pos = Vec3{ 0.0005, 0.0, 0.0};  //  0.5mm
        
        ions[0].vel = Vec3{0.0, 0.0, 0.0};
        ions[1].vel = Vec3{0.0, 0.0, 0.0};
        
        ions[0].ion_charge_C = ELEM_CHARGE_C;  // 1.602e-19 C
        ions[1].ion_charge_C = ELEM_CHARGE_C;
        
        ions[0].mass_kg = AMU_TO_KG;
        ions[1].mass_kg = AMU_TO_KG;
        
        ions[0].active = true;
        ions[1].active = true;
        ions[0].born = true;
        ions[1].born = true;
        
        // Analytical Coulomb force: F = k_e * q1 * q2 / r²
        double r = 0.001;  // 1mm separation
        double F_analytical = COULOMB_CONST * ELEM_CHARGE_C * ELEM_CHARGE_C / (r * r);
        
        INFO("Analytical Coulomb force: " << F_analytical << " N");
        INFO("Expected: ~2.3e-22 N for two protons at 1mm");
        
        // ----------------------------------------------------------------------
        // Method 1: Direct Coulomb (SpaceChargeDirect)
        // ----------------------------------------------------------------------
        
        SpaceChargeDirect direct_force(1e-10);  // 0.1nm softening
        ForceContext ctx_direct;
        ctx_direct.all_ions = &ions;
        
        Vec3 F_direct_ion0 = direct_force.compute(ions[0], 0.0, ctx_direct);
        Vec3 F_direct_ion1 = direct_force.compute(ions[1], 0.0, ctx_direct);
        
        INFO("Direct method:");
        INFO("  F_ion0 = (" << F_direct_ion0.x << ", " << F_direct_ion0.y << ", " << F_direct_ion0.z << ")");
        INFO("  F_ion1 = (" << F_direct_ion1.x << ", " << F_direct_ion1.y << ", " << F_direct_ion1.z << ")");
        
        // Ion 0 should be pushed LEFT (-x direction)
        REQUIRE(F_direct_ion0.x < 0.0);
        // Ion 1 should be pushed RIGHT (+x direction)
        REQUIRE(F_direct_ion1.x > 0.0);
        
        // Magnitudes should be equal (symmetry)
        double F_direct_mag0 = std::sqrt(F_direct_ion0.x * F_direct_ion0.x + 
                                        F_direct_ion0.y * F_direct_ion0.y + 
                                        F_direct_ion0.z * F_direct_ion0.z);
        double F_direct_mag1 = std::sqrt(F_direct_ion1.x * F_direct_ion1.x + 
                                        F_direct_ion1.y * F_direct_ion1.y + 
                                        F_direct_ion1.z * F_direct_ion1.z);
        
        REQUIRE_THAT(F_direct_mag0, WithinRel(F_direct_mag1, 1e-10));
        
        // Should match analytical result (Direct is exact)
        REQUIRE_THAT(F_direct_mag0, WithinRel(F_analytical, 0.01));  // 1% tolerance (softening)
        
        // ----------------------------------------------------------------------
        // Method 2: Grid Poisson (SpaceChargeGrid)
        // ----------------------------------------------------------------------
        
        // Grid setup: 64³ cells, 4mm domain (±2mm in each direction)
        // Larger domain with good resolution for 1mm ion separation
        int nx = 64, ny = 64, nz = 64;
        double domain_size = 0.004;  // 4mm
        double dx = domain_size / nx;  // 62.5μm cell size
        double dy = domain_size / ny;
        double dz = domain_size / nz;
        Vec3 grid_origin{-0.002, -0.002, -0.002};  // Center at (0,0,0)
        
        auto solver = std::make_shared<SpaceChargeSolver>(nx, ny, nz, dx, dy, dz, grid_origin);
        SpaceChargeGrid grid_force(solver);
        
        ForceContext ctx_grid;
        ctx_grid.all_ions = &ions;
        
        Vec3 F_grid_ion0 = grid_force.compute(ions[0], 0.0, ctx_grid);
        Vec3 F_grid_ion1 = grid_force.compute(ions[1], 0.0, ctx_grid);
        
        INFO("Grid method:");
        INFO("  F_ion0 = (" << F_grid_ion0.x << ", " << F_grid_ion0.y << ", " << F_grid_ion0.z << ")");
        INFO("  F_ion1 = (" << F_grid_ion1.x << ", " << F_grid_ion1.y << ", " << F_grid_ion1.z << ")");
        
        // Grid should also show repulsion
        REQUIRE(F_grid_ion0.x < 0.0);
        REQUIRE(F_grid_ion1.x > 0.0);
        
        double F_grid_mag0 = std::sqrt(F_grid_ion0.x * F_grid_ion0.x + 
                                      F_grid_ion0.y * F_grid_ion0.y + 
                                      F_grid_ion0.z * F_grid_ion0.z);
        double F_grid_mag1 = std::sqrt(F_grid_ion1.x * F_grid_ion1.x + 
                                      F_grid_ion1.y * F_grid_ion1.y + 
                                      F_grid_ion1.z * F_grid_ion1.z);
        
        INFO("Grid vs Analytical error: " << std::abs(F_grid_mag0 - F_analytical) / F_analytical * 100 << "%");
        
        // Grid method should be within 30% of analytical for 2-ion case
        // Note: Grid methods have larger errors for low charge densities (N=2)
        // Accuracy improves with more ions (N>100 typically <10% error)
        REQUIRE_THAT(F_grid_mag0, WithinRel(F_analytical, 0.30));
        
        // ----------------------------------------------------------------------
        // Direct vs Grid Comparison
        // ----------------------------------------------------------------------
        
        double direct_vs_grid_error = std::abs(F_direct_mag0 - F_grid_mag0) / F_direct_mag0;
        INFO("Direct vs Grid relative error: " << direct_vs_grid_error * 100 << "%");
        
        // Direct and Grid should agree within 30% for 2-ion case
        // Known limitation: Grid accuracy improves significantly with more ions (N>100 → <10%)
        REQUIRE(direct_vs_grid_error < 0.30);
    }
    
    SECTION("Ten ions in random cluster") {
        // More complex case: 10 ions in 1mm³ volume
        std::vector<IonState> ions(10);
        
        // Random positions in 1mm³ cube
        std::vector<Vec3> positions = {
            {-0.0004, -0.0003,  0.0002}, { 0.0003, -0.0004, -0.0001},
            {-0.0002,  0.0004,  0.0003}, { 0.0001,  0.0002, -0.0004},
            { 0.0004,  0.0001,  0.0001}, {-0.0003,  0.0000,  0.0004},
            { 0.0000, -0.0002, -0.0003}, { 0.0002,  0.0003,  0.0002},
            {-0.0001, -0.0001,  0.0000}, { 0.0003,  0.0004,  0.0003}
        };
        
        for (size_t i = 0; i < ions.size(); ++i) {
            ions[i].pos = positions[i];
            ions[i].vel = Vec3{0.0, 0.0, 0.0};
            ions[i].ion_charge_C = ELEM_CHARGE_C;
            ions[i].mass_kg = AMU_TO_KG;
            ions[i].active = true;
            ions[i].born = true;
        }
        
        // Direct method
        SpaceChargeDirect direct_force(1e-10);
        ForceContext ctx_direct;
        ctx_direct.all_ions = &ions;
        
        std::vector<Vec3> forces_direct;
        for (const auto& ion : ions) {
            forces_direct.push_back(direct_force.compute(ion, 0.0, ctx_direct));
        }
        
        // Grid method
        auto solver = std::make_shared<SpaceChargeSolver>(64, 64, 64, 
                                                          3e-5, 3e-5, 3e-5,  // 30μm cells
                                                          Vec3{-0.001, -0.001, -0.001});
        SpaceChargeGrid grid_force(solver);
        ForceContext ctx_grid;
        ctx_grid.all_ions = &ions;
        
        std::vector<Vec3> forces_grid;
        for (const auto& ion : ions) {
            forces_grid.push_back(grid_force.compute(ion, 0.0, ctx_grid));
        }
        
        // Compare forces
        int agreements = 0;
        for (size_t i = 0; i < ions.size(); ++i) {
            double F_direct_mag = std::sqrt(forces_direct[i].x * forces_direct[i].x + 
                                           forces_direct[i].y * forces_direct[i].y + 
                                           forces_direct[i].z * forces_direct[i].z);
            double F_grid_mag = std::sqrt(forces_grid[i].x * forces_grid[i].x + 
                                         forces_grid[i].y * forces_grid[i].y + 
                                         forces_grid[i].z * forces_grid[i].z);
            
            if (F_direct_mag > 1e-24) {  // Skip negligible forces
                double rel_error = std::abs(F_grid_mag - F_direct_mag) / F_direct_mag;
                INFO("Ion " << i << ": F_direct=" << F_direct_mag << " N, F_grid=" << F_grid_mag << " N, error=" << rel_error * 100 << "%");
                
                if (rel_error < 0.30) {  // 30% tolerance for complex geometry
                    agreements++;
                }
            }
        }
        
        // At least 70% of ions should have good agreement
        INFO("Agreements: " << agreements << " / " << ions.size());
        REQUIRE(agreements >= 7);
    }
}

// ====================================================================================
// TEST CASE 2: E-field smoothness validation
// ====================================================================================

TEST_CASE("SpaceCharge Integration: E-field smoothness near ion", "[spacecharge][integration][field]") {
    SECTION("E-field smoothness near ion") {
        // Test: E-field should be smooth around the ion (no grid noise)
        
        IonState ion;
        ion.pos = Vec3{0.0, 0.0, 0.0};  // Center
        ion.vel = Vec3{0.0, 0.0, 0.0};
        ion.ion_charge_C = ELEM_CHARGE_C;
        ion.mass_kg = AMU_TO_KG;
        ion.active = true;
        ion.born = true;
        
        std::vector<IonState> ions = {ion};
        
        auto solver = std::make_shared<SpaceChargeSolver>(64, 64, 64, 
                                                          3e-5, 3e-5, 3e-5,  // 30μm cells
                                                          Vec3{-0.001, -0.001, -0.001});
        solver->update(ions);
        
        // Sample E-field along x-axis at different distances
        std::vector<double> distances = {0.0001, 0.0002, 0.0003, 0.0004, 0.0005};  // 0.1-0.5mm
        std::vector<double> E_magnitudes;
        
        for (double r : distances) {
            Vec3 E = solver->fieldAt(Vec3{r, 0.0, 0.0});
            double E_mag = std::sqrt(E.x * E.x + E.y * E.y + E.z * E.z);
            E_magnitudes.push_back(E_mag);
            
            INFO("r = " << r * 1000 << " mm: |E| = " << E_mag << " V/m");
        }
        
        // E-field should decrease monotonically (no oscillations from grid noise)
        for (size_t i = 1; i < E_magnitudes.size(); ++i) {
            REQUIRE(E_magnitudes[i] < E_magnitudes[i-1]);
        }
        
        // E-field should roughly follow 1/r² scaling
        for (size_t i = 1; i < E_magnitudes.size(); ++i) {
            double r_ratio = distances[i] / distances[i-1];
            double E_ratio = E_magnitudes[i-1] / E_magnitudes[i];
            double expected_ratio = r_ratio * r_ratio;
            
            INFO("r_ratio = " << r_ratio << ", E_ratio = " << E_ratio << ", expected = " << expected_ratio);
            
            // Should be within 30% of 1/r² (grid discretization + boundary effects)
            REQUIRE_THAT(E_ratio, WithinRel(expected_ratio, 0.30));
        }
    }
}

// ====================================================================================
// TEST CASE 3: Force symmetry and conservation
// ====================================================================================

TEST_CASE("SpaceCharge Integration: Force symmetry", "[spacecharge][integration]") {
    SECTION("Newton's third law: F_ij = -F_ji") {
        // Two ions should exert equal and opposite forces
        
        std::vector<IonState> ions(2);
        ions[0].pos = Vec3{-0.0003, 0.0, 0.0};
        ions[1].pos = Vec3{ 0.0003, 0.0, 0.0};
        
        ions[0].vel = Vec3{0.0, 0.0, 0.0};
        ions[1].vel = Vec3{0.0, 0.0, 0.0};
        
        ions[0].ion_charge_C = ELEM_CHARGE_C;
        ions[1].ion_charge_C = ELEM_CHARGE_C;
        
        ions[0].mass_kg = AMU_TO_KG;
        ions[1].mass_kg = AMU_TO_KG;
        
        ions[0].active = true;
        ions[1].active = true;
        ions[0].born = true;
        ions[1].born = true;
        
        // Direct method (should be exact)
        SpaceChargeDirect direct_force(1e-10);
        ForceContext ctx;
        ctx.all_ions = &ions;
        
        Vec3 F0 = direct_force.compute(ions[0], 0.0, ctx);
        Vec3 F1 = direct_force.compute(ions[1], 0.0, ctx);
        
        // Forces should be equal and opposite
        REQUIRE_THAT(F0.x, WithinAbs(-F1.x, 1e-30));
        REQUIRE_THAT(F0.y, WithinAbs(-F1.y, 1e-30));
        REQUIRE_THAT(F0.z, WithinAbs(-F1.z, 1e-30));
        
        // Net force should be zero (momentum conservation)
        Vec3 F_total = F0 + F1;
        REQUIRE_THAT(F_total.x, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(F_total.y, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(F_total.z, WithinAbs(0.0, 1e-30));
    }
}

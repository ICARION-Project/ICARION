// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file test_charge_deposition.cpp
 * @brief Unit tests for charge deposition (NGP method)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/spacecharge/depositCharge.h"
#include "core/types/Grid3D.h"
#include "core/types/IonState.h"

#include <cmath>
#include <vector>
#include <numeric>
#include <random>

using Catch::Approx;
using ICARION::core::IonState;
using ICARION::core::Vec3;

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * @brief Create ion cloud with random positions
 */
std::vector<IonState> create_random_ion_cloud(int N, double radius = 0.001, unsigned seed = 42) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<double> dist(-radius, radius);
    
    std::vector<IonState> ions;
    ions.reserve(N);
    
    for (int i = 0; i < N; ++i) {
        IonState ion;
        ion.active = true;
        ion.born = true;
        ion.ion_charge_C = 1.6e-19;  // Elementary charge
        ion.pos = {dist(rng), dist(rng), dist(rng)};
        ion.vel = {0, 0, 0};
        ion.t = 0.0;
        
        ions.push_back(ion);
    }
    
    return ions;
}

// ============================================================================
// TEST CASES
// ============================================================================

TEST_CASE("depositCharge: Charge conservation", "[spacecharge][deposition][conservation]") {
    // Total charge on grid should equal sum of ion charges
    
    const int N = 32;
    const double cell_size = 2e-4;  // 200 μm cells
    const double domain_size = N * cell_size;  // 6.4 mm domain
    
    Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
    grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    SECTION("100 ions, random positions") {
        auto ions = create_random_ion_cloud(100, 0.002);  // Within ±2mm
        
        // Compute total ion charge
        double Q_ions = 0.0;
        for (const auto& ion : ions) {
            if (ion.active && ion.born) {
                Q_ions += ion.ion_charge_C;
            }
        }
        
        // Deposit charges
        auto rho = deposit_charge(ions, grid);
        
        // Compute total charge on grid
        double Q_grid = 0.0;
        const double cell_volume = cell_size * cell_size * cell_size;
        for (double rho_cell : rho) {
            Q_grid += rho_cell * cell_volume;
        }
        
        INFO("Total ion charge: " << Q_ions << " C");
        INFO("Total grid charge: " << Q_grid << " C");
        INFO("Relative error: " << std::abs(Q_grid - Q_ions) / Q_ions);
        
        // Charge should be conserved to within 1%
        REQUIRE(Q_grid == Approx(Q_ions).epsilon(0.01));
    }
    
    SECTION("1000 ions, various charges") {
        auto ions = create_random_ion_cloud(1000, 0.002);
        
        // Mix of different charges
        for (size_t i = 0; i < ions.size(); ++i) {
            if (i % 3 == 0) ions[i].ion_charge_C = 1.6e-19;   // +1e
            else if (i % 3 == 1) ions[i].ion_charge_C = 3.2e-19;  // +2e
            else ions[i].ion_charge_C = -1.6e-19;  // -1e (electrons)
        }
        
        double Q_ions = std::accumulate(ions.begin(), ions.end(), 0.0,
            [](double sum, const IonState& ion) {
                return sum + (ion.active && ion.born ? ion.ion_charge_C : 0.0);
            });
        
        auto rho = deposit_charge(ions, grid);
        
        const double cell_volume = cell_size * cell_size * cell_size;
        double Q_grid = std::accumulate(rho.begin(), rho.end(), 0.0) * cell_volume;
        
        INFO("Total ion charge (mixed): " << Q_ions << " C");
        INFO("Total grid charge: " << Q_grid << " C");
        
        REQUIRE(Q_grid == Approx(Q_ions).epsilon(0.01));
    }
}

TEST_CASE("depositCharge: Single ion at grid point", "[spacecharge][deposition][ngp]") {
    // Ion exactly at grid point → all charge in one cell
    
    const int N = 10;
    const double cell_size = 1e-4;  // 100 μm cells
    const double domain_size = N * cell_size;
    
    Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
    grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    IonState ion;
    ion.active = true;
    ion.born = true;
    ion.ion_charge_C = 1.6e-19;
    ion.pos = {0, 0, 0};  // At grid origin (center of domain)
    
    auto rho = deposit_charge({ion}, grid);
    
    // Find cell with charge
    int center_idx = grid.index(N/2, N/2, N/2);
    const double cell_volume = cell_size * cell_size * cell_size;
    double expected_rho = ion.ion_charge_C / cell_volume;
    
    INFO("Expected charge density: " << expected_rho << " C/m³");
    INFO("Computed charge density at center: " << rho[center_idx] << " C/m³");
    
    // All charge should be in center cell
    REQUIRE(rho[center_idx] == Approx(expected_rho).epsilon(1e-10));
    
    // All other cells should be zero (or nearly zero)
    double sum_other = 0.0;
    for (size_t i = 0; i < rho.size(); ++i) {
        if (static_cast<int>(i) != center_idx) {
            sum_other += std::abs(rho[i]);
        }
    }
    
    REQUIRE(sum_other < 1e-30);  // Essentially zero
}

TEST_CASE("depositCharge: NGP method validation", "[spacecharge][deposition][ngp]") {
    // Test NGP (Nearest Grid Point) assignment
    
    const int N = 16;
    const double cell_size = 1e-4;
    const double domain_size = N * cell_size;
    
    Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
    grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    SECTION("Ion between grid points snaps to nearest") {
        // Place ion slightly off-center
        IonState ion;
        ion.active = true;
        ion.born = true;
        ion.ion_charge_C = 1.6e-19;
        ion.pos = {0.00003, 0, 0};  // 30 μm from center (closer to i=8 than i=7)
        
        auto rho = deposit_charge({ion}, grid);
        
        // Find which cell got the charge
        int count_nonzero = 0;
        int nonzero_idx = -1;
        for (size_t i = 0; i < rho.size(); ++i) {
            if (std::abs(rho[i]) > 1e-30) {
                count_nonzero++;
                nonzero_idx = static_cast<int>(i);
            }
        }
        
        INFO("Non-zero cells: " << count_nonzero);
        
        // NGP: Charge should be in exactly ONE cell
        REQUIRE(count_nonzero == 1);
        REQUIRE(nonzero_idx >= 0);
        
        // Charge conservation
        const double cell_volume = cell_size * cell_size * cell_size;
        double total_charge = rho[nonzero_idx] * cell_volume;
        REQUIRE(total_charge == Approx(ion.ion_charge_C).epsilon(1e-10));
    }
}

TEST_CASE("depositCharge: OpenMP thread safety", "[spacecharge][deposition][parallel]") {
    // Test with many ions + parallel deposition
    
    const int N = 64;
    const double cell_size = 2e-4;
    const double domain_size = N * cell_size;
    
    Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
    grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    auto ions = create_random_ion_cloud(10000, 0.005);  // ±5mm
    
    // Run multiple times - should give same result (deterministic)
    auto rho1 = deposit_charge(ions, grid);
    auto rho2 = deposit_charge(ions, grid);
    
    for (size_t i = 0; i < rho1.size(); ++i) {
        INFO("Cell " << i << ": rho1=" << rho1[i] << " rho2=" << rho2[i]);
        REQUIRE(rho2[i] == Approx(rho1[i]).margin(1e-20));  // Exact match expected
    }
}

TEST_CASE("depositCharge: Inactive ions excluded", "[spacecharge][deposition][active]") {
    // Only active ions should contribute
    
    const int N = 16;
    const double cell_size = 1e-4;
    const double domain_size = N * cell_size;
    
    Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
    grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    std::vector<IonState> ions;
    
    // Active ion
    IonState ion1;
    ion1.active = true;
    ion1.born = true;
    ion1.ion_charge_C = 1.6e-19;
    ion1.pos = {0, 0, 0};
    ions.push_back(ion1);
    
    // Inactive ion (should be ignored)
    IonState ion2;
    ion2.active = false;
    ion2.born = true;
    ion2.ion_charge_C = 1.6e-19;
    ion2.pos = {0, 0, 0};  // Same position
    ions.push_back(ion2);
    
    // Not born ion (should be ignored)
    IonState ion3;
    ion3.active = true;
    ion3.born = false;
    ion3.ion_charge_C = 1.6e-19;
    ion3.pos = {0, 0, 0};
    ions.push_back(ion3);
    
    auto rho = deposit_charge(ions, grid);
    
    // Compute total charge
    const double cell_volume = cell_size * cell_size * cell_size;
    double Q_grid = std::accumulate(rho.begin(), rho.end(), 0.0) * cell_volume;
    
    INFO("Total grid charge: " << Q_grid << " C");
    INFO("Expected (1 active ion): " << ion1.ion_charge_C << " C");
    
    // Only the active+born ion should contribute
    REQUIRE(Q_grid == Approx(ion1.ion_charge_C).epsilon(1e-10));
}

TEST_CASE("depositCharge: Boundary handling", "[spacecharge][deposition][boundary]") {
    // Ions outside grid should be ignored
    
    const int N = 16;
    const double cell_size = 1e-4;
    const double domain_size = N * cell_size;
    
    Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
    grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    std::vector<IonState> ions;
    
    // Ion inside grid
    IonState ion_inside;
    ion_inside.active = true;
    ion_inside.born = true;
    ion_inside.ion_charge_C = 1.6e-19;
    ion_inside.pos = {0, 0, 0};
    ions.push_back(ion_inside);
    
    // Ion outside grid (should be ignored)
    IonState ion_outside;
    ion_outside.active = true;
    ion_outside.born = true;
    ion_outside.ion_charge_C = 1.6e-19;
    ion_outside.pos = {0.01, 0, 0};  // 10mm (way outside 1.6mm domain)
    ions.push_back(ion_outside);
    
    auto rho = deposit_charge(ions, grid);
    
    const double cell_volume = cell_size * cell_size * cell_size;
    double Q_grid = std::accumulate(rho.begin(), rho.end(), 0.0) * cell_volume;
    
    INFO("Total grid charge: " << Q_grid << " C");
    INFO("Expected (1 inside ion): " << ion_inside.ion_charge_C << " C");
    
    // Only the inside ion should contribute
    REQUIRE(Q_grid == Approx(ion_inside.ion_charge_C).epsilon(1e-10));
}

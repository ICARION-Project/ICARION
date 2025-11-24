// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file test_poisson_solver.cpp
 * @brief Unit tests for PoissonSolver (legacy space charge code validation)
 * 
 * Tests the legacy Poisson solver against analytical solutions and
 * validates numerical properties (convergence, solver agreement).
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/spacecharge/poissonSolver.h"
#include "core/physics/spacecharge/depositCharge.h"
#include "core/types/Grid3D.h"
#include "core/types/IonState.h"
#include "utils/constants.h"

#include <cmath>
#include <vector>
#include <iostream>

using Catch::Approx;
using ICARION::core::IonState;
using ICARION::core::Vec3;

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * @brief Compute analytical potential for point charge at origin
 * @param r Distance from charge [m]
 * @param Q Charge [C]
 * @return Potential [V]
 */
double analytical_point_charge_potential(double r, double Q) {
    if (r < 1e-12) return 0.0;  // Avoid singularity at origin
    return Q / (4.0 * M_PI * EPSILON_0 * r);
}

/**
 * @brief Compute analytical potential for uniform charged sphere
 * @param r Distance from center [m]
 * @param R Sphere radius [m]
 * @param Q Total charge [C]
 * @return Potential [V]
 * 
 * Inside sphere (r < R): φ(r) = (Q/(8πε₀R)) * (3 - r²/R²)
 * Outside sphere (r ≥ R): φ(r) = Q/(4πε₀r)
 */
double analytical_sphere_potential(double r, double R, double Q) {
    const double k = 1.0 / (4.0 * M_PI * EPSILON_0);
    
    if (r < R) {
        // Inside sphere
        return (Q / (8.0 * M_PI * EPSILON_0 * R)) * (3.0 - (r*r) / (R*R));
    } else {
        // Outside sphere
        return k * Q / r;
    }
}

/**
 * @brief Interpolate potential at arbitrary position using trilinear interpolation
 */
double interpolate_potential(const Grid3D& grid, const Vec3& pos) {
    // Convert position to grid coordinates
    double x_rel = (pos.x - grid.origin_m.x) / grid.dx;
    double y_rel = (pos.y - grid.origin_m.y) / grid.dy;
    double z_rel = (pos.z - grid.origin_m.z) / grid.dz;
    
    // Get cell indices
    int i0 = static_cast<int>(std::floor(x_rel));
    int j0 = static_cast<int>(std::floor(y_rel));
    int k0 = static_cast<int>(std::floor(z_rel));
    
    // Clamp to grid bounds
    i0 = std::max(0, std::min(i0, grid.Nx - 2));
    j0 = std::max(0, std::min(j0, grid.Ny - 2));
    k0 = std::max(0, std::min(k0, grid.Nz - 2));
    
    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;
    
    // Fractional positions within cell
    double fx = x_rel - i0;
    double fy = y_rel - j0;
    double fz = z_rel - k0;
    
    fx = std::max(0.0, std::min(1.0, fx));
    fy = std::max(0.0, std::min(1.0, fy));
    fz = std::max(0.0, std::min(1.0, fz));
    
    // Trilinear interpolation
    double c000 = grid.phi[grid.index(i0, j0, k0)];
    double c100 = grid.phi[grid.index(i1, j0, k0)];
    double c010 = grid.phi[grid.index(i0, j1, k0)];
    double c110 = grid.phi[grid.index(i1, j1, k0)];
    double c001 = grid.phi[grid.index(i0, j0, k1)];
    double c101 = grid.phi[grid.index(i1, j0, k1)];
    double c011 = grid.phi[grid.index(i0, j1, k1)];
    double c111 = grid.phi[grid.index(i1, j1, k1)];
    
    double c00 = c000 * (1 - fx) + c100 * fx;
    double c10 = c010 * (1 - fx) + c110 * fx;
    double c01 = c001 * (1 - fx) + c101 * fx;
    double c11 = c011 * (1 - fx) + c111 * fx;
    
    double c0 = c00 * (1 - fy) + c10 * fy;
    double c1 = c01 * (1 - fy) + c11 * fy;
    
    return c0 * (1 - fz) + c1 * fz;
}

// ============================================================================
// TEST CASES
// ============================================================================

TEST_CASE("PoissonSolver: Point charge in vacuum", "[poisson][analytical][point-charge]") {
    // Test analytical solution: φ(r) = Q/(4πε₀r)
    
    SECTION("64³ grid, single electron at center") {
        const int N = 64;
        const double cell_size = 1e-4;  // 100 μm cells
        const double domain_size = N * cell_size;  // 6.4 mm domain
        
        Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
        grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
        
        PoissonSolver solver(grid);
        
        // Create charge density with single electron at center
        std::vector<double> rho(grid.size(), 0.0);
        int center_idx = grid.index(N/2, N/2, N/2);
        double Q = 1.6e-19;  // Elementary charge
        rho[center_idx] = Q / (cell_size * cell_size * cell_size);  // C/m³
        
        solver.setSourceTerm(rho);
        solver.solve(EPSILON_0, 1e-6, 5000);
        
        // Test points at various radii
        std::vector<double> test_radii = {0.0005, 0.001, 0.0015, 0.002};  // 0.5-2.0 mm
        
        for (double r : test_radii) {
            Vec3 test_pos = {r, 0, 0};  // Test along x-axis
            
            double phi_analytical = analytical_point_charge_potential(r, Q);
            double phi_computed = interpolate_potential(grid, test_pos);
            
            INFO("Testing at radius r = " << r*1e3 << " mm");
            INFO("Analytical potential: " << phi_analytical << " V");
            INFO("Computed potential: " << phi_computed << " V");
            
            // Allow 15% error due to grid discretization
            REQUIRE(phi_computed == Approx(phi_analytical).epsilon(0.15));
        }
    }
    
    SECTION("Convergence with grid refinement") {
        // Test that error decreases with finer grid (O(h²) expected)
        const double Q = 1.6e-19;
        const double test_radius = 0.001;  // 1 mm
        const double phi_analytical = analytical_point_charge_potential(test_radius, Q);
        
        std::vector<int> grid_sizes = {16, 32, 64};
        std::vector<double> errors;
        
        for (int N : grid_sizes) {
            const double cell_size = 0.01 / N;  // Keep domain size constant
            const double domain_size = 0.01;
            
            Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
            grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
            
            PoissonSolver solver(grid);
            
            // Single charge at center
            std::vector<double> rho(grid.size(), 0.0);
            int center_idx = grid.index(N/2, N/2, N/2);
            rho[center_idx] = Q / (cell_size * cell_size * cell_size);
            
            solver.setSourceTerm(rho);
            solver.solve(EPSILON_0, 1e-6, 5000);
            
            // Compute error at test point
            Vec3 test_pos = {test_radius, 0, 0};
            double phi_computed = interpolate_potential(grid, test_pos);
            double error = std::abs(phi_computed - phi_analytical);
            
            errors.push_back(error);
            
            INFO("Grid size: " << N << "³, Error: " << error << " V");
        }
        
        // Check error reduction with refinement
        // For 2x refinement, error should reduce by ~4x (O(h²))
        for (size_t i = 1; i < errors.size(); ++i) {
            double reduction_ratio = errors[i-1] / errors[i];
            INFO("Error reduction ratio (2x refinement): " << reduction_ratio);
            
            // Should be at least 2x reduction (conservative test)
            REQUIRE(reduction_ratio > 2.0);
        }
    }
}

TEST_CASE("PoissonSolver: Uniform charged sphere", "[poisson][analytical][sphere]") {
    // Test with uniform charge distribution (sphere)
    
    const int N = 64;
    const double cell_size = 1e-4;  // 100 μm cells
    const double domain_size = N * cell_size;
    
    Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
    grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    PoissonSolver solver(grid);
    
    // Create uniform sphere of charge
    const double sphere_radius = 0.001;  // 1 mm radius
    const double total_charge = 1e-15;   // 1 fC
    const double sphere_volume = (4.0/3.0) * M_PI * std::pow(sphere_radius, 3);
    const double rho_uniform = total_charge / sphere_volume;  // C/m³
    
    std::vector<double> rho(grid.size(), 0.0);
    
    // Fill sphere with uniform charge density
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            for (int k = 0; k < N; ++k) {
                Vec3 pos = {
                    grid.origin_m.x + i * cell_size + cell_size/2,
                    grid.origin_m.y + j * cell_size + cell_size/2,
                    grid.origin_m.z + k * cell_size + cell_size/2
                };
                
                double r = std::sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z);
                
                if (r < sphere_radius) {
                    int idx = grid.index(i, j, k);
                    rho[idx] = rho_uniform;
                }
            }
        }
    }
    
    solver.setSourceTerm(rho);
    solver.solve(EPSILON_0, 1e-5, 8000);  // Tighter tolerance for this test
    
    SECTION("Inside sphere") {
        // Test points inside sphere
        std::vector<double> test_radii = {0.0003, 0.0005, 0.0008};  // 0.3-0.8 mm
        
        for (double r : test_radii) {
            Vec3 test_pos = {r, 0, 0};
            
            double phi_analytical = analytical_sphere_potential(r, sphere_radius, total_charge);
            double phi_computed = interpolate_potential(grid, test_pos);
            
            INFO("Inside sphere at r = " << r*1e3 << " mm");
            INFO("Analytical: " << phi_analytical << " V, Computed: " << phi_computed << " V");
            
            // 20% tolerance inside sphere (harder to resolve)
            REQUIRE(phi_computed == Approx(phi_analytical).epsilon(0.20));
        }
    }
    
    SECTION("Outside sphere") {
        // Test points outside sphere (should look like point charge)
        std::vector<double> test_radii = {0.0015, 0.002, 0.0025};  // 1.5-2.5 mm
        
        for (double r : test_radii) {
            Vec3 test_pos = {r, 0, 0};
            
            double phi_analytical = analytical_sphere_potential(r, sphere_radius, total_charge);
            double phi_computed = interpolate_potential(grid, test_pos);
            
            INFO("Outside sphere at r = " << r*1e3 << " mm");
            INFO("Analytical: " << phi_analytical << " V, Computed: " << phi_computed << " V");
            
            // 15% tolerance outside sphere
            REQUIRE(phi_computed == Approx(phi_analytical).epsilon(0.15));
        }
    }
}

TEST_CASE("PoissonSolver: Charge conservation", "[poisson][conservation]") {
    // Verify that charge deposition conserves total charge
    
    const int N = 32;
    const double cell_size = 2e-4;  // 200 μm cells
    const double domain_size = N * cell_size;
    
    Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
    grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    // Create ion cloud
    std::vector<IonState> ions;
    const int num_ions = 100;
    const double charge_per_ion = 1.6e-19;
    
    for (int i = 0; i < num_ions; ++i) {
        IonState ion;
        ion.active = true;
        ion.born = true;
        ion.ion_charge_C = charge_per_ion;
        
        // Random position within central region
        ion.pos.x = (static_cast<double>(rand()) / RAND_MAX - 0.5) * 0.002;  // ±1 mm
        ion.pos.y = (static_cast<double>(rand()) / RAND_MAX - 0.5) * 0.002;
        ion.pos.z = (static_cast<double>(rand()) / RAND_MAX - 0.5) * 0.002;
        
        ions.push_back(ion);
    }
    
    // Deposit charges
    auto rho = deposit_charge(ions, grid);
    
    // Compute total charge on grid
    double Q_grid = 0.0;
    for (double rho_cell : rho) {
        Q_grid += rho_cell * (cell_size * cell_size * cell_size);
    }
    
    // Compute total ion charge
    double Q_ions = num_ions * charge_per_ion;
    
    INFO("Total ion charge: " << Q_ions << " C");
    INFO("Total grid charge: " << Q_grid << " C");
    INFO("Relative error: " << std::abs(Q_grid - Q_ions) / Q_ions);
    
    // Charge should be conserved to within 1%
    REQUIRE(Q_grid == Approx(Q_ions).epsilon(0.01));
}

TEST_CASE("PoissonSolver: Solver method comparison", "[poisson][solvers][performance]") {
    // All solver methods should give same result (within tolerance)
    
    const int N = 32;
    const double cell_size = 2e-4;
    const double domain_size = N * cell_size;
    
    // Create identical grids for each solver
    Grid3D grid_default(N, N, N, cell_size, cell_size, cell_size);
    Grid3D grid_hp(N, N, N, cell_size, cell_size, cell_size);
    
    grid_default.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    grid_hp.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    // Create test charge distribution
    std::vector<double> rho(grid_default.size(), 0.0);
    int center_idx = grid_default.index(N/2, N/2, N/2);
    rho[center_idx] = 1e-15 / (cell_size * cell_size * cell_size);
    
    // Solve with default method (Gauss-Seidel)
    PoissonSolver solver_default(grid_default);
    solver_default.setPerformanceMode(false);
    solver_default.setSourceTerm(rho);
    solver_default.solve(EPSILON_0, 1e-6, 5000);
    
    // Solve with high-performance method (Red-Black SOR)
    PoissonSolver solver_hp(grid_hp);
    solver_hp.setPerformanceMode(true);
    solver_hp.setSourceTerm(rho);
    solver_hp.solve(EPSILON_0, 1e-6, 5000);
    
    // Compare solutions
    double max_diff = 0.0;
    double max_val = 0.0;
    
    for (size_t i = 0; i < grid_default.size(); ++i) {
        double diff = std::abs(grid_hp.phi[i] - grid_default.phi[i]);
        max_diff = std::max(max_diff, diff);
        max_val = std::max(max_val, std::abs(grid_default.phi[i]));
    }
    
    double relative_diff = max_diff / max_val;
    
    INFO("Max absolute difference: " << max_diff << " V");
    INFO("Max potential value: " << max_val << " V");
    INFO("Relative difference: " << relative_diff);
    
    // Solvers should agree within 1%
    REQUIRE(relative_diff < 0.01);
}

TEST_CASE("PoissonSolver: Field computation from potential", "[poisson][field][gradient]") {
    // Test that E = -∇φ is computed correctly
    
    const int N = 32;
    const double cell_size = 2e-4;
    const double domain_size = N * cell_size;
    
    Grid3D grid(N, N, N, cell_size, cell_size, cell_size);
    grid.origin_m = {-domain_size/2, -domain_size/2, -domain_size/2};
    
    PoissonSolver solver(grid);
    
    // Create point charge at center
    std::vector<double> rho(grid.size(), 0.0);
    int center_idx = grid.index(N/2, N/2, N/2);
    double Q = 1.6e-19;
    rho[center_idx] = Q / (cell_size * cell_size * cell_size);
    
    solver.setSourceTerm(rho);
    solver.solve(EPSILON_0, 1e-6, 5000);
    solver.computeElectricField();
    
    // For point charge: E(r) = Q/(4πε₀r²) * r̂
    // Test at a few points
    
    SECTION("Field magnitude scales as 1/r²") {
        // Test closer to center (within resolved region)
        Vec3 pos1 = {0.0004, 0, 0};  // 0.4 mm from center (2 cells away)
        Vec3 pos2 = {0.0008, 0, 0};  // 0.8 mm from center (4 cells away)
        
        // Get grid indices (converting from world space to grid space)
        int i1 = static_cast<int>((pos1.x - grid.origin_m.x) / cell_size);
        int j1 = N/2;
        int k1 = N/2;
        
        int i2 = static_cast<int>((pos2.x - grid.origin_m.x) / cell_size);
        int j2 = N/2;
        int k2 = N/2;
        
        INFO("Grid indices: i1=" << i1 << " i2=" << i2 << " (N/2=" << N/2 << ")");
        
        // Check potential first (debug)
        double phi1 = grid.phi[grid.index(i1, j1, k1)];
        double phi2 = grid.phi[grid.index(i2, j2, k2)];
        double phi_center = grid.phi[grid.index(N/2, N/2, N/2)];
        
        INFO("Potentials: center=" << phi_center << " V, r=1mm=" << phi1 << " V, r=2mm=" << phi2 << " V");
        
        Vec3 E1 = grid.E[grid.index(i1, j1, k1)];
        Vec3 E2 = grid.E[grid.index(i2, j2, k2)];
        
        INFO("E1 components: Ex=" << E1.x << " Ey=" << E1.y << " Ez=" << E1.z);
        INFO("E2 components: Ex=" << E2.x << " Ey=" << E2.y << " Ez=" << E2.z);
        
        double E1_mag = std::sqrt(E1.x*E1.x + E1.y*E1.y + E1.z*E1.z);
        double E2_mag = std::sqrt(E2.x*E2.x + E2.y*E2.y + E2.z*E2.z);
        
        INFO("E at r=1mm: " << E1_mag << " V/m");
        INFO("E at r=2mm: " << E2_mag << " V/m");
        
        // Both fields should be non-zero
        REQUIRE(E1_mag > 0.0);
        REQUIRE(E2_mag > 0.0);
        
        // E1/E2 should be ≈ (r2/r1)² = 4
        double ratio = E1_mag / E2_mag;
        INFO("Ratio E1/E2: " << ratio << " (expected ~4)");
        
        // Allow 30% tolerance (coarse grid)
        REQUIRE(ratio == Approx(4.0).epsilon(0.30));
    }
}

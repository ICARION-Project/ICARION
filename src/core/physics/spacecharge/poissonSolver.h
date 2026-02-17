// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once
#include "core/types/Grid3D.h"
#include "fieldsolver/utils/grid/BoundaryConditions.h"
#include "utils/constants.h"
#include <vector>
#include <string>

/**
 * @class PoissonSolver
 * @brief Solves the Poisson or Laplace equation on a uniform 3D grid.
 *
 * Equation:
 *     ∇²φ = -ρ/ε₀       (Poisson)
 *     ∇²φ = 0            (Laplace)
 *
 * Current CPU implementation supports Dirichlet/Neumann boundary conditions via
 * explicit masks and uses Gauss-Seidel or Red-Black SOR iterations. Periodic/GPU/FFT
 * paths are stubbed for future work; callers should expect the GS/SOR pipeline.
 * Used by the CPU space-charge solver; geometry awareness is provided externally by
 * building grids from IDomainGeometry bounding boxes and applying Dirichlet masks.
 */
class PoissonSolver {
public:
    /**
     * @brief Construct Poisson solver for given 3D grid
     * @param grid Reference to Grid3D containing potential and field data
     */
    explicit PoissonSolver(Grid3D& grid);

    // --- Setup ---
    
    /**
     * @brief Set boundary conditions for solver
     * @param bc Boundary condition specification (Dirichlet, Neumann, periodic)
     */
    void setBoundaryConditions(const BoundaryConditions& bc);
    
    /**
     * @brief Set charge density source term
     * @param rho Charge density [C/m³] at each grid point (size = Nx×Ny×Nz)
     * 
     * For space-charge calculations, rho is computed from ion positions.
     * For static problems, rho = 0 (Laplace equation).
     */
    void setSourceTerm(const std::vector<double>& rho);

    // --- Solve Poisson or Laplace equation ---
    
    /**
     * @brief Solve Poisson equation: ∇²φ = -ρ/ε₀
     * @param eps0 Permittivity of free space [F/m] (default: EPSILON_0)
     * @param tol Convergence tolerance (relative residual)
     * @param max_iter Maximum number of iterations
     * 
     * Uses Gauss-Seidel or Red-Black SOR iteration depending on performance mode.
     * Updates grid.phi with computed potential values.
     */
    void solve(double eps0 = EPSILON_0, double tol = 1e-6, int max_iter = 5000);
    
    /**
     * @brief Adaptive solver with tolerance ramping for dynamic simulations
     * @param eps0 Permittivity of free space [F/m]
     * @param final_tol Final convergence tolerance
     * @param max_iter Maximum iterations
     * @param time_step Current simulation timestep (for adaptive tolerance)
     * 
     * Starts with looser tolerance for early timesteps, gradually tightening
     * to final_tol. Reduces computational cost for transient problems.
     */
    void solveAdaptive(double eps0 = EPSILON_0, double final_tol = 1e-6, int max_iter = 5000, int time_step = 0);

    // --- Performance settings ---
    
    /** @brief Enable high-performance mode (Red-Black SOR, parallelization) */
    void setPerformanceMode(bool high_performance) { m_high_performance = high_performance; }
    
    /** @brief Set number of OpenMP threads for parallel solve */
    void setParallelThreads(int threads) { m_threads = threads; }

    /**
     * @brief Set Dirichlet boundary mask
     * @param mask Binary mask (0=free node, 1=fixed potential). Size = Nx×Ny×Nz
     * 
     * Nodes with mask[i]==1 will not be updated during iteration.
     * Used for electrode surfaces and boundary conditions.
     */
    void setDirichletMask(const std::vector<char>& mask) { m_dirichlet_mask = mask; }
    
    /**
     * @brief Set fixed potential values for Dirichlet nodes
     * @param values Potential [V] at fixed nodes. Size = Nx×Ny×Nz
     * 
     * Only used at positions where dirichlet_mask==1.
     */
    void setDirichletValues(const std::vector<double>& values) { m_dirichlet_values = values; }

    // --- Compute E field from potential ---
    
    /**
     * @brief Compute electric field from potential: E = -∇φ
     * 
     * Uses finite differences to compute Ex, Ey, Ez from phi.
     * Updates grid.Ex, grid.Ey, grid.Ez arrays.
     * Should be called after solve() to obtain field values.
     */
    void computeElectricField();

    // --- Utility ---
    
    /** @brief Reset potential and field to zero */
    void reset();
    
    /** @brief Get final residual norm from last solve */
    double residualNorm() const { return m_lastResidual; }

private:
    Grid3D& m_grid;
    std::vector<double> m_rho;
    BoundaryConditions m_bc;
    double m_lastResidual = 0.0;
    bool m_high_performance = false;
    int m_threads = 1;
    // mask for Dirichlet fixed potential nodes (0 = free, 1 = fixed)
    std::vector<char> m_dirichlet_mask;
    // Values for Dirichlet nodes (only used where mask==1)
    std::vector<double> m_dirichlet_values;

    void solveGaussSeidel(double eps0, double tol, int max_iter);
    void solveRedBlack(double eps0, double tol, int max_iter);
    void solveMultiGrid(double eps0, double tol, int max_iter);
    void solveConjugateGradient(double eps0, double tol, int max_iter);
    void solveFFT(double eps0); // for periodic BCs
};

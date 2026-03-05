// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once
#include "core/types/Grid3D.h"
#include "core/physics/spacecharge/poissonSolver.h"
#include "core/physics/spacecharge/depositCharge.h"
#include "core/io/fieldArrayLoader.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include <vector>

namespace ICARION::config {
class IDomainGeometry;
}

/**
 * @class SpaceChargeSolver
 * @brief Poisson-based solver for self-consistent space-charge fields (CPU)
 *
 * CPU-only grid solver for box-shaped domains. Suitable for rough estimates; not
 * geometry-aware for cylindrical/Orbitrap setups and not validated for accuracy
 * beyond simple boxes. Falls back to zero outside the grid.
 * 
 * Workflow:
 * 1. Deposit ion charges onto 3D grid (cloud-in-cell or NGP)
 * 2. Solve Poisson equation for potential: ∇²φ = -ρ/ε₀
 * 3. Compute electric field: E = -∇φ
 * 4. Interpolate field at ion positions
 * 
 * Used by ICARION_Simulator for self-consistent trajectory integration.
 */
class SpaceChargeSolver {
public:
    /**
     * @brief Construct space-charge solver with specified grid resolution (box domain)
     * @param Nx Number of grid points in x direction
     * @param Ny Number of grid points in y direction
     * @param Nz Number of grid points in z direction
     * @param dx Grid spacing in x [m]
     * @param dy Grid spacing in y [m]
     * @param dz Grid spacing in z [m]
     * @param origin Grid origin position [m]
     */
    SpaceChargeSolver(int Nx, int Ny, int Nz,
                      double dx, double dy, double dz,
                      const Vec3& origin);

    /**
     * @brief Update electric field based on current ion positions
     * @param ions Vector of all ion states
     * 
     * Performs complete update cycle:
     * 1. Clear charge density grid
     * 2. Deposit ion charges (CIC or NGP)
     * 3. Solve Poisson equation
     * 4. Compute field from potential
     *
     * Box-grid only; geometry boundaries are not imposed beyond the grid extents.
     */
    void update(const std::vector<IonState>& ions);

    /**
     * @brief Update electric field using SoA ensemble
     */
    void update(const ICARION::core::IonEnsemble& ions);

    /**
     * @brief Interpolate space-charge field at given position
     * @param pos Position [m] in simulation domain
     * @return Electric field [V/m] from space charge
     * 
     * Uses trilinear interpolation on field grid.
     * Returns zero field if position is outside grid bounds.
     */
    Vec3 fieldAt(const Vec3& pos) const;

    // --- Performance settings ---
    
    /** @brief Enable high-performance mode (Red-Black SOR, parallelization) */
    void setPerformanceMode(bool high_performance_mode) { m_high_performance = high_performance_mode; }
    
    /** @brief Update field only every N timesteps (0 = every step) */
    void setUpdateFrequency(int every_n_steps) { m_update_frequency = every_n_steps; }
    
    /** @brief Cache field values to reduce interpolation overhead */
    void enableFieldCaching(bool enable) { m_cache_fields = enable; }
    
    /** @brief Set ion movement threshold for adaptive updates */
    void setAdaptiveUpdateThreshold(double ion_movement_threshold) { m_movement_threshold = ion_movement_threshold; }
    
    /**
     * @brief Check if field update is needed based on ion movement
     * @param ions Current ion ensemble
     * @return true if ions have moved sufficiently to warrant update
     * 
     * Compares current positions with cached positions from last update.
     * Useful for adaptive timestep schemes and sparse update strategies.
     */
    bool needsUpdate(const std::vector<IonState>& ions) const;
    bool needsUpdate(const ICARION::core::IonEnsemble& ions) const;

    // --- Expose PoissonSolver wrappers for external use ---
    
    /** @brief Set Dirichlet boundary mask (for electrode boundaries) */
    void setDirichletMask(const std::vector<char>& mask) { m_solver.setDirichletMask(mask); }
    
    /** @brief Set fixed potential values at Dirichlet boundaries */
    void setDirichletValues(const std::vector<double>& values) { m_solver.setDirichletValues(values); }
    
    /** @brief Set custom charge density (advanced use) */
    void setSourceTerm(const std::vector<double>& rho) { m_solver.setSourceTerm(rho); }
    
    /** @brief Manually trigger Poisson solve (for field building tool) */
    void solvePoisson(double eps0 = EPSILON_0, double tol = 1e-6, int max_iter = 5000) { m_solver.solve(eps0, tol, max_iter); }
    
    /** @brief Manually compute field from current potential */
    void computeField() { m_solver.computeElectricField(); }

    /** @brief Use geometry mask to reject out-of-domain ions during deposition */
    void setGeometryMask(const ICARION::config::IDomainGeometry* geometry) { m_geometry_mask = geometry; }

    // --- Preset configurations ---
    
    /** @brief Optimize for long simulations with frequent field updates */
    void configureForManyTimesteps(int typical_ion_count = 500);
    
    /** @brief Optimize for large ion ensembles (use coarser grid, faster methods) */
    void configureForManyIons(int typical_ion_count = 100000);
    
    /** @brief Access underlying grid for debugging/visualization */
    const Grid3D& grid() const { return m_grid; }

private:
    Grid3D m_grid;
    PoissonSolver m_solver;
    
    // Performance optimizations for many timesteps
    bool m_high_performance = false;
    int m_update_frequency = 10;  // Default: update every 10 steps for small ion counts
    bool m_cache_fields = true;
    double m_movement_threshold = 1e-5; // meters - update if ions move more than this
    
    mutable FieldArray m_cached_field;
    mutable bool m_field_cache_valid = false;
    size_t m_last_ion_count = 0;
    int m_step_counter = 0;
    std::vector<Vec3> m_last_ion_positions;
    const ICARION::config::IDomainGeometry* m_geometry_mask = nullptr;
};

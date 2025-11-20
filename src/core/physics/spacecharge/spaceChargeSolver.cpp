/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        spaceChargeSolver.cpp
 *   @brief       Fast Poisson-based solver for self-consistent space-charge fields.
 *
 *   @details
 *   Implements a SpaceChargeSolver class that uses a Poisson solver on a 3D grid to
 *  compute the electric field due to space charge from ion distributions. 
 *
 *   @date        <2025-10-18>
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#include "spaceChargeSolver.h"
#include <iostream>
#include "core/io/logger.h"

/* Constructor */
SpaceChargeSolver::SpaceChargeSolver(int Nx, int Ny, int Nz,
                                     double dx, double dy, double dz,
                                     const Vec3& origin)
    : m_grid(Nx, Ny, Nz, dx, dy, dz),
      m_solver(m_grid)
{
    m_grid.origin_m = origin;
    
    // Optimize defaults for few ions, many timesteps
    size_t grid_size = Nx * Ny * Nz;
    if (grid_size < 100000) {  // Small to medium grids
        m_update_frequency = 20;  // Update every 20 timesteps
        m_movement_threshold = 5e-6;  // 5 μm movement threshold
        ICARION::io::debug_log(std::string("[SpaceChargeSolver] Optimized for many timesteps mode (update every ") + std::to_string(m_update_frequency) + " steps)");
    }
}

/** 
* @brief Update electric field based on current ion positions
* @details
*  Computes the charge density from the current ion positions and updates the
*  electric field accordingly.
* @param ions Vector of current ions in the simulation
*/
bool SpaceChargeSolver::needsUpdate(const std::vector<IonState>& ions) const {
    // For first update or when ion count changes
    if (m_last_ion_positions.empty() || m_last_ion_positions.size() != ions.size()) {
        return true;
    }
    
    // Check if ions have moved significantly
    double max_movement = 0.0;
    for (size_t i = 0; i < ions.size() && i < m_last_ion_positions.size(); ++i) {
        if (!ions[i].active || !ions[i].born) continue;
        
        Vec3 displacement = ions[i].pos - m_last_ion_positions[i];
        double movement = norm(displacement);
        max_movement = std::max(max_movement, movement);
    }
    
    return max_movement > m_movement_threshold;
}

void SpaceChargeSolver::update(const std::vector<IonState>& ions)
{
    const size_t num_ions = ions.size();
    m_step_counter++;
    
    // Intelligent update logic for few ions, many timesteps
    bool should_update = false;
    
    if (num_ions <= 1000) {
        // Few ions: Use adaptive updating based on movement
        if (m_high_performance) {
            should_update = (m_step_counter % m_update_frequency == 0) || needsUpdate(ions);
        } else {
            should_update = needsUpdate(ions);
        }
    } else {
        // Many ions: Use fixed frequency
        should_update = (m_step_counter % m_update_frequency == 0);
    }
    
    if (!should_update) {
        return; // Use cached field from previous update
    }
    
    // 1. Charge deposition
    auto rho = deposit_charge(ions, m_grid);
    m_solver.setSourceTerm(rho);

    // 2. Optimized solve for few ions, many timesteps
    double tolerance;
    int max_iter;
    
    if (num_ions <= 1000) {
        // Small ion count - optimize for speed per timestep
        if (m_high_performance) {
            tolerance = 1e-4;     // Faster convergence
            max_iter = 200;       // Lower iteration limit
        } else {
            tolerance = 1e-6;     // Good balance
            max_iter = 400;       
        }
        // Silent operation for many timesteps
    } else if (num_ions <= 10000) {
        tolerance = 1e-6;
        max_iter = 500;
    } else if (num_ions <= 100000) {
        tolerance = 1e-4;
        max_iter = 300;
        if (num_ions != m_last_ion_count) {
            ICARION::io::debug_log(std::string("[SpaceCharge] Medium solve for ") + std::to_string(num_ions) + " ions");
        }
    } else {
        // Many ions: ultra-fast mode
        tolerance = 1e-3;
        max_iter = 100;
        if (num_ions != m_last_ion_count) {
            ICARION::io::debug_log(std::string("[SpaceCharge] Fast solve for ") + std::to_string(num_ions) + " ions");
        }
    }
    
    m_solver.solve(EPSILON_0, tolerance, max_iter);

    // 3. Compute E = -∇φ
    m_solver.computeElectricField();
    
    // 4. Invalidate field cache and store ion positions
    m_field_cache_valid = false;
    m_last_ion_count = num_ions;
    
    // Store current ion positions for next update check
    m_last_ion_positions.clear();
    m_last_ion_positions.reserve(ions.size());
    for (const auto& ion : ions) {
        m_last_ion_positions.push_back(ion.pos);
    }
}

/**
 * @brief Interpolate the field at a given ion position
 * @param pos Position at which to evaluate the field
 * @return Vec3 Electric field [V/m] at the given position
 */
Vec3 SpaceChargeSolver::fieldAt(const Vec3& pos) const
{
    // Use cached field array for performance (especially important for many ions)
    if (m_cache_fields && !m_field_cache_valid) {
        m_cached_field = grid_to_fieldarray(m_grid);
        m_field_cache_valid = true;
    }
    
    FieldArray* field_ptr = m_cache_fields ? &m_cached_field : nullptr;
    
    if (!m_cache_fields || !field_ptr->is_valid()) {
        // Fallback: direct conversion (slower but always works)
        FieldArray temp_field = grid_to_fieldarray(m_grid);
        return interpolate_field(temp_field, pos);
    }
    
    return interpolate_field(*field_ptr, pos);
}

void SpaceChargeSolver::configureForManyTimesteps(int typical_ion_count) {
    ICARION::io::debug_log(std::string("[SpaceChargeSolver] Configuring for many timesteps with ~") + std::to_string(typical_ion_count) + " ions");
    
    m_high_performance = true;
    m_cache_fields = true;
    
    if (typical_ion_count < 100) {
        m_update_frequency = 50;        // Very infrequent updates
        m_movement_threshold = 10e-6;   // 10 μm threshold
    } else if (typical_ion_count < 500) {
        m_update_frequency = 25;        // Update every 25 steps
        m_movement_threshold = 5e-6;    // 5 μm threshold  
    } else if (typical_ion_count < 2000) {
        m_update_frequency = 15;        // Update every 15 steps
        m_movement_threshold = 3e-6;    // 3 μm threshold
    } else {
        m_update_frequency = 10;        // Update every 10 steps
        m_movement_threshold = 2e-6;    // 2 μm threshold
    }
    
    ICARION::io::debug_log(std::string("  - Update frequency: every ") + std::to_string(m_update_frequency) + " timesteps");
    ICARION::io::debug_log(std::string("  - Movement threshold: ") + std::to_string(m_movement_threshold*1e6) + " μm");
}

void SpaceChargeSolver::configureForManyIons(int typical_ion_count) {
    ICARION::io::debug_log(std::string("[SpaceChargeSolver] Configuring for many ions (~") + std::to_string(typical_ion_count) + ")");
              
    m_high_performance = true;
    m_cache_fields = true;
    m_update_frequency = 1;         // Update every timestep
    m_movement_threshold = 0.0;     // Always update
    
    ICARION::io::debug_log("  - Update frequency: every timestep");
}

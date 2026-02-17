// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include "spaceChargeSolver.h"
#include <iostream>
#include "core/log/Logger.h"

// Adaptive update configuration constants
namespace {
    // Ion count thresholds for different optimization strategies
    constexpr int SMALL_ION_COUNT = 100;       ///< Very small ion ensembles
    constexpr int MEDIUM_ION_COUNT = 500;      ///< Small-medium ion ensembles
    constexpr int LARGE_ION_COUNT = 2000;      ///< Large ion ensembles
    constexpr int ADAPTIVE_THRESHOLD = 1000;   ///< Threshold for adaptive method selection
    
    // Update frequency settings (timesteps between updates)
    constexpr int UPDATE_FREQ_VERY_INFREQUENT = 50;  ///< For <100 ions, many timesteps
    constexpr int UPDATE_FREQ_INFREQUENT = 25;       ///< For <500 ions
    constexpr int UPDATE_FREQ_MODERATE = 15;         ///< For <2000 ions
    constexpr int UPDATE_FREQ_FREQUENT = 10;         ///< For >2000 ions
    constexpr int UPDATE_FREQ_DEFAULT_SMALL = 20;    ///< Default for small grids
    constexpr int UPDATE_FREQ_EVERY_STEP = 1;        ///< Update every timestep
    
    // Movement thresholds (meters) - trigger update if ions move more than this
    constexpr double MOVEMENT_THRESHOLD_LARGE = 10e-6;   ///< 10 μm for very small ion counts
    constexpr double MOVEMENT_THRESHOLD_MEDIUM = 5e-6;   ///< 5 μm for small ion counts
    constexpr double MOVEMENT_THRESHOLD_SMALL = 3e-6;    ///< 3 μm for medium ion counts
    constexpr double MOVEMENT_THRESHOLD_TINY = 2e-6;     ///< 2 μm for large ion counts
    constexpr double MOVEMENT_THRESHOLD_DEFAULT = 5e-6;  ///< Default 5 μm
    constexpr double MOVEMENT_THRESHOLD_ALWAYS = 0.0;    ///< Always update (no threshold)
    
    // Grid size threshold for performance optimization
    constexpr size_t SMALL_MEDIUM_GRID_THRESHOLD = 100000;  ///< 100k cells
}

/* Constructor */
SpaceChargeSolver::SpaceChargeSolver(int Nx, int Ny, int Nz,
                                     double dx, double dy, double dz,
                                     const Vec3& origin)
    : m_grid(Nx, Ny, Nz, dx, dy, dz),
      m_solver(m_grid)
{
    m_grid.origin_m = origin;
    
    // Validate grid parameters
    if (Nx < 8 || Ny < 8 || Nz < 8) {
        throw std::runtime_error("[SpaceChargeSolver] Grid too small (" + std::to_string(Nx) + "x" + 
                               std::to_string(Ny) + "x" + std::to_string(Nz) + "). Minimum: 8x8x8");
    }
    
    if (dx <= 0.0 || dy <= 0.0 || dz <= 0.0) {
        throw std::runtime_error("[SpaceChargeSolver] Invalid grid spacing: dx=" + std::to_string(dx) +
                               ", dy=" + std::to_string(dy) + ", dz=" + std::to_string(dz));
    }
    
    // Warning for very coarse grids (accuracy may be limited)
    if (Nx < 16 || Ny < 16 || Nz < 16) {
        ICARION::log::debug_log(std::string("[SpaceChargeSolver] WARNING: Low resolution grid (") + 
                               std::to_string(Nx) + "x" + std::to_string(Ny) + "x" + std::to_string(Nz) + 
                               "). Accuracy may be limited. Consider 32x32x32 or higher.");
    }
    
    // Optimize defaults for few ions, many timesteps
    size_t grid_size = Nx * Ny * Nz;
    if (grid_size < SMALL_MEDIUM_GRID_THRESHOLD) {
        m_update_frequency = UPDATE_FREQ_DEFAULT_SMALL;
        m_movement_threshold = MOVEMENT_THRESHOLD_DEFAULT;
        ICARION::log::debug_log(std::string("[SpaceChargeSolver] Optimized for many timesteps mode (update every ") + std::to_string(m_update_frequency) + " steps)");
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

bool SpaceChargeSolver::needsUpdate(const ICARION::core::IonEnsemble& ions) const {
    if (m_last_ion_positions.empty() || m_last_ion_positions.size() != ions.size()) {
        return true;
    }

    double max_movement = 0.0;
    const auto* active = ions.active_data();
    const auto* born = ions.born_data();
    const auto* pos_x = ions.pos_x_data();
    const auto* pos_y = ions.pos_y_data();
    const auto* pos_z = ions.pos_z_data();

    for (size_t i = 0; i < ions.size() && i < m_last_ion_positions.size(); ++i) {
        if (active[i] == 0 || born[i] == 0) continue;
        Vec3 displacement{pos_x[i] - m_last_ion_positions[i].x,
                          pos_y[i] - m_last_ion_positions[i].y,
                          pos_z[i] - m_last_ion_positions[i].z};
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
    
    if (num_ions <= ADAPTIVE_THRESHOLD) {
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
    
    // 1. Charge deposition (box grid; no geometry masking)
    auto rho = deposit_charge(ions, m_grid, m_geometry_mask);
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
            ICARION::log::debug_log(std::string("[SpaceCharge] Medium solve for ") + std::to_string(num_ions) + " ions");
        }
    } else {
        // Many ions: ultra-fast mode
        tolerance = 1e-3;
        max_iter = 100;
        if (num_ions != m_last_ion_count) {
            ICARION::log::debug_log(std::string("[SpaceCharge] Fast solve for ") + std::to_string(num_ions) + " ions");
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

void SpaceChargeSolver::update(const ICARION::core::IonEnsemble& ions)
{
    const size_t num_ions = ions.size();
    m_step_counter++;

    bool should_update = false;

    if (num_ions <= static_cast<size_t>(ADAPTIVE_THRESHOLD)) {
        if (m_high_performance) {
            should_update = (m_step_counter % m_update_frequency == 0) || needsUpdate(ions);
        } else {
            should_update = needsUpdate(ions);
        }
    } else {
        should_update = (m_step_counter % m_update_frequency == 0);
    }

    if (!should_update) {
        return;
    }

    auto rho = deposit_charge(ions, m_grid, m_geometry_mask);
    m_solver.setSourceTerm(rho);

    double tolerance;
    int max_iter;

    if (num_ions <= 1000) {
        if (m_high_performance) {
            tolerance = 1e-4;
            max_iter = 200;
        } else {
            tolerance = 1e-6;
            max_iter = 400;
        }
    } else if (num_ions <= 10000) {
        tolerance = 1e-6;
        max_iter = 500;
    } else if (num_ions <= 100000) {
        tolerance = 1e-4;
        max_iter = 300;
        if (num_ions != m_last_ion_count) {
            ICARION::log::debug_log(std::string("[SpaceCharge] Medium solve for ") + std::to_string(num_ions) + " ions");
        }
    } else {
        tolerance = 1e-3;
        max_iter = 100;
        if (num_ions != m_last_ion_count) {
            ICARION::log::debug_log(std::string("[SpaceCharge] Fast solve for ") + std::to_string(num_ions) + " ions");
        }
    }

    m_solver.solve(EPSILON_0, tolerance, max_iter);
    m_solver.computeElectricField();

    m_field_cache_valid = false;
    m_last_ion_count = num_ions;

    m_last_ion_positions.clear();
    m_last_ion_positions.reserve(num_ions);
    const auto* pos_x = ions.pos_x_data();
    const auto* pos_y = ions.pos_y_data();
    const auto* pos_z = ions.pos_z_data();
    for (size_t i = 0; i < num_ions; ++i) {
        m_last_ion_positions.push_back(Vec3{pos_x[i], pos_y[i], pos_z[i]});
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
    ICARION::log::debug_log(std::string("[SpaceChargeSolver] Configuring for many timesteps with ~") + std::to_string(typical_ion_count) + " ions");
    
    m_high_performance = true;
    m_cache_fields = true;
    
    if (typical_ion_count < SMALL_ION_COUNT) {
        m_update_frequency = UPDATE_FREQ_VERY_INFREQUENT;
        m_movement_threshold = MOVEMENT_THRESHOLD_LARGE;
    } else if (typical_ion_count < MEDIUM_ION_COUNT) {
        m_update_frequency = UPDATE_FREQ_INFREQUENT;
        m_movement_threshold = MOVEMENT_THRESHOLD_MEDIUM;
    } else if (typical_ion_count < LARGE_ION_COUNT) {
        m_update_frequency = UPDATE_FREQ_MODERATE;
        m_movement_threshold = MOVEMENT_THRESHOLD_SMALL;
    } else {
        m_update_frequency = UPDATE_FREQ_FREQUENT;
        m_movement_threshold = MOVEMENT_THRESHOLD_TINY;
    }
    
    ICARION::log::debug_log(std::string("  - Update frequency: every ") + std::to_string(m_update_frequency) + " timesteps");
    ICARION::log::debug_log(std::string("  - Movement threshold: ") + std::to_string(m_movement_threshold*1e6) + " μm");
}

void SpaceChargeSolver::configureForManyIons(int typical_ion_count) {
    ICARION::log::debug_log(std::string("[SpaceChargeSolver] Configuring for many ions (~") + std::to_string(typical_ion_count) + ")");
              
    m_high_performance = true;
    m_cache_fields = true;
    m_update_frequency = UPDATE_FREQ_EVERY_STEP;
    m_movement_threshold = MOVEMENT_THRESHOLD_ALWAYS;
    
    ICARION::log::debug_log("  - Update frequency: every timestep");
}

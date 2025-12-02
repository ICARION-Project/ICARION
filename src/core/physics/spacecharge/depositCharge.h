// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once
#include "core/types/IonState.h"
#include "core/types/Grid3D.h"
#include "utils/constants.h"
#include <vector>

/**
 * @brief Deposit ion charges onto 3D grid for space-charge calculations (box grid)
 * 
 * @param ions Vector of ion states with positions and charges
 * @param grid 3D grid specification (dimensions, spacing, origin)
 * 
 * @return Charge density [C/m³] at each grid point (size = Nx×Ny×Nz)
 * 
 * @note Current implementation uses CIC on a box grid. No geometry masking or boundary\n+ *       conditions are applied beyond the grid extents.\n*** End Patch
 * 
 * Maps discrete ion positions to continuous charge density field using
 * particle-in-cell (PIC) methods with second-order accuracy.
 * 
 * **Implemented deposition schemes:**
 * - **CIC (Cloud-In-Cell)**: Trilinear interpolation over 8 surrounding nodes
 *   O(h²) convergence, smooth fields, standard for publication-quality results.
 *   Charge conserving, parallel-safe with OpenMP atomics.
 *   Recommended for all production simulations.
 * 
 * **Legacy methods (deprecated):**
 * - **NGP (Nearest Grid Point)**: Only used in ultra-high performance mode (N>1M ions)
 *   O(h) convergence, causes grid noise, not recommended for publication.
 * 
 * **Future enhancement (v1.1+):**
 * - TSC (Triangular Shaped Cloud): Higher-order interpolation (27 nodes)
 *   O(h³) convergence, best smoothness but 3x slower than CIC
 * 
 * **Grid Resolution Requirements:**
 * - Minimum: Grid spacing < ion cloud size (otherwise under-resolved)
 * - Recommended: Grid spacing ≤ Debye length / 2
 * - Optimal: Grid spacing ≈ smallest feature size / 5
 * 
 * Charge density normalized by grid cell volume (dx×dy×dz).
 * Used as source term for Poisson solver in space-charge calculations.
 * 
 * **Algorithm (CIC):**
 * 1. For each ion, find lower-left grid node (i0, j0, k0)
 * 2. Compute fractional position within cell: fx, fy, fz ∈ [0,1]
 * 3. Calculate trilinear weights for 8 surrounding nodes (weights sum to 1.0)
 * 4. Distribute charge to all 8 nodes proportionally (atomic operations for thread-safety)
 * 5. Normalize by cell volume to get density [C/m³]
 * 
 * **Trilinear weighting formula:**
 *   w[i][j][k] = (1-fx or fx) × (1-fy or fy) × (1-fz or fz)
 * 
 * Example: w[0][0][0] = (1-fx)(1-fy)(1-fz) → weight for lower-left-back node
 * 
 * @throws std::runtime_error if grid spacing is zero or negative
 * @throws std::runtime_error if grid size is insufficient for ion distribution
 */
std::vector<double> deposit_charge(const std::vector<IonState>& ions,
                                   const Grid3D& grid);

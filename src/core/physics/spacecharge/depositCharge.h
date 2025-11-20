// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        depositCharge.h
 *   @brief       Utility functions for charge deposition onto a 3D grid in ICARION.
 *
 * @details
 * Provides:
 * - Function to deposit ion charges onto a 3D grid.
 * - Support for different deposition methods (NGP, CIC, TSC).
 * - Charge density computation from ion positions and charges.
 *
 *   @date        2025-10-18
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#pragma once
#include "core/types/IonState.h"
#include "core/types/Grid3D.h"
#include "utils/constants.h"
#include <vector>

/**
 * @brief Deposit ion charges onto 3D grid for space-charge calculations
 * 
 * @param ions Vector of ion states with positions and charges
 * @param grid 3D grid specification (dimensions, spacing, origin)
 * 
 * @return Charge density [C/m³] at each grid point (size = Nx×Ny×Nz)
 * 
 * Maps discrete ion positions to continuous charge density field using
 * particle-in-cell (PIC) methods. Supports multiple deposition schemes:
 * 
 * - NGP (Nearest Grid Point): Assigns full charge to nearest grid node
 *   Fast but can cause numerical noise
 * 
 * - CIC (Cloud-In-Cell): Distributes charge over 8 surrounding nodes
 *   using trilinear weighting. Smoother, reduces grid artifacts.
 * 
 * - TSC (Triangular Shaped Cloud): Higher-order interpolation (27 nodes)
 *   Best smoothness but more expensive
 * 
 * Charge density normalized by grid cell volume (dx×dy×dz).
 * Used as source term for Poisson solver in space-charge calculations.
 * 
 * Algorithm:
 * 1. For each ion, find enclosing grid cell
 * 2. Compute interpolation weights based on fractional position
 * 3. Distribute charge to surrounding nodes proportionally
 * 4. Normalize by cell volume to get density [C/m³]
 */
std::vector<double> deposit_charge(const std::vector<IonState>& ions,
                                   const Grid3D& grid);

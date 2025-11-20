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
 *   @file        defineSpaceChargeField.cpp
 *   @brief       Calculates the electric fields due to space charge.
 *
 *   @details
 *   Provides the function to compute the Coulombic electric field due to space charge.
 *
 *
 *   @date        2025-10-17
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#pragma once
#include "core/types/IonState.h"
#include <vector>

// Compute Coulombic space-charge field acting on one ion
Vec3 SpaceChargeField(const IonState& ion,
                      const std::vector<IonState>& ions,
                      double eps0);

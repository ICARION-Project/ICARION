// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        SimulationResult.h
 *   @brief       Defines struct to store simulation results.
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 * 
 * =====================================================================
 */

#pragma once
#include <vector>
#include "core/types/IonState.h"

/**
 * @struct SimulationResult
 * @brief Stores the final states of all ions after a simulation.
 *
 * @details
 * Holds a vector of IonState objects, including positions, velocities,
 * and optionally arrival times or domain indices. Can be used to print
 * results, save to file, or analyze trajectories.
 */
struct SimulationResult {
    std::vector<IonState> ions;
};
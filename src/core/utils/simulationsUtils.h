// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file simulationsUtils.h
 * @brief Simulation initialization and output utilities
 */
#pragma once

#include <string>
#include <utility>
#include <vector>
#include "core/param/paramUtils.h"
#include "core/integrator/integrator.h"
#include "core/types/SimulationResult.h"
#include "core/io/speciesLoader.h"

namespace ICARION {
namespace utils {

// [REMOVED] init_ions() - replaced by FullConfig::generate_ions() in config system

/**
 * @brief Print instrument domain configuration summary
 * @param domains Vector of instrument domains
 * 
 * Outputs to console:
 * - Domain index and type (LQIT, SIFDT, IMS, etc.)
 * - Geometry dimensions (length, radius, bounds)
 * - Voltages (DC, RF, AC) and frequencies
 * - Environment (pressure, temperature, gas type)
 * - Collision model and field source
 */
void print_domain_summary(std::vector<InstrumentDomain>& domains);

/**
 * @brief Print simulation results summary
 * @param result Simulation result with final ion states
 * @param max_nr_ions Maximum number of ions to print (limits output)
 * 
 * Outputs:
 * - Number of active vs lost ions
 * - Final positions and velocities (for first max_nr_ions)
 * - Arrival times at detector (if applicable)
 * - Domain indices
 * - Species distribution
 */
void print_results(const SimulationResult& result, size_t max_nr_ions);

}  // namespace utils
}  // namespace ICARION

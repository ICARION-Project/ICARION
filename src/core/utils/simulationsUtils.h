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

/**
 * @brief Initialize ion ensemble from configuration
 * @param gParams Global parameters (ion count, species, positions)
 * @param speciesDB Species database with physical properties
 * @param domains Instrument domain configurations
 * @return Vector of initialized IonState objects
 * 
 * Initialization modes:
 * 1. From ion list file: Load positions, velocities, species from HDF5/CSV
 * 2. Random initialization: Generate ions with random positions in domain
 * 3. Single ion: Create one ion at specified position
 * 
 * Sets initial properties:
 * - Position and velocity from config or random
 * - Mass, charge, mobility, CCS from species database
 * - Domain index from position
 * - Birth time = 0 for initial ions
 */
std::vector<IonState> init_ions(GlobalParams& gParams,
                                const ICARION::io::SpeciesDatabase& speciesDB,
                                const std::vector<InstrumentDomain>& domains);

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

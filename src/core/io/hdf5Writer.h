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
 *   @file        hdf5Writer.h
 *   @brief       Writes simulation parameters to HDF5 output files.
 *
 *   @details
 *   Provides functions to serialize global simulation parameters and
 *   instrument domain configurations into HDF5 files. Supports scalar,
 *   string, and array data types. Each instrument domain (geometry,
 *   fields, environment) is saved in a hierarchical structure.
 *
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#pragma once

#include <string>
#include <vector>
#include "core/param/paramUtils.h"
#include "core/types/IonState.h"
#include "H5Cpp.h"

namespace ICARION {
namespace io {

// -----------------------------
// Helper functions
// -----------------------------
template <typename T>
const H5::PredType& getH5Type();
template <typename T>
void writeScalar(H5::Group& group, const std::string& name, const T& value);
void writeBool(H5::Group& group, const std::string& name, bool value);
void writeString(H5::Group& group, const std::string& name, const std::string& value);
template <typename T>
void writeArray(H5::Group& group, const std::string& name, const std::vector<T>& values);

// -----------------------------
// HDF Writer
// -----------------------------

/**
 * @brief Write simulation parameters to HDF5 file
 * @param file HDF5 file handle (must be open)
 * @param gParams Global simulation parameters
 * @param domains Vector of instrument domain configurations
 */
void write_params_to_HDF5(H5::H5File& file, const GlobalParams& gParams,
                          const std::vector<InstrumentDomain>& domains);

/**
 * @brief Append trajectory data to existing HDF5 file
 * @param filename Path to HDF5 file
 * @param times_buffer Time points for trajectory snapshots [s]
 * @param trajectory_buffer Ion states at each time point
 */
void append_to_HDF5(const std::string& filename,
                    const std::vector<double>& times_buffer,
                    const std::vector<std::vector<IonState>>& trajectory_buffer);

// -----------------------------
// Continue mode support
// -----------------------------

/**
 * @brief Write species metadata to /metadata/species group
 * @param file HDF5 file handle
 * @param ions Vector of ion states (to extract unique species)
 * 
 * Writes species-constant properties (mass, charge, mobility, CCS) once per species.
 * This avoids redundant storage of these values at every timestep.
 */
void write_species_metadata(H5::H5File& file, const std::vector<IonState>& ions);

/**
 * @brief Write simulation metadata (completion status, RNG seed, git hash)
 * @param file HDF5 file handle
 * @param complete True if simulation finished normally
 * @param active_ions Number of ions still active at end
 * @param final_time_s Final simulation time [s]
 * @param rng_seed Random number generator seed for reproducibility
 * @param git_hash Git commit hash of build
 * @param config_json Full configuration JSON string (optional)
 */
void write_simulation_metadata(H5::H5File& file, bool complete, int active_ions, double final_time_s, 
                               unsigned int rng_seed = 0, const std::string& git_hash = "",
                               const std::string& config_json = "");

/**
 * @brief Load final ion states from HDF5 checkpoint for continue mode
 * @param filename Path to HDF5 checkpoint file
 * @return Vector of ion states loaded from file
 */
std::vector<IonState> load_final_state_from_HDF5(const std::string& filename);

}  // namespace io
}  // namespace ICARION

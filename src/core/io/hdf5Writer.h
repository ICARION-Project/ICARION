// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file hdf5Writer.h
 * @brief Modern HDF5 writer using FullConfig
 * 
 * Clean implementation based on new configuration system.
 * No legacy GlobalParams support.
 * 
 * @see docs/HDF5_OUTPUT_STRUCTURE.md for file format specification
 */

#pragma once

#include "core/config/types/FullConfig.h"
#include "core/types/IonEnsemble.h"
#include <H5Cpp.h>
#include <vector>
#include <string>
#include <cstdint>

namespace ICARION::io {

/**
 * @brief Modern HDF5 writer for FullConfig-based simulations
 * 
 * Creates HDF5 files with the following structure:
 * - /metadata/          Configuration, reproducibility, system info
 * - /trajectory/        Time-series data (positions, velocities, etc.)
 * - /ions/              Per-ion metadata (initial conditions)
 * - /domains/           Domain configurations
 */
class HDF5Writer {
public:
    /**
     * @brief Create HDF5 file and write metadata
     * 
     * @param filename Output HDF5 file path
     * @param config Complete configuration
     * @param ions Initial ion states
     * @param git_hash Git commit hash
     * @param build_info Compiler/build information
     * 
     * Creates file structure and writes all metadata groups.
     * Trajectory data is written later via append_trajectory().
     */
    static void create_file(
        const std::string& filename,
        const config::FullConfig& config,
        const core::IonEnsemble& ions,
        const std::string& git_hash,
        const std::string& build_info
    );
    
    /**
     * @brief Append trajectory snapshot
     * 
     * @param filename HDF5 file to append to
     * @param time Current simulation time [s]
     * @param ions Current ion states
     * 
     * Appends one timestep to trajectory datasets.
     * 
     * @note For better performance, use append_trajectory_batch() when writing multiple timesteps.
     */
    static void append_trajectory(
        const std::string& filename,
        double time,
        const core::IonEnsemble& ions
    );
    
    /**
     * @brief Append multiple trajectory snapshots in batch
     * 
     * @param filename HDF5 file to append to
     * @param times Vector of simulation times [s]
     * @param trajectories Vector of ion state vectors (one per timestep)
     * 
     * Much more efficient than calling append_trajectory() in a loop.
     * Opens file once, extends datasets once, writes all data, then closes.
     * Provides ~100x speedup for buffered writes.
     */
    static void append_trajectory_batch(
        const std::string& filename,
        const std::vector<double>& times,
        const std::vector<core::IonEnsemble>& trajectories
    );

    /**
     * @brief Append multiple trajectory snapshots using flattened buffers (SoA)
     *
     * @param filename HDF5 file to append to
     * @param times Vector of simulation times [s] (size = n_steps)
     * @param n_ions Number of ions per step (constant across steps)
     * @param positions Flattened positions [n_steps * n_ions * 3]
     * @param velocities Flattened velocities [n_steps * n_ions * 3]
     * @param domain_indices Flattened domain indices [n_steps * n_ions]
     * @param species_indices Flattened species pool indices [n_steps * n_ions]
     *
     * Faster path for buffered SoA output without IonEnsemble copies.
     * Species indices correspond to the pool written in /ions/.
     */
    static void append_trajectory_batch_flat(
        const std::string& filename,
        const std::vector<double>& times,
        size_t n_ions,
        const std::vector<double>& positions,
        const std::vector<double>& velocities,
        const std::vector<int>& domain_indices,
        const std::vector<uint32_t>& species_indices,
        const std::vector<std::string>* species_pool = nullptr,
        const std::vector<double>& per_ion_times = {}
    );
    
    /**
     * @brief Update death_time_s dataset with final ion states
     * 
     * @param filename HDF5 file to update
     * @param final_ions Final ion states at end of simulation
     * 
     * Called before finalize() to ensure death_time_s reflects
     * boundary absorption times. Initial write in write_ion_metadata()
     * sets all to -1 (alive). This updates with actual death times.
     */
    static void update_death_times(
        const std::string& filename,
        const core::IonEnsemble& final_ensemble
    );
    
    /**
     * @brief Write simulation completion metadata
     * 
     * @param filename HDF5 file
     * @param success Simulation completed successfully?
     * @param final_time Final simulation time [s]
     * @param active_ions Number of active ions at end
     */
    static void finalize(
        const std::string& filename,
        bool success,
        double final_time,
        size_t active_ions
    );

private:
    // === Metadata writers ===
    
    /**
     * @brief Write /metadata/config/
     * 
     * Stores complete FullConfig as JSON string + key parameters as datasets
     */
    static void write_config_metadata(
        H5::H5File& file,
        const config::FullConfig& config
    );
    
    /**
     * @brief Write /metadata/reproducibility/
     * 
     * Git hash, RNG seed, compiler version, build flags
     */
    static void write_reproducibility_metadata(
        H5::H5File& file,
        const config::FullConfig& config,
        const std::string& git_hash,
        const std::string& build_info
    );
    
    /**
     * @brief Write /metadata/system/
     * 
     * Hostname, OS, CPU, GPU, memory, timestamp
     */
    static void write_system_metadata(H5::H5File& file);
    
    /**
     * @brief Write /metadata/species/
     * 
     * Tabular format: names, masses, charges, mobilities, CCS
     * Only writes species referenced by ions (filters large databases)
     */
    static void write_species_metadata(
        H5::H5File& file,
        const config::SpeciesDatabase& species_db,
        const core::IonEnsemble& ions
    );
    
    /**
     * @brief Write /metadata/reactions/
     * 
     * Tabular format: id, reactants, products, rate constants, types
     * Only writes reactions involving species present in ions (filters large databases)
     */
    static void write_reactions_metadata(
        H5::H5File& file,
        const config::ReactionDatabase& reaction_db,
        const core::IonEnsemble& ions,
        const config::SpeciesDatabase& species_db
    );
    
    // === Domain writers ===
    
    /**
     * @brief Write /domains/ hierarchy
     * 
     * One group per domain with geometry, environment, fields
     */
    static void write_domains(
        H5::H5File& file,
        const std::vector<config::DomainConfig>& domains
    );
    
    /**
     * @brief Write single domain group
     */
    static void write_domain(
        H5::Group& parent,
        const config::DomainConfig& domain,
        size_t index
    );
    
    /**
     * @brief Write waveform library to HDF5 group
     * 
     * Serializes waveforms for reproducibility.
     * Each waveform becomes a subgroup with type and parameters.
     */
    static void write_waveform_library(
        H5::Group& parent,
        const std::map<std::string, config::Waveform>& library
    );
    
    /**
     * @brief Write single waveform to HDF5 group
     */
    static void write_waveform(
        H5::Group& parent,
        const std::string& name,
        const config::Waveform& waveform
    );
    
    // === Ion metadata ===
    
    /**
     * @brief Write /ions/ initial conditions
     * 
     * Per-ion metadata: species, initial pos/vel, birth time
     */
    static void write_ion_metadata(
        H5::H5File& file,
        const core::IonEnsemble& ions
    );
    
    // === Helpers ===
    
    static void write_scalar(H5::Group& group, const std::string& name, double value);
    static void write_scalar(H5::Group& group, const std::string& name, int value);
    static void write_scalar(H5::Group& group, const std::string& name, unsigned int value);
    static void write_scalar(H5::Group& group, const std::string& name, bool value);
    static void write_string(H5::Group& group, const std::string& name, const std::string& value);
    static void write_blob(H5::Group& group, const std::string& name, const std::string& data);
    static void write_array(H5::Group& group, const std::string& name, const std::vector<double>& data);
    static void write_vec3(H5::Group& group, const std::string& name, const Vec3& vec);
};

} // namespace ICARION::io

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file OutputManager.h
 * @brief Manages simulation output: HDF5 trajectories + text logging
 * 
 * Part of Phase 5A refactoring: SimulationEngine
 * 
 * Responsibilities:
 * - Buffer trajectory data in RAM (reduce HDF5 I/O overhead)
 * - Periodic HDF5 writes (time-based or size-based)
 * - Metadata export (species, parameters, git hash, version)
 * - Text logging (progress, statistics, completion summary)
 * - Config snapshot written separately by main.cpp (not handled here)
 * 
 * Wraps both HDF5Writer and RunLogger for unified output API.
 */

#pragma once

#include "core/config/types/FullConfig.h"
#include "core/types/IonEnsemble.h"
#include "core/log/Logger.h"
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <mutex>
#include <cstdint>

namespace ICARION::io { class HDF5Writer; }

namespace ICARION {
namespace integrator {

/**
 * @class OutputManager
 * @brief Unified output manager for HDF5 + text logging
 * 
 * Handles both binary trajectory output (HDF5) and human-readable
 * text logs (progress, statistics, completion summary).
 * 
 * **Buffering Strategy:**
 * - Trajectory snapshots buffered in RAM (configurable size)
 * - Periodic flush to HDF5 (time-based or buffer-full)
 * - Reduces I/O overhead for long simulations
 * 
 * **Memory Usage Warning:**
 * Current implementation stores full IonState vectors in RAM.
 * For large ensembles (>100k ions), this can consume significant memory.
 * 
 * **Text Logging (optional):**
 * - Progress messages ("50% completed")
 * - Ion statistics (active/lost counts)
 * - Completion summary
 * - Can be disabled by passing empty log filename
 * 
 * **Thread Safety:**
 * - HDF5 writes are serialized (external to this class)
 * - Text logging protected by internal mutex (std::lock_guard)
 * - Safe for concurrent log_progress() calls from multiple threads
 */
class OutputManager {
public:
    /**
     * @brief Construct output manager
     * @param hdf5_filename Path to HDF5 trajectory file (required)
     * @param log_filename Path to text log file (empty = no text log)
     * @param write_interval_dt Time interval between HDF5 writes [s]
     * @param buffer_max Max timesteps in RAM before forced flush
     * 
     * Creates HDF5 file on initialization (via initialize()).
     * Text log file is created on first log_progress() call.
     */
    OutputManager(const std::string& hdf5_filename,
                  const std::string& log_filename = "",
                  double write_interval_dt = 0.001,  // Default: 1 ms
                  size_t buffer_max = 50);
    
    /**
     * @brief Destructor - ensures final flush
     */
    ~OutputManager();
    
    /**
     * @brief Initialize output system
     * @param config Simulation configuration (SSOT)
     * @param ions Initial ion ensemble (for species metadata)
     * 
     * Creates HDF5 file and writes metadata:
     * - Simulation parameters (timestep, duration, collision model)
     * - Domain configurations
     * - Species list (mass, charge, CCS)
     * - Git hash and version
     * 
     * Also writes text log header if text logging enabled.
     * 
     * Must be called before log_step().
     */
    void initialize(const config::FullConfig& config,
                    const core::IonEnsemble& ensemble);
    
    /**
     * @brief Log timestep snapshot (SoA variant - Phase 5)
     * @param t Current simulation time [s]
     * @param ensemble Current ion ensemble (SoA format)
     * 
     * Direct SoA→HDF5 writing without to_legacy() conversion.
     * Appends snapshot to buffer. Same trigger logic as log_step().
     */
    void log_step(double t, const core::IonEnsemble& ensemble);
    
    /**
     * @brief Log progress message (to text log only)
     * @param message Progress message (e.g., "50% completed (t = 0.5 ms)")
     * 
     * Does nothing if text logging disabled (log_filename empty).
     * Thread-safe (uses mutex for text log access).
     */
    void log_progress(const std::string& message);
    
    /**
     * @brief Check if HDF5 write is needed
     * @param t_current Current simulation time [s]
     * @return true if flush should be triggered
     * 
     * Returns true if:
     * - Buffer is full (times_buffer_.size() >= buffer_max_)
     * - Time interval exceeded (t_current >= next_write_time_)
     */
    bool should_write(double t_current) const;
    
    /**
     * @brief Set maximum buffer size in bytes to guard against OOM.
     */
    void set_buffer_byte_cap(size_t bytes) { buffer_byte_cap_ = bytes; }
    
private:
    /**
     * @brief Internal check before adding to buffer
     * @param t_current Current simulation time [s]
     * @return true if flush should be triggered before adding new element
     * 
     * Same logic as should_write(), but used internally in log_step()
     * to allow buffer to fill to buffer_max exactly.
     */
    bool should_write_before_add(double t_current) const;
    
public:
    /**
     * @brief Flush buffers to HDF5 file
     * 
     * Writes all buffered trajectory snapshots to HDF5.
     * Clears buffers after successful write.
     * Updates next_write_time_.
     * 
     * Thread-safe (HDF5 library handles file locking).
     */
    void flush();
    
    /**
     * @brief Finalize and close output files (SoA variant)
     * @param t_final Final simulation time [s]
     * @param final_ensemble Final ion ensemble states
     */
    void finalize(double t_final, const core::IonEnsemble& final_ensemble);
    
    /**
     * @brief Get HDF5 filename
     */
    const std::string& get_hdf5_filename() const { return hdf5_filename_; }
    
    /**
     * @brief Get text log filename (empty if disabled)
     */
    const std::string& get_log_filename() const { return log_filename_; }
    
    /**
     * @brief Check if text logging is enabled
     */
    bool has_text_log() const { return !log_filename_.empty(); }
    
    /**
     * @brief Get number of buffered snapshots
     */
    size_t buffer_size() const { return times_buffer_.size(); }
    
private:
    // HDF5 output
    std::string hdf5_filename_;
    bool trajectory_enabled_ = true;            ///< full mode=true, minimal mode=false
    double write_interval_dt_;
    double next_write_time_;
    size_t buffer_max_;
    
    // HDF5 buffers
    // Flattened buffers (avoid full IonEnsemble copies)
    size_t n_ions_ = 0;
    std::vector<double> times_buffer_;
    std::vector<double> per_ion_time_buffer_;   // [steps * n_ions]
    std::vector<double> positions_buffer_;   // contiguous [steps * n_ions * 3]
    std::vector<double> velocities_buffer_;  // contiguous [steps * n_ions * 3]
    std::vector<int> domain_buffer_;         // [steps * n_ions]
    std::vector<uint32_t> species_buffer_;   // [steps * n_ions] (pool indices)
    const std::vector<std::string>* species_pool_ = nullptr; // non-owning; set on init
    size_t buffer_byte_cap_ = 0;  ///< Optional byte cap; 0 disables
    
    // Text logging (optional)
    std::string log_filename_;
    std::ofstream text_log_file_;
    std::mutex text_log_mutex_;  ///< Thread-safety for text logging
    
    // State tracking
    bool initialized_ = false;
    bool finalized_ = false;     ///< Has finalize() been called?
    double last_time_ = 0.0;     ///< Last logged time (for incomplete metadata)
    size_t total_writes_ = 0;    ///< Total number of HDF5 writes
};

}  // namespace integrator
}  // namespace ICARION

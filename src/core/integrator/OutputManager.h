// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

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
 * 
 * Wraps both HDF5Writer v2 and RunLogger for unified output API.
 */

#pragma once

#include "core/config/types/FullConfig.h"
#include "core/types/IonState.h"
#include "core/log/Logger.h"
#include <string>
#include <vector>
#include <memory>
#include <fstream>
#include <mutex>

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
 * TODO(v1.1): Implement memory-efficient output modes:
 * - positions_only mode (skip velocity, skip inactive ions)
 * - reduced_precision mode (float instead of double)
 * - sparse_logging mode (log every N-th ion)
 * - streaming mode (extendible datasets, no buffering)
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
                    const std::vector<IonState>& ions);
    
    /**
     * @brief Log trajectory snapshot (buffers in RAM)
     * @param t Current simulation time [s]
     * @param ions Current ion states
     * 
     * Appends snapshot to buffer. Triggers flush if:
     * - Buffer is full (size >= buffer_max)
     * - Time interval exceeded (t >= next_write_time)
     */
    void log_step(double t, const std::vector<IonState>& ions);
    
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
     * @brief Finalize output (write remaining data, close files)
     * @param t_final Final simulation time [s]
     * @param final_ions Final ion states
     * 
     * Performs:
     * 1. Flush remaining HDF5 buffers
     * 2. Write completion metadata to HDF5 (success flag, active ions count)
     * 3. Write text log completion summary (if enabled)
     * 
     * Should be called at end of simulation (after last log_step()).
     */
    void finalize(double t_final, const std::vector<IonState>& final_ions);
    
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
    double write_interval_dt_;
    double next_write_time_;
    size_t buffer_max_;
    
    // HDF5 buffers
    std::vector<double> times_buffer_;
    std::vector<std::vector<IonState>> trajectory_buffer_;
    
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

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "OutputManager.h"
#include "core/io/hdf5Writer.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>

// Git hash from CMake (via compile definition)
#ifndef GIT_HASH
#define GIT_HASH "unknown"
#endif

// Build info
#ifdef NDEBUG
#define BUILD_TYPE "Release"
#else
#define BUILD_TYPE "Debug"
#endif

namespace ICARION {
namespace integrator {

OutputManager::OutputManager(
    const std::string& hdf5_filename,
    const std::string& log_filename,
    double write_interval_dt,
    size_t buffer_max
) : hdf5_filename_(hdf5_filename),
    log_filename_(log_filename),
    write_interval_dt_(write_interval_dt),
    next_write_time_(0.0),
    buffer_max_(buffer_max)
{
    if (hdf5_filename_.empty()) {
        throw std::invalid_argument("OutputManager: HDF5 filename cannot be empty");
    }
    
    // Reserve buffer space
    times_buffer_.reserve(buffer_max_);
    trajectory_buffer_.reserve(buffer_max_);
}

OutputManager::~OutputManager() {
    // Ensure final flush on destruction (if not finalized)
    if (!times_buffer_.empty()) {
        try {
            flush();
        } catch (const std::exception& e) {
            std::cerr << "Warning: OutputManager destructor failed to flush: " 
                      << e.what() << std::endl;
        }
    }
}

void OutputManager::initialize(
    const config::FullConfig& config,
    const std::vector<IonState>& ions
) {
    if (initialized_) {
        throw std::runtime_error("OutputManager: Already initialized");
    }
    
    // 1. Create HDF5 file with metadata using HDF5Writer
    try {
        io::HDF5Writer::create_file(
            hdf5_filename_,
            config,
            ions,
            GIT_HASH,
            BUILD_TYPE
        );
    } catch (const std::exception& e) {
        throw std::runtime_error("OutputManager: Failed to create HDF5 file: " + 
                                 std::string(e.what()));
    }
    
    // 2. Initialize text logger (if enabled)
    if (!log_filename_.empty()) {
        text_logger_ = std::make_unique<io::RunLogger>(log_filename_);
        text_logger_->writeHeader();
        
        // Note: writeGlobalParams() and writeInstrumentDomains() removed (legacy)
        // Text logger now only used for progress messages via log()
    }
    
    initialized_ = true;
    next_write_time_ = write_interval_dt_;
}

void OutputManager::log_step(double t, const std::vector<IonState>& ions) {
    if (!initialized_) {
        throw std::runtime_error("OutputManager: Not initialized (call initialize() first)");
    }
    
    // Check if flush needed BEFORE adding (allows buffer to fill to buffer_max)
    if (should_write_before_add(t)) {
        flush();
    }
    
    // Buffer snapshot
    times_buffer_.push_back(t);
    trajectory_buffer_.push_back(ions);
}

void OutputManager::log_progress(const std::string& message) {
    if (text_logger_) {
        text_logger_->log(message);
    }
}

bool OutputManager::should_write(double t_current) const {
    // Check if write is needed (used by external callers)
    return (times_buffer_.size() >= buffer_max_) || (t_current >= next_write_time_);
}

bool OutputManager::should_write_before_add(double t_current) const {
    // Check BEFORE adding to buffer (allows buffer to reach buffer_max exactly)
    return (times_buffer_.size() >= buffer_max_) || (t_current >= next_write_time_);
}

void OutputManager::flush() {
    if (times_buffer_.empty()) {
        return;  // Nothing to flush
    }
    
    try {
        // Write buffered snapshots using HDF5Writer static methods
        for (size_t i = 0; i < times_buffer_.size(); ++i) {
            io::HDF5Writer::append_trajectory(
                hdf5_filename_,
                times_buffer_[i],
                trajectory_buffer_[i]
            );
        }
        
        // Clear buffers
        times_buffer_.clear();
        trajectory_buffer_.clear();
        
        // Update next write time
        next_write_time_ += write_interval_dt_;
        total_writes_++;
        
    } catch (const std::exception& e) {
        throw std::runtime_error("OutputManager: HDF5 flush failed: " + 
                                 std::string(e.what()));
    }
}

void OutputManager::finalize(double t_final, const std::vector<IonState>& final_ions) {
    if (!initialized_) {
        return;  // Nothing to finalize
    }
    
    // 1. Ensure last snapshot is included
    if (times_buffer_.empty() || times_buffer_.back() < t_final) {
        times_buffer_.push_back(t_final);
        trajectory_buffer_.push_back(final_ions);
    }
    
    // 2. Flush remaining HDF5 data
    if (!times_buffer_.empty()) {
        flush();
    }
    
    // 3. Finalize HDF5 file (writes completion metadata)
    try {
        size_t active_count = std::count_if(
            final_ions.begin(), final_ions.end(),
            [](const IonState& ion) { return ion.active && ion.born; }
        );
        
        io::HDF5Writer::finalize(
            hdf5_filename_,
            true,  // success
            t_final,
            active_count
        );
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to finalize HDF5 file: " 
                  << e.what() << std::endl;
    }
    
    // 4. Write text log completion summary
    if (text_logger_) {
        text_logger_->finalize(final_ions, hdf5_filename_);
    }
    
    initialized_ = false;
}

}  // namespace integrator
}  // namespace ICARION

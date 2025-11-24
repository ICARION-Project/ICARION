// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "OutputManager.h"
#include "core/io/hdf5Writer.h"
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <iomanip>

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
        text_log_file_.open(log_filename_, std::ios::out | std::ios::trunc);
        if (!text_log_file_.is_open()) {
            throw std::runtime_error("Failed to create text log: " + log_filename_);
        }
        
        // Write header banner
        text_log_file_ << "============================================================\n";
        text_log_file_ << "       ICARION Ion Collision And Reaction IntegratiON       \n";
        text_log_file_ << "============================================================\n";
        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        std::tm tm = *std::localtime(&now_c);
        text_log_file_ << "Simulation started: ";
        text_log_file_ << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << "\n\n";
        text_log_file_.flush();
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
    if (text_log_file_.is_open()) {
        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        std::tm tm = *std::localtime(&now_c);
        text_log_file_ << "[";
        text_log_file_ << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
        text_log_file_ << "] " << message << "\n";
        text_log_file_.flush();
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
    
    // Store last time before clearing buffer (for drift-free next_write_time_)
    double last_flushed_time = times_buffer_.back();
    
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
        
        // Update next write time (avoid drift with adaptive timesteps)
        // Use last flushed time as anchor, not incremental addition
        next_write_time_ = last_flushed_time + write_interval_dt_;
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
    if (text_log_file_.is_open()) {
        // Compute ion statistics
        size_t total = final_ions.size();
        size_t active = std::count_if(final_ions.begin(), final_ions.end(),
                                      [](const IonState& i){ return i.active; });
        size_t lost = total - active;
        double frac_lost = total > 0 ? 100.0 * static_cast<double>(lost) / total : 0.0;
        
        // Compute HDF5 file size
        double file_size_MB = 0.0;
        std::ifstream f(hdf5_filename_, std::ifstream::ate | std::ifstream::binary);
        if (f.is_open()) {
            file_size_MB = static_cast<double>(f.tellg()) / (1024.0 * 1024.0);
        }
        
        // Write summary
        text_log_file_ << "------------------------------------------------------------\n";
        text_log_file_ << "Summary:\n";
        text_log_file_ << "  Active ions remaining    : " << active << "\n";
        text_log_file_ << "  Lost ions (boundary)     : " << lost
                       << " (" << std::fixed << std::setprecision(1) << frac_lost << " %)\n";
        text_log_file_ << "  Output file size         : "
                       << std::setprecision(1) << file_size_MB << " MB\n";
        text_log_file_ << "------------------------------------------------------------\n\n";
        
        // Write completion timestamp
        auto now = std::chrono::system_clock::now();
        std::time_t now_c = std::chrono::system_clock::to_time_t(now);
        std::tm tm = *std::localtime(&now_c);
        text_log_file_ << "------------------------------------------------------------\n";
        text_log_file_ << "Simulation finished: ";
        text_log_file_ << std::put_time(&tm, "%Y-%m-%d %H:%M:%S") << "\n";
        text_log_file_ << "------------------------------------------------------------\n";
        text_log_file_ << "Project repository:\n";
        text_log_file_ << "  https://github.com/ICARION-Project/ICARION\n";
        text_log_file_ << "Please cite ICARION once the reference paper is available.\n";
        text_log_file_ << "------------------------------------------------------------\n";
        
        text_log_file_.close();
        
        auto logger = ICARION::log::Logger::main();
        if (logger) {
            logger->info("Text log written to {}", log_filename_);
        } else {
            std::cout << "Log written to " << log_filename_ << std::endl;
        }
    }
    
    initialized_ = false;
}

}  // namespace integrator
}  // namespace ICARION

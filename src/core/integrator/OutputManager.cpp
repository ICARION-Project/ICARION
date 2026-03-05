// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "OutputManager.h"
#include "core/io/hdf5Writer.h"
#include <H5Cpp.h>
#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <chrono>
#include <iomanip>
#include <sstream>
#ifdef _OPENMP
#include <omp.h>
#endif

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

namespace {
// Compose a compact build-info string for reproducibility metadata.
std::string make_build_info() {
    std::ostringstream oss;
    oss << "type=" << BUILD_TYPE;

#ifdef _OPENMP
    const bool openmp_enabled = true;
    const int openmp_threads = omp_get_max_threads();
#else
    const bool openmp_enabled = false;
    const int openmp_threads = 1;
#endif
    oss << "; openmp=" << (openmp_enabled ? "on" : "off");
    oss << "; openmp_threads=" << openmp_threads;

#if defined(ICARION_ENABLE_CUDA) || defined(USE_CUDA)
    oss << "; cuda=on";
#else
    oss << "; cuda=off";
#endif

#ifdef ICARION_BUILD_CORE_ONLY
    oss << "; mode=core-only";
#else
    oss << "; mode=full";
#endif

    return oss.str();
}
}  // namespace

namespace ICARION {
namespace integrator {

OutputManager::OutputManager(
    const std::string& hdf5_filename,
    const std::string& log_filename,
    double write_interval_dt,
    size_t buffer_max
) : hdf5_filename_(hdf5_filename),
    write_interval_dt_(write_interval_dt),
    next_write_time_(0.0),
    buffer_max_(buffer_max),
    log_filename_(log_filename)
{
    if (hdf5_filename_.empty()) {
        throw std::invalid_argument("OutputManager: HDF5 filename cannot be empty");
    }
    
    // Reserve buffer space
    times_buffer_.reserve(buffer_max_);
}

OutputManager::~OutputManager() {
    // Ensure cleanup even if finalize() was never called (crash/exception)
    if (initialized_ && !finalized_ && total_writes_ > 0) {
        try {
            // Flush any remaining buffered data
            if (!times_buffer_.empty()) {
                flush();
            }
            
            // Write incomplete simulation metadata
            // Only if we actually wrote trajectory data (total_writes_ > 0)
            try {
                io::HDF5Writer::finalize(
                    hdf5_filename_,
                    false,      // success = false (incomplete)
                    last_time_, // last known time
                    0           // active_ions unknown (set to 0)
                );
            } catch (const H5::Exception&) {
                // Silently ignore - HDF5 file may not be in writable state
                // This is acceptable since we're in destructor cleanup
            } catch (const std::exception&) {
                // Silently ignore - HDF5 file may not be in writable state
                // This is acceptable since we're in destructor cleanup
            } catch (...) {
                // Silently ignore - destructor must not throw
            }
            
            // Note: Text log may be incomplete, but that's acceptable
        } catch (const std::exception& e) {
            std::cerr << "Warning: OutputManager destructor cleanup failed: " 
                      << e.what() << std::endl;
        }
    }
}

void OutputManager::initialize(
    const config::FullConfig& config,
    const core::IonEnsemble& ensemble
) {
    if (initialized_) {
        throw std::runtime_error("OutputManager: Already initialized");
    }
    
    // 1. Create HDF5 file with metadata using HDF5Writer (SoA)
    try {
        io::HDF5Writer::create_file(
            hdf5_filename_,
            config,
            ensemble,
            GIT_HASH,
            make_build_info()
        );
    } catch (const std::exception& e) {
        throw std::runtime_error("OutputManager: Failed to create HDF5 file (SoA): " + 
                                 std::string(e.what()));
    }
    
    // 2. Initialize text logger (if enabled)
    if (!log_filename_.empty()) {
        std::lock_guard<std::mutex> lock(text_log_mutex_);
        text_log_file_.open(log_filename_, std::ios::out | std::ios::trunc);
        if (!text_log_file_.is_open()) {
            throw std::runtime_error("Failed to create text log: " + log_filename_);
        }
        
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
    species_pool_ = ensemble.species_pool();
}

void OutputManager::log_step(double t, const core::IonEnsemble& ensemble) {
    if (!initialized_) {
        throw std::runtime_error("OutputManager: Must call initialize() before log_step()");
    }
    
    if (finalized_) {
        throw std::runtime_error("OutputManager: Cannot log after finalize()");
    }
    
    last_time_ = t;
    
    // Check if flush needed BEFORE adding 
    if (should_write_before_add(t)) {
        flush();
    }
    
    // Initialize shape on first use
    if (n_ions_ == 0) {
        n_ions_ = ensemble.size();
    } else if (ensemble.size() != n_ions_) {
        throw std::runtime_error("OutputManager: Ion count changed during simulation");
    }

    // Buffer snapshot (flattened)
    times_buffer_.push_back(t);
    const auto* pos_x = ensemble.pos_x_data();
    const auto* pos_y = ensemble.pos_y_data();
    const auto* pos_z = ensemble.pos_z_data();
    const auto* vel_x = ensemble.vel_x_data();
    const auto* vel_y = ensemble.vel_y_data();
    const auto* vel_z = ensemble.vel_z_data();
    const auto* t_ptr = ensemble.time_data();
    const auto* domain_idx = ensemble.domain_index_data();
    const auto* species_idx = ensemble.species_id_indices();

    positions_buffer_.reserve(positions_buffer_.size() + n_ions_ * 3);
    velocities_buffer_.reserve(velocities_buffer_.size() + n_ions_ * 3);
    domain_buffer_.reserve(domain_buffer_.size() + n_ions_);
    species_buffer_.reserve(species_buffer_.size() + n_ions_);
    per_ion_time_buffer_.reserve(per_ion_time_buffer_.size() + n_ions_);

    for (size_t i = 0; i < n_ions_; ++i) {
        positions_buffer_.push_back(pos_x[i]);
        positions_buffer_.push_back(pos_y[i]);
        positions_buffer_.push_back(pos_z[i]);

        velocities_buffer_.push_back(vel_x[i]);
        velocities_buffer_.push_back(vel_y[i]);
        velocities_buffer_.push_back(vel_z[i]);

        domain_buffer_.push_back(domain_idx[i]);
        species_buffer_.push_back(species_idx[i]);
        per_ion_time_buffer_.push_back(t_ptr ? t_ptr[i] : t);
    }
}

void OutputManager::log_progress(const std::string& message) {
    std::lock_guard<std::mutex> lock(text_log_mutex_);
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
    if ((times_buffer_.size() >= buffer_max_) || (t_current >= next_write_time_)) {
        return true;
    }

    if (buffer_byte_cap_ == 0) {
        return false;
    }

    // Rough byte estimate: flattened buffers
    size_t bytes = sizeof(double) * times_buffer_.capacity();
    bytes += sizeof(double) * positions_buffer_.capacity();
    bytes += sizeof(double) * velocities_buffer_.capacity();
    bytes += sizeof(double) * per_ion_time_buffer_.capacity();
    bytes += sizeof(int) * domain_buffer_.capacity();
    bytes += sizeof(uint32_t) * species_buffer_.capacity();
    return bytes >= buffer_byte_cap_;
}

void OutputManager::flush() {
    if (times_buffer_.empty()) {
        return;  // Nothing to flush
    }
    if (buffer_byte_cap_ > 0) {
        size_t bytes = sizeof(double) * times_buffer_.capacity();
        bytes += sizeof(double) * positions_buffer_.capacity();
        bytes += sizeof(double) * velocities_buffer_.capacity();
        bytes += sizeof(double) * per_ion_time_buffer_.capacity();
        bytes += sizeof(int) * domain_buffer_.capacity();
        bytes += sizeof(uint32_t) * species_buffer_.capacity();
        if (bytes > buffer_byte_cap_) {
            throw std::runtime_error("OutputManager: Buffer byte cap exceeded (" + std::to_string(bytes) +
                                     " > " + std::to_string(buffer_byte_cap_) + ")");
        }
    }
    
    // Store last time before clearing buffer (for drift-free next_write_time_)
    double last_flushed_time = times_buffer_.back();
    
    try {
        // Write all buffered snapshots in ONE batch operation
        io::HDF5Writer::append_trajectory_batch_flat(
            hdf5_filename_,
            times_buffer_,
            n_ions_,
            positions_buffer_,
            velocities_buffer_,
            domain_buffer_,
            species_buffer_,
            species_pool_,
            per_ion_time_buffer_
        );
        
        // Clear buffers
        times_buffer_.clear();
        positions_buffer_.clear();
        velocities_buffer_.clear();
        per_ion_time_buffer_.clear();
        domain_buffer_.clear();
        species_buffer_.clear();
        
        // Update next write time (avoid drift with adaptive timesteps)
        // Use last flushed time as anchor, not incremental addition
        next_write_time_ = last_flushed_time + write_interval_dt_;
        total_writes_++;
        
    } catch (const std::exception& e) {
        throw std::runtime_error("OutputManager: HDF5 flush failed: " + 
                                 std::string(e.what()));
    }
}

void OutputManager::finalize(double t_final, const core::IonEnsemble& final_ensemble) {
    if (!initialized_) {
        return;
    }

    // Ensure last snapshot present
    if (times_buffer_.empty() || times_buffer_.back() < t_final) {
        log_step(t_final, final_ensemble);
    }

    // Flush remaining data
    if (!times_buffer_.empty()) {
        flush();
    }

    // Update death times (uses single conversion inside writer)
    try {
        io::HDF5Writer::update_death_times(hdf5_filename_, final_ensemble);
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to update death times in HDF5 file: " 
                  << e.what() << std::endl;
    }

    // Finalize HDF5
    try {
        size_t active_count = 0;
        const auto* active = final_ensemble.active_data();
        const auto* born = final_ensemble.born_data();
        for (size_t i = 0; i < final_ensemble.size(); ++i) {
            if (active[i] && born[i]) {
                ++active_count;
            }
        }

        io::HDF5Writer::finalize(
            hdf5_filename_,
            true,
            t_final,
            active_count
        );
    } catch (const H5::Exception& e) {
        std::cerr << "Warning: Failed to finalize HDF5 file: "
                  << e.getCDetailMsg() << std::endl;
    } catch (const std::exception& e) {
        std::cerr << "Warning: Failed to finalize HDF5 file: "
                  << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Warning: Failed to finalize HDF5 file (unknown error)" << std::endl;
    }

    // Text log completion summary
    {
        std::lock_guard<std::mutex> lock(text_log_mutex_);
        if (text_log_file_.is_open()) {
            size_t total = final_ensemble.size();
            size_t active = 0;
            const auto* active_ptr = final_ensemble.active_data();
            const auto* born_ptr = final_ensemble.born_data();
            for (size_t i = 0; i < total; ++i) {
                if (active_ptr[i] && born_ptr[i]) active++;
            }
            size_t lost = total - active;
            
            text_log_file_ << "============================================================\n";
            text_log_file_ << "Simulation completed successfully\n";
            text_log_file_ << "Final time: " << t_final << " s\n";
            text_log_file_ << "Active ions: " << active << " / " << total << "\n";
            text_log_file_ << "Lost ions:   " << lost << "\n";
            text_log_file_ << "============================================================\n";
            
            text_log_file_.close();
        }
    }

    finalized_ = true;
    initialized_ = false;
}

}  // namespace integrator
}  // namespace ICARION

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "OutputManager.h"
#include "core/io/hdf5Writer.h"  // Legacy HDF5Writer (TODO: migrate to v2 in Phase 5B)
#include "core/config/adapter/LegacyAdapter.h"  // Temporary for GlobalParams conversion
#include "H5Cpp.h"
#include <iostream>
#include <stdexcept>

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
    
    // 1. Create HDF5 file and write metadata
    try {
        H5::H5File file(hdf5_filename_, H5F_ACC_TRUNC);
        
        // Convert FullConfig → GlobalParams (temporary until Phase 5B)
        // TODO: Extend write_params_to_HDF5 to accept FullConfig directly
        auto legacy_params = config::LegacyAdapter::to_global_params(config);
        auto legacy_domains = config::LegacyAdapter::to_instrument_domains(config);
        
        ICARION::io::write_params_to_HDF5(file, legacy_params, legacy_domains);
        
        // Write species metadata
        if (!ions.empty()) {
            ICARION::io::write_species_metadata(file, ions);
        }
        
        file.close();
        
    } catch (const H5::Exception& e) {
        throw std::runtime_error("OutputManager: Failed to create HDF5 file: " + 
                                 std::string(e.getDetailMsg()));
    }
    
    // 2. Initialize text logger (if enabled)
    if (!log_filename_.empty()) {
        text_logger_ = std::make_unique<io::RunLogger>(log_filename_);
        text_logger_->writeHeader();
        
        // Write params (requires SpeciesDatabase + ReactionList - TODO: integrate in Phase 5B)
        // For now, skip writeGlobalParams() (requires full species/reaction data)
        // Only log progress messages via log()
        
        // Write domain info
        auto legacy_domains = config::LegacyAdapter::to_instrument_domains(config);
        text_logger_->writeInstrumentDomains(legacy_domains);
    }
    
    initialized_ = true;
    next_write_time_ = write_interval_dt_;
}

void OutputManager::log_step(double t, const std::vector<IonState>& ions) {
    if (!initialized_) {
        throw std::runtime_error("OutputManager: Not initialized (call initialize() first)");
    }
    
    // Buffer snapshot
    times_buffer_.push_back(t);
    trajectory_buffer_.push_back(ions);
    
    // Auto-flush if needed
    if (should_write(t)) {
        flush();
    }
}

void OutputManager::log_progress(const std::string& message) {
    if (text_logger_) {
        text_logger_->log(message);
    }
}

bool OutputManager::should_write(double t_current) const {
    // Flush if buffer full OR time interval exceeded
    return (times_buffer_.size() >= buffer_max_) || (t_current >= next_write_time_);
}

void OutputManager::flush() {
    if (times_buffer_.empty()) {
        return;  // Nothing to flush
    }
    
    try {
        // Write to HDF5
        ICARION::io::append_to_HDF5(hdf5_filename_, times_buffer_, trajectory_buffer_);
        
        // Clear buffers
        times_buffer_.clear();
        trajectory_buffer_.clear();
        
        // Update next write time
        next_write_time_ += write_interval_dt_;
        total_writes_++;
        
    } catch (const H5::Exception& e) {
        throw std::runtime_error("OutputManager: HDF5 flush failed: " + 
                                 std::string(e.getDetailMsg()));
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
    
    // 3. Write HDF5 completion metadata
    try {
        H5::H5File file(hdf5_filename_, H5F_ACC_RDWR);
        
        // Count active ions
        size_t active_count = std::count_if(
            final_ions.begin(), final_ions.end(),
            [](const IonState& ion) { return ion.active && ion.born; }
        );
        
        // Write completion flag
        hsize_t dims[1] = {1};
        H5::DataSpace scalar_space(1, dims);
        H5::Attribute completion_attr = file.createAttribute(
            "simulation_completed", H5::PredType::NATIVE_HBOOL, scalar_space);
        hbool_t completed = true;
        completion_attr.write(H5::PredType::NATIVE_HBOOL, &completed);
        
        // Write active ions count
        H5::Attribute active_attr = file.createAttribute(
            "final_active_ions", H5::PredType::NATIVE_HSIZE, scalar_space);
        active_attr.write(H5::PredType::NATIVE_HSIZE, &active_count);
        
        file.close();
        
    } catch (const H5::Exception& e) {
        std::cerr << "Warning: Failed to write HDF5 completion metadata: " 
                  << e.getDetailMsg() << std::endl;
    }
    
    // 4. Write text log completion summary
    if (text_logger_) {
        text_logger_->finalize(final_ions, hdf5_filename_);
    }
    
    initialized_ = false;
}

}  // namespace integrator
}  // namespace ICARION

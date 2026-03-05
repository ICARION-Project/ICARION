// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_CONFIG_OUTPUT_CONFIG_H
#define ICARION_CONFIG_OUTPUT_CONFIG_H

#include "../validation/ValidationResult.h"
#include <string>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief Output configuration
 * 
 * Controls what data is written and where.
 * Replaces GlobalParams (output-related parts only).
 */
struct OutputConfig {
    // === Output paths ===
    std::string folder = "./results";           ///< Output directory
    std::string trajectory_file = "trajectories.h5";  ///< HDF5 trajectory file
    
    // === Output control ===
    bool print_progress = true;                 ///< Print simulation progress
    size_t buffer_byte_cap = 0;                 ///< Optional memory cap for trajectory buffer (bytes, 0=disabled)
    
    // Future extensions:
    // bool save_reactions = false;
    // bool save_collision_events = false;
    // std::string log_file = "icarion.log";
    
    /**
     * @brief Validate output configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        if (folder.empty()) {
            result.add_error("Output folder cannot be empty");
        }
        if (trajectory_file.empty()) {
            result.add_error("Trajectory file name cannot be empty");
        }
        
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_OUTPUT_CONFIG_H

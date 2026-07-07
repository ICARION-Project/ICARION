// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_CONFIG_OUTPUT_CONFIG_H
#define ICARION_CONFIG_OUTPUT_CONFIG_H

#include "OutputEnums.h"
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
    struct DeepAnalysisConfig {
        DeepAnalysisMode mode_type = DeepAnalysisMode::Off;
        std::string mode = "off";
        int domain_filter_index = -1;           ///< -1 = all domains
        size_t sample_every_n = 10;             ///< sampled_events: keep every Nth collision event
        size_t max_events_per_ion = 0;          ///< 0 = unlimited
    };

    // === Output paths ===
    std::string folder = "./results";           ///< Output directory
    std::string trajectory_file = "trajectories.h5";  ///< HDF5 trajectory file
    
    // === Output control ===
    std::string trajectory_mode = "full";      ///< Trajectory output mode: full|minimal
    bool print_progress = true;                 ///< Print simulation progress
    size_t buffer_byte_cap = 0;                 ///< Optional memory cap for trajectory buffer (bytes, 0=disabled)
    DeepAnalysisConfig deep_analysis;           ///< Optional collision diagnostics output
    
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
        if (trajectory_mode != "full" && trajectory_mode != "minimal") {
            result.add_error("output.trajectory_mode must be one of: full, minimal");
        }
        if (deep_analysis.sample_every_n == 0) {
            result.add_error("output.deep_analysis.sample_every_n must be >= 1");
        }
        if (deep_analysis.domain_filter_index < -1) {
            result.add_error("output.deep_analysis.domain_filter_index must be -1 or >= 0");
        }
        
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_OUTPUT_CONFIG_H

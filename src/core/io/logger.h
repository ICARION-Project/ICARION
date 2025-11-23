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
 *   @file        logger.h
 *   @brief       Logs simulation parameters.
 *
 *   @details
 *   Provides functions to log simulation parameters separated by 
 *   global and instrument domain settings. Supports runtime appends.
 *
 *
 *   @date        2025-10-16
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#pragma once
#include <fstream>
#include <string>
#include <vector>
#include "core/types/IonState.h"
#include "core/io/speciesLoader.h"

namespace ICARION {
namespace io {

/**
 * @brief Simulation run logger for parameter documentation and timestamping
 * 
 * Logs simulation configuration, ion states, and final results to text file.
 */
class RunLogger {
public:
    explicit RunLogger(const std::string& output_file_base);
    
    /** @brief Write header banner and timestamp */
    void writeHeader();
    
    // Note: writeGlobalParams() and writeInstrumentDomains() removed (legacy)
    
    /** @brief Log arbitrary message with timestamp */
    void log(const std::string& msg);
    
    /** @brief Write final simulation results and close log file */
    void finalize(const std::vector<IonState>& ions, const std::string& output_file);

private:
    std::ofstream file_;
    std::string filepath_;
    std::string timestamp() const;
};

/**
 * @brief Debug logger respecting log level configuration
 * @param msg Message to log
 * 
 * Used by core computation modules for debugging output (e.g., geometry loading,
 * collision handler creation). Only outputs when log level is DEBUG.
 * 
 * **Usage:** Enable with `--verbose` or `--log-level DEBUG` CLI flags.
 * 
 * **Example:**
 * ```bash
 * icarion --verbose config.json        # Shows debug output
 * icarion --log-level INFO config.json # Hides debug output (default)
 * ```
 */
void debug_log(const std::string& msg);

}  // namespace io
}  // namespace ICARION

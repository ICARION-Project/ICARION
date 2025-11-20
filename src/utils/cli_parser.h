// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   Command-line argument parser
 *
 *   @file        cli_parser.h
 *   @brief       Parse and validate command-line arguments
 *
 *   @date        2025-11-10
 *   @version     1.0.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#pragma once

#include <string>
#include <optional>
#include <map>
#include <vector>

namespace ICARION {
namespace cli {

/**
 * @brief Command-line options for ICARION
 * 
 * Extended with logging, output control, and config overrides.
 */
struct CLIOptions {
    // === Core options ===
    std::string config_file;              ///< Path to JSON configuration file
    std::optional<unsigned int> seed;     ///< RNG seed override (--seed)
    bool dry_run{false};                  ///< Validate config only (--dry-run)
    bool no_reactions{false};             ///< Disable reactions (--no-reactions)
    bool show_help{false};                ///< Show help message (--help)
    bool show_version{false};             ///< Show version (--version)
    bool validate_config{false};          ///< Validate config and exit (--validate-config)
    
    // === Logging options (Phase 1: ACTIVE) ===
    std::string log_level{"INFO"};        ///< Log level: DEBUG, INFO, WARN, ERROR
    std::optional<std::string> log_file;  ///< Log to file instead of console
    bool verbose{false};                  ///< Verbose mode (alias for --log-level DEBUG)
    
    // === Output control (Phase 1: ACTIVE) ===
    std::optional<std::string> output_file;  ///< Override output HDF5 filename
    std::optional<std::string> output_dir;   ///< Override output directory
    
    // === Config overrides (Phase 1: ACTIVE) ===
    std::map<std::string, std::string> overrides;  ///< Config key-value overrides (--set)
    
    // === Performance/Debug (TODO - Phase 2) ===
    bool benchmark{false};       ///< Print timing statistics (TODO)
    bool profile{false};         ///< Enable profiling (TODO)
    bool check_nan{false};       ///< Enable NaN/Inf checks (TODO)
    
    // === Information flags (Phase 1: ACTIVE) ===
    bool dump_build_info{false};      ///< Show detailed build information
    bool dump_hdf5_schema{false};     ///< Show HDF5 output schema
    bool dump_config_schema{false};   ///< Export JSON config schema
    bool list_collision_models{false}; ///< List available collision models
    bool list_integrators{false};     ///< List available integrators
};

/**
 * @brief Parse command-line arguments using cxxopts
 * 
 * @param argc Argument count
 * @param argv Argument vector
 * @return Parsed options
 */
CLIOptions parse_arguments(int argc, char* argv[]);

/**
 * @brief Print help message
 * 
 * @param program_name Name of the executable (argv[0])
 */
void print_help(const std::string& program_name);

/**
 * @brief Print version information
 */
void print_version();

/**
 * @brief Print detailed build information
 */
void print_build_info();

/**
 * @brief Print HDF5 output schema documentation
 */
void print_hdf5_schema();

/**
 * @brief Print JSON config schema (or path to schemas)
 */
void print_config_schema();

/**
 * @brief List available collision models
 */
void list_collision_models();

/**
 * @brief List available integrators
 */
void list_integrators();

}  // namespace cli
}  // namespace ICARION

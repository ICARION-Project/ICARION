// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once

#include <string>
#include <optional>
#include <map>
#include <vector>
#include <cstdint>

namespace ICARION {
namespace cli {

/**
 * @brief Command-line options for ICARION
 * 
 * CLI surface primarily covers config selection, logging, and a few info flags. GPU use\n+ * is driven by config; there is no CLI toggle beyond the JSON settings.
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
    
    // === Logging options ===
    std::string log_level{"INFO"};        ///< Log level: DEBUG, INFO, WARN, ERROR
    std::optional<std::string> log_file;  ///< Log to file instead of console
    std::string log_format{"text"};       ///< Log format: text or json
    bool verbose{false};                  ///< Verbose mode (alias for --log-level DEBUG)
    
    // === Output control ===
    std::optional<std::string> output_file;  ///< Override output HDF5 filename
    std::optional<std::string> output_dir;   ///< Override output directory
    std::optional<uint64_t> buffer_byte_cap; ///< Cap for trajectory buffer (bytes, 0 = unlimited)
    
    // === Config overrides ===
    std::map<std::string, std::string> overrides;  ///< Config key-value overrides (--set)
    
    // === Information flags ===
    bool dump_build_info{false};      ///< Show detailed build information
    bool dump_hdf5_schema{false};     ///< Show HDF5 output schema
    bool dump_config_schema{false};   ///< Export JSON config schema
    bool list_collision_models{false}; ///< List available collision models
    bool list_integrators{false};     ///< List available integrators
    bool validate_schema{false};      ///< Validate config against schema using validator.py
    bool check_deps{false};           ///< Verify all dependencies and versions
    
    // === Performance/Profiling options ===
    bool benchmark{false};                      ///< Enable detailed timing statistics (--benchmark)
    bool profile{false};                        ///< Enable profiling instrumentation (--profile)
    std::optional<std::string> profile_output;  ///< Profile output file (--profile-output)
    std::optional<int> num_threads;             ///< Number of OpenMP threads (--threads)
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

/**
 * @brief Validate config file against JSON schema
 */
void validate_schema(const std::string& config_file);

/**
 * @brief Check all dependencies and their versions
 */
void check_dependencies();

}  // namespace cli
}  // namespace ICARION

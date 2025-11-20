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

namespace ICARION {
namespace cli {

/**
 * @brief Command-line options for ICARION
 */
struct CLIOptions {
    std::string config_file;              ///< Path to JSON configuration file
    std::optional<unsigned int> seed;     ///< RNG seed override (--seed)
    bool dry_run{false};                  ///< Validate config only (--dry-run)
    bool no_reactions{false};             ///< Disable reactions (--no-reactions)
    bool show_help{false};                ///< Show help message (--help)
    bool show_version{false};             ///< Show version (--version)
};

/**
 * @brief Parse command-line arguments
 * 
 * @param argc Argument count
 * @param argv Argument vector
 * @return Parsed options
 * 
 * Supported flags:
 * - --help: Show usage information
 * - --version: Show ICARION version
 * - --seed N: Override RNG seed from config
 * - --dry-run: Validate configuration without running simulation
 * - --no-reactions: Disable chemical reactions
 * 
 * Positional argument: configuration file path
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

}  // namespace cli
}  // namespace ICARION

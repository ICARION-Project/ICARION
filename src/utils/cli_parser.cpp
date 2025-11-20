/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   Command-line argument parser - Implementation
 *
 *   @file        cli_parser.cpp
 *   @brief       Parse and validate command-line arguments
 *
 *   @date        2025-11-10
 *   @version     1.0.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include "cli_parser.h"
#include <iostream>
#include <cstring>
#include <cstdlib>

#ifndef ICARION_VERSION
#define ICARION_VERSION "1.0.0"
#endif

#ifndef GIT_HASH
#define GIT_HASH "unknown"
#endif

namespace ICARION {
namespace cli {

void print_help(const std::string& program_name) {
    std::cout << "ICARION - Ion Collision And Reaction IntegratiON\n";
    std::cout << "Version " << ICARION_VERSION << " (commit " << GIT_HASH << ")\n";
    std::cout << "\n";
    std::cout << "Usage: " << program_name << " [OPTIONS] <config.json>\n";
    std::cout << "\n";
    std::cout << "Positional Arguments:\n";
    std::cout << "  config.json          Path to JSON configuration file\n";
    std::cout << "\n";
    std::cout << "Options:\n";
    std::cout << "  --help               Show this help message and exit\n";
    std::cout << "  --version            Show version information and exit\n";
    std::cout << "  --seed N             Override RNG seed from config (for reproducibility)\n";
    std::cout << "  --dry-run            Validate configuration without running simulation\n";
    std::cout << "  --no-reactions       Disable chemical reactions (collisions only)\n";
    std::cout << "  --validate-config    Validate config file and display warnings/errors\n";
    std::cout << "\n";
    std::cout << "Examples:\n";
    std::cout << "  " << program_name << " params.json\n";
    std::cout << "  " << program_name << " --seed 12345 params.json\n";
    std::cout << "  " << program_name << " --dry-run params.json\n";
    std::cout << "\n";
    std::cout << "For more information, see README.md or visit:\n";
    std::cout << "  https://github.com/chsch95/numerical_model\n";
}

void print_version() {
    std::cout << "ICARION v" << ICARION_VERSION << "\n";
    std::cout << "Git commit: " << GIT_HASH << "\n";
    std::cout << "Build type: ";
#ifdef ICARION_BUILD_CORE_ONLY
    std::cout << "Core-Only (no FieldSolver)\n";
#else
    std::cout << "Full (with FieldSolver & Optimizer)\n";
#endif
#ifdef USE_CUDA
    std::cout << "GPU acceleration: Enabled\n";
#else
    std::cout << "GPU acceleration: Disabled\n";
#endif
}

CLIOptions parse_arguments(int argc, char* argv[]) {
    CLIOptions opts;
    
    // Parse flags
    int config_arg_index = -1;
    for (int i = 1; i < argc; ++i) {
        std::string arg(argv[i]);
        
        if (arg == "--help" || arg == "-h") {
            opts.show_help = true;
            return opts;  // Early return for help
        }
        else if (arg == "--version" || arg == "-v") {
            opts.show_version = true;
            return opts;  // Early return for version
        }
        else if (arg == "--seed") {
            if (i + 1 < argc) {
                opts.seed = static_cast<unsigned int>(std::atoi(argv[i + 1]));
                ++i;  // Skip next argument
            } else {
                std::cerr << "Error: --seed requires an argument\n";
                std::exit(1);
            }
        }
        else if (arg == "--dry-run") {
            opts.dry_run = true;
        }
        else if (arg == "--no-reactions") {
            opts.no_reactions = true;
        }
        else if (arg == "--validate-config") {
            opts.validate_config = true;
        }
        else if (arg[0] == '-') {
            std::cerr << "Error: Unknown option '" << arg << "'\n";
            std::cerr << "Use --help for usage information\n";
            std::exit(1);
        }
        else {
            // Positional argument (config file)
            if (config_arg_index == -1) {
                config_arg_index = i;
                opts.config_file = arg;
            } else {
                std::cerr << "Error: Multiple configuration files specified\n";
                std::cerr << "  First: " << argv[config_arg_index] << "\n";
                std::cerr << "  Second: " << arg << "\n";
                std::exit(1);
            }
        }
    }
    
    // Validate: config file required (unless --help or --version)
    if (opts.config_file.empty() && !opts.show_help && !opts.show_version) {
        std::cerr << "Error: Missing configuration file\n";
        std::cerr << "Usage: " << argv[0] << " [OPTIONS] <config.json>\n";
        std::cerr << "Use --help for more information\n";
        std::exit(1);
    }
    
    return opts;
}

}  // namespace cli
}  // namespace ICARION

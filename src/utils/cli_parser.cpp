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
#include <cxxopts.hpp>
#include <iostream>
#include <cstdlib>

#ifndef ICARION_VERSION
#define ICARION_VERSION "1.0.0"
#endif

#ifndef GIT_HASH
#define GIT_HASH "unknown"
#endif

namespace ICARION {
namespace cli {

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

void print_help(const std::string& program_name) {
    // Delegated to cxxopts (called in parse_arguments)
    std::cout << "Use --help to see available options\n";
}

CLIOptions parse_arguments(int argc, char* argv[]) {
    CLIOptions opts;
    
    // Create cxxopts parser
    cxxopts::Options parser("icarion", "Ion Collision And Reaction IntegratiON");
    parser.allow_unrecognised_options();
    
    // === Define all options ===
    parser.add_options("Core")
        ("h,help", "Show this help message and exit")
        ("v,version", "Show version information and exit")
        ("c,config", "Path to JSON configuration file", 
         cxxopts::value<std::string>(), "FILE")
        ("seed", "Override RNG seed from config (for reproducibility)", 
         cxxopts::value<unsigned int>(), "N")
        ("dry-run", "Validate configuration without running simulation")
        ("validate-config", "Validate config file and display warnings/errors")
        ("no-reactions", "Disable chemical reactions (collisions only)");
    
    // === Logging options (Phase 1: ACTIVE) ===
    parser.add_options("Logging")
        ("log-level", "Set log level (DEBUG|INFO|WARN|ERROR) [default: INFO]", 
         cxxopts::value<std::string>()->default_value("INFO"), "LEVEL")
        ("log-file", "Write logs to FILE instead of console", 
         cxxopts::value<std::string>(), "FILE")
        ("verbose", "Enable verbose output (alias for --log-level DEBUG)");
    
    // === Output options (Phase 1: ACTIVE) ===
    parser.add_options("Output")
        ("o,output", "Override output HDF5 filename", 
         cxxopts::value<std::string>(), "FILE")
        ("output-dir", "Override output directory", 
         cxxopts::value<std::string>(), "DIR");
    
    // === Config overrides (Phase 1: ACTIVE) ===
    parser.add_options("Advanced")
        ("set", "Override config value (e.g., --set simulation.dt_s=1e-10)", 
         cxxopts::value<std::vector<std::string>>(), "KEY=VALUE")
        ("benchmark", "[TODO] Print timing statistics for each simulation phase")
        ("profile", "[TODO] Enable profiling (requires profiler build)")
        ("check-nan", "[TODO] Enable NaN/Inf checks in integrator");
    
    // Parse
    cxxopts::ParseResult result;
    try {
        result = parser.parse(argc, argv);
    } catch (const cxxopts::exceptions::exception& e) {
        std::cerr << "Error parsing arguments: " << e.what() << "\n";
        std::cerr << "Use --help for usage information\n";
        std::exit(1);
    }
    
    // === Handle --help and --version early ===
    if (result.count("help")) {
        std::cout << parser.help({"Core", "Logging", "Output", "Advanced"}) << "\n";
        std::cout << "\nExamples:\n";
        std::cout << "  icarion config.json\n";
        std::cout << "  icarion --verbose config.json\n";
        std::cout << "  icarion --log-level DEBUG --log-file debug.log config.json\n";
        std::cout << "  icarion --output custom.h5 config.json\n";
        std::cout << "  icarion --set simulation.dt_s=1e-10 config.json\n";
        std::cout << "  icarion --dry-run --validate-config config.json\n";
        std::cout << "\nFor more information, see README.md\n";
        opts.show_help = true;
        return opts;
    }
    
    if (result.count("version")) {
        print_version();
        opts.show_version = true;
        return opts;
    }
    
    // === Parse core options ===
    if (result.count("config")) {
        opts.config_file = result["config"].as<std::string>();
    } else {
        // Try positional argument (first non-flag argument)
        if (result.unmatched().size() > 0) {
            opts.config_file = result.unmatched()[0];
        } else {
            std::cerr << "Error: Missing configuration file\n";
            std::cerr << "Usage: icarion [OPTIONS] <config.json>\n";
            std::cerr << "Use --help for more information\n";
            std::exit(1);
        }
    }
    
    if (result.count("seed")) {
        opts.seed = result["seed"].as<unsigned int>();
    }
    
    opts.dry_run = result.count("dry-run") > 0;
    opts.validate_config = result.count("validate-config") > 0;
    opts.no_reactions = result.count("no-reactions") > 0;
    
    // === Parse logging options (Phase 1: ACTIVE) ===
    if (result.count("log-level")) {
        opts.log_level = result["log-level"].as<std::string>();
        
        // Validate log level
        if (opts.log_level != "DEBUG" && opts.log_level != "INFO" &&
            opts.log_level != "WARN" && opts.log_level != "ERROR") {
            std::cerr << "Error: Invalid log level '" << opts.log_level << "'\n";
            std::cerr << "Valid options: DEBUG, INFO, WARN, ERROR\n";
            std::exit(1);
        }
    }
    
    if (result.count("log-file")) {
        opts.log_file = result["log-file"].as<std::string>();
    }
    
    if (result.count("verbose")) {
        opts.verbose = true;
        opts.log_level = "DEBUG";  // Verbose implies DEBUG
    }
    
    // === Parse output options (Phase 1: ACTIVE) ===
    if (result.count("output")) {
        opts.output_file = result["output"].as<std::string>();
    }
    
    if (result.count("output-dir")) {
        opts.output_dir = result["output-dir"].as<std::string>();
    }
    
    // === Parse config overrides (Phase 1: ACTIVE) ===
    if (result.count("set")) {
        auto overrides = result["set"].as<std::vector<std::string>>();
        for (const auto& override : overrides) {
            size_t eq_pos = override.find('=');
            if (eq_pos == std::string::npos) {
                std::cerr << "Error: Invalid --set format (expected KEY=VALUE): " << override << "\n";
                std::exit(1);
            }
            std::string key = override.substr(0, eq_pos);
            std::string value = override.substr(eq_pos + 1);
            opts.overrides[key] = value;
        }
    }
    
    // === Parse performance/debug flags (Phase 2: TODO) ===
    opts.benchmark = result.count("benchmark") > 0;
    opts.profile = result.count("profile") > 0;
    opts.check_nan = result.count("check-nan") > 0;
    
    if (opts.benchmark || opts.profile || opts.check_nan) {
        std::cerr << "Warning: Performance/debug flags are not yet implemented (TODO Phase 2)\n";
    }
    
    return opts;
}

}  // namespace cli
}  // namespace ICARION

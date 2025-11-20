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
#include <fstream>
#include <cstdlib>
#include <hdf5.h>
#include <unistd.h>

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
    
    // === Information flags (Phase 1: ACTIVE) ===
    parser.add_options("Information")
        ("dump-build-info", "Show detailed build configuration and features")
        ("dump-hdf5-schema", "Show HDF5 output schema documentation")
        ("dump-config-schema", "Show path to JSON config schemas")
        ("list-collision-models", "List available collision models")
        ("list-integrators", "List available integrators")
        ("validate-schema", "Validate config against schema (uses validator.py)")
        ("check-deps", "Verify all dependencies and their versions");
    
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
        std::cout << parser.help({"Core", "Logging", "Output", "Advanced", "Information"}) << "\n";
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
    
    // === Handle information flags early (don't require config file) ===
    opts.dump_build_info = result.count("dump-build-info") > 0;
    opts.dump_hdf5_schema = result.count("dump-hdf5-schema") > 0;
    opts.dump_config_schema = result.count("dump-config-schema") > 0;
    opts.list_collision_models = result.count("list-collision-models") > 0;
    opts.list_integrators = result.count("list-integrators") > 0;
    opts.check_deps = result.count("check-deps") > 0;
    
    // validate-schema needs config file, so parse it separately
    opts.validate_schema = result.count("validate-schema") > 0;
    
    if (opts.dump_build_info || opts.dump_hdf5_schema || opts.dump_config_schema ||
        opts.list_collision_models || opts.list_integrators || opts.check_deps) {
        // These info flags don't require config file
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

void print_build_info() {
    std::cout << "=== ICARION Build Information ===\n";
    std::cout << "Version:           " << ICARION_VERSION << "\n";
    std::cout << "Git Commit:        " << GIT_HASH << "\n";
    std::cout << "Build Date:        " << __DATE__ << " " << __TIME__ << "\n";
    
#ifdef CMAKE_BUILD_TYPE
    std::cout << "Build Type:        " << CMAKE_BUILD_TYPE << "\n";
#else
    std::cout << "Build Type:        Unknown\n";
#endif
    
    std::cout << "\n=== Features ===\n";
#ifdef USE_CUDA
    std::cout << "GPU Acceleration:  Enabled (CUDA)\n";
#else
    std::cout << "GPU Acceleration:  Disabled\n";
#endif

#ifdef _OPENMP
    std::cout << "OpenMP:            Enabled\n";
#else
    std::cout << "OpenMP:            Disabled\n";
#endif

#ifdef ICARION_BUILD_CORE_ONLY
    std::cout << "Build Mode:        Core-Only (no FieldSolver)\n";
#else
    std::cout << "Build Mode:        Full (with FieldSolver & Optimizer)\n";
#endif

    std::cout << "\n=== Dependencies ===\n";
    std::cout << "HDF5:              Available\n";
    std::cout << "JsonCpp:           Available\n";
    std::cout << "cxxopts:           3.1.1\n";
    
    std::cout << "\n=== Compiler ===\n";
#if defined(__GNUC__) && !defined(__clang__)
    std::cout << "Compiler:          GCC " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
#elif defined(__clang__)
    std::cout << "Compiler:          Clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__ << "\n";
#elif defined(_MSC_VER)
    std::cout << "Compiler:          MSVC " << _MSC_VER << "\n";
#else
    std::cout << "Compiler:          Unknown\n";
#endif
    
#if __cplusplus == 202002L
    std::cout << "C++ Standard:      C++20\n";
#elif __cplusplus == 201703L
    std::cout << "C++ Standard:      C++17\n";
#elif __cplusplus == 201402L
    std::cout << "C++ Standard:      C++14\n";
#else
    std::cout << "C++ Standard:      " << __cplusplus << "\n";
#endif
    
    std::cout << "\n=== Paths ===\n";
    std::cout << "Schema Directory:  ./schema/\n";
    std::cout << "Examples:          ./examples/\n";
}

void print_hdf5_schema() {
    std::cout << R"(
ICARION HDF5 Output Schema v1.0
================================

Structure:
----------

/metadata/
  version         string      ICARION version
  git_hash        string      Git commit hash
  timestamp       string      ISO8601 timestamp (YYYY-MM-DDTHH:MM:SSZ)
  rng_seed        int64       Random number generator seed
  hostname        string      Machine hostname
  username        string      Username
  config_file     string      Path to configuration file

/system_info/
  os              string      Operating system
  cpu_model       string      CPU model name
  cpu_cores       int32       Number of CPU cores
  memory_gb       float64     Total system memory [GB]

/trajectory/
  time            float64[N]      Time points [s]
  position        float64[N,3]    Ion positions [x,y,z] [m]
  velocity        float64[N,3]    Ion velocities [vx,vy,vz] [m/s]
  species_id      int32[N]        Species identifier
  active          bool[N]         Ion active flag (true=active, false=lost)
  domain_index    int32[N]        Current domain index
  ion_id          int32[N]        Unique ion identifier

/species/<name>/
  mass_kg         float64     Ion mass [kg]
  charge_C        float64     Ion charge [C]
  mobility_m2Vs   float64     Ion mobility [m²/(V·s)]
  ccs_m2          float64     Collision cross-section [m²]

/config/
  full_config     string      Complete JSON configuration (serialized)

/domains/
  domain_<N>/
    name          string      Domain name
    instrument    string      Instrument type (IMS, LQIT, TOF, etc.)
    length_m      float64     Domain length [m]
    radius_m      float64     Domain radius [m]

Notes:
------
  - N = number of time snapshots × number of ions
  - All quantities in SI units unless specified
  - Use h5dump, HDFView, or Python h5py to read files
  - Trajectory data stored in chunked + compressed format

Example (Python):
-----------------
  import h5py
  with h5py.File('output.h5', 'r') as f:
      time = f['/trajectory/time'][:]
      pos = f['/trajectory/position'][:]
      metadata = dict(f['/metadata'].attrs)

For detailed documentation, see: docs/OUTPUT_SCHEMA.md
)";
}

void print_config_schema() {
    std::cout << "JSON Configuration Schemas:\n";
    std::cout << "===========================\n\n";
    std::cout << "Location: ./schema/\n\n";
    std::cout << "Available schemas:\n";
    std::cout << "  - icarion-config.schema.json    Main configuration schema\n";
    std::cout << "  - simulation.schema.json        Simulation parameters\n";
    std::cout << "  - physics.schema.json           Physics settings\n";
    std::cout << "  - domain.schema.json            Domain configuration\n";
    std::cout << "  - fields.schema.json            Field definitions\n";
    std::cout << "  - environment.schema.json       Environment parameters\n";
    std::cout << "  - geometry.schema.json          Geometry specification\n";
    std::cout << "  - ions.schema.json              Ion initialization\n";
    std::cout << "  - output.schema.json            Output configuration\n";
    std::cout << "  - common-types.schema.json      Common type definitions\n";
    std::cout << "\nValidation:\n";
    std::cout << "  python schema/validator.py config.json\n";
    std::cout << "\nDocumentation:\n";
    std::cout << "  See: docs/CONFIG_GUIDE.md\n";
    std::cout << "       docs/INPUT_SCHEMA.md\n";
}

void list_collision_models() {
    std::cout << "Available Collision Models:\n";
    std::cout << "===========================\n\n";
    std::cout << "  NoCollisions   - No collision modeling (free flight)\n";
    std::cout << "  HardSphere     - Hard sphere collisions (deterministic)\n";
    std::cout << "  Langevin       - Langevin collisions (polarization)\n";
    std::cout << "  Friction       - Frictional drag\n";
    std::cout << "  HSS            - Hard sphere stochastical\n";
    std::cout << "  EHSS           - Enhanced hard sphere statistical\n";
    std::cout << "\nUsage:\n";
    std::cout << "  In config:   \"collision_model\": \"EHSS\"\n";
    std::cout << "  Via CLI:     --set physics.collision_model=EHSS\n";
    std::cout << "\nDocumentation:\n";
    std::cout << "  See: docs/COLLISION_MODELS.md (if exists)\n";
}

void list_integrators() {
    std::cout << "Available Integrators:\n";
    std::cout << "======================\n\n";
    std::cout << "  RK4      - 4th-order Runge-Kutta (fixed timestep)\n";
    std::cout << "           - Fast, deterministic\n";
    std::cout << "           - Best for: Regular fields, known dynamics\n";
    std::cout << "\n";
    std::cout << "  RK45     - Adaptive Runge-Kutta (4/5 order)\n";
    std::cout << "           - Variable timestep, error control\n";
    std::cout << "           - Best for: Stiff systems, Orbitrap, varying dynamics\n";
    std::cout << "\n";
    std::cout << "  Boris    - Boris pusher (magnetic fields)\n";
    std::cout << "           - Specialized for E×B fields\n";
    std::cout << "           - Best for: ICR, magnetic confinement\n";
    std::cout << "\nUsage:\n";
    std::cout << "  In config:   \"integrator\": \"RK4\"\n";
    std::cout << "  Via CLI:     --set simulation.integrator=RK45\n";
    std::cout << "\nNotes:\n";
    std::cout << "  - RK4 is default for most instrument types\n";
    std::cout << "  - Can be overridden per-domain in config\n";
}

void check_dependencies() {
    std::cout << "Checking Dependencies:\n";
    std::cout << "======================\n\n";
    
    // Core dependencies
    std::cout << "Core Libraries:\n";
    std::cout << "  HDF5         : " << H5_VERS_INFO << "\n";
    
#ifdef JSONCPP_VERSION_STRING
    std::cout << "  JsonCpp      : " << JSONCPP_VERSION_STRING << "\n";
#else
    std::cout << "  JsonCpp      : Found (version unknown)\n";
#endif
    
    std::cout << "  cxxopts      : 3.1.1\n";
    
    // Optional features
    std::cout << "\nOptional Features:\n";
#ifdef USE_CUDA
    std::cout << "  CUDA         : Enabled\n";
#else
    std::cout << "  CUDA         : Disabled\n";
#endif
    
#ifdef _OPENMP
    std::cout << "  OpenMP       : Enabled (version " << _OPENMP << ")\n";
#else
    std::cout << "  OpenMP       : Disabled\n";
#endif
    
    // Build configuration
    std::cout << "\nBuild Configuration:\n";
#ifdef __GNUC__
    std::cout << "  Compiler     : GCC " << __GNUC__ << "." << __GNUC_MINOR__ << "." << __GNUC_PATCHLEVEL__ << "\n";
#elif defined(__clang__)
    std::cout << "  Compiler     : Clang " << __clang_major__ << "." << __clang_minor__ << "." << __clang_patchlevel__ << "\n";
#elif defined(_MSC_VER)
    std::cout << "  Compiler     : MSVC " << _MSC_VER << "\n";
#else
    std::cout << "  Compiler     : Unknown\n";
#endif
    std::cout << "  C++ Standard : C++" << __cplusplus / 100 % 100 << "\n";
    
#ifdef NDEBUG
    std::cout << "  Build Type   : Release\n";
#else
    std::cout << "  Build Type   : Debug\n";
#endif
    
    std::cout << "\nStatus: ✓ All required dependencies found\n";
}

void validate_schema(const std::string& config_file) {
    std::cout << "Validating configuration against schema...\n";
    std::cout << "Config file: " << config_file << "\n\n";
    
    // Get absolute paths (validator.py needs them for URI resolution)
    char* cwd = getcwd(nullptr, 0);
    std::string abs_cwd = cwd ? std::string(cwd) : ".";
    free(cwd);
    
    std::string validator_path = abs_cwd + "/schema/validator.py";
    std::string schema_dir = abs_cwd + "/schema";
    
    // Check if validator exists
    std::ifstream validator_check(validator_path);
    
    if (!validator_check.good()) {
        std::cerr << "Error: Schema validator not found at " << validator_path << "\n";
        std::cerr << "Please ensure schema/validator.py exists\n";
        std::exit(1);
    }
    
    // Run validator
    std::string cmd = "python3 " + validator_path + " --schema-dir " + schema_dir + " " + config_file;
    std::cout << "Running: " << cmd << "\n";
    std::cout << "----------------------------------------\n";
    
    int result = std::system(cmd.c_str());
    
    if (result == 0) {
        std::cout << "\n✓ Configuration is valid\n";
    } else {
        std::cout << "\n✗ Configuration validation failed\n";
        std::exit(1);
    }
}

}  // namespace cli
}  // namespace ICARION

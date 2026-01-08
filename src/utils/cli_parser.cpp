// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "cli_parser.h"
#include <cxxopts.hpp>
#include <exception>
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
    (void)program_name;
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
        ("config", "Path to JSON configuration file", 
         cxxopts::value<std::string>(), "FILE")
        ("seed", "Override RNG seed from config (for reproducibility)", 
         cxxopts::value<unsigned int>(), "N")
        ("dry-run", "Validate configuration without running simulation")
        ("validate-config", "Validate config file and display warnings/errors")
        ("no-reactions", "Disable chemical reactions (collisions only)");
    
    // Configure positional arguments (allows: icarion config.json)
    parser.parse_positional({"config"});
    parser.positional_help("<config.json>");
    
    // === Logging options (Phase 1: ACTIVE) ===
    parser.add_options("Logging")
        ("log-level", "Set log level (DEBUG|INFO|WARN|ERROR) [default: INFO]", 
         cxxopts::value<std::string>()->default_value("INFO"), "LEVEL")
        ("log-file", "Write logs to FILE instead of console", 
         cxxopts::value<std::string>(), "FILE")
        ("log-format", "Log output format: text or json [default: text]",
         cxxopts::value<std::string>()->default_value("text"), "FORMAT")
        ("verbose", "Enable verbose output (alias for --log-level DEBUG)");
    
    // === Output options (Phase 1: ACTIVE) ===
    parser.add_options("Output")
        ("o,output", "Override output HDF5 filename", 
         cxxopts::value<std::string>(), "FILE")
        ("output-dir", "Override output directory", 
         cxxopts::value<std::string>(), "DIR")
        ("buffer-byte-cap", "Cap in-memory trajectory buffer (bytes, 0 = unlimited)", 
         cxxopts::value<uint64_t>(), "BYTES");
    
    // === Config overrides (Phase 1: ACTIVE) ===
    parser.add_options("Advanced")
        ("set", "Override config value (e.g., --set simulation.dt_s=1e-10)", 
         cxxopts::value<std::vector<std::string>>(), "KEY=VALUE");
    
    // === Information flags (Phase 1: ACTIVE) ===
    parser.add_options("Information")
        ("dump-build-info", "Show detailed build configuration and features")
        ("dump-hdf5-schema", "Show HDF5 output schema documentation")
        ("dump-config-schema", "Show path to JSON config schemas")
        ("list-collision-models", "List available collision models")
        ("list-integrators", "List available integrators")
        ("validate-schema", "Validate config against schema (uses validator.py)")
        ("check-deps", "Verify all dependencies and their versions");
    
    // === Performance/Profiling options ===
    parser.add_options("Performance")
        ("benchmark", "Print detailed timing statistics for simulation phases")
        ("profile", "Enable profiling instrumentation (for performance analysis)")
        ("profile-output", "Write profiling data to FILE (json or csv) [default: profile.json]",
         cxxopts::value<std::string>(), "FILE")
        ("threads", "Number of OpenMP threads for CPU parallelization (overrides OMP_NUM_THREADS)",
         cxxopts::value<int>(), "N");
    
    // Parse
    cxxopts::ParseResult result;
    try {
        result = parser.parse(argc, argv);
    } catch (const std::exception& e) {
        // Note: Logger not yet initialized at this point, use stderr
        std::cerr << "Error parsing arguments: " << e.what() << "\n";
        std::cerr << "Use --help for usage information\n";
        std::exit(1);
    }
    
    // === Handle --help and --version early ===
    if (result.count("help")) {
        std::cout << parser.help({"Core", "Logging", "Output", "Advanced", "Performance", "Information"}) << "\n";
        std::cout << "\nExamples:\n";
        std::cout << "  icarion config.json\n";
        std::cout << "  icarion --verbose config.json\n";
        std::cout << "  icarion --threads 8 config.json  # Use 8 CPU threads\n";
        std::cout << "  icarion --log-level DEBUG --log-file debug.log config.json\n";
        std::cout << "  icarion --output custom.h5 config.json\n";
        std::cout << "  icarion --set simulation.dt_s=1e-10 config.json\n";
        std::cout << "  icarion --dry-run --validate-config config.json\n";
        std::cout << "  icarion --benchmark config.json  # Show timing breakdown\n";
        std::cout << "  icarion --profile --profile-output profile.csv config.json\n";
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
    // Config file is required (either as --config FILE or positional)
    if (!result.count("config")) {
        // Note: Logger not yet initialized at this point, use stderr
        std::cerr << "Error: Missing configuration file\n";
        std::cerr << "Usage: icarion [OPTIONS] <config.json>\n";
        std::cerr << "Use --help for more information\n";
        std::exit(1);
    }
    
    opts.config_file = result["config"].as<std::string>();
    
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
            // Note: Logger not yet initialized at this point, use stderr
            std::cerr << "Error: Invalid log level '" << opts.log_level << "'\n";
            std::cerr << "Valid options: DEBUG, INFO, WARN, ERROR\n";
            std::exit(1);
        }
    }
    
    if (result.count("log-file")) {
        opts.log_file = result["log-file"].as<std::string>();
    }
    
    if (result.count("log-format")) {
        opts.log_format = result["log-format"].as<std::string>();
        
        // Validate log format
        if (opts.log_format != "text" && opts.log_format != "json") {
            // Note: Logger not yet initialized at this point, use stderr
            std::cerr << "Error: Invalid log format '" << opts.log_format << "'\n";
            std::cerr << "Valid options: text, json\n";
            std::exit(1);
        }
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

    if (result.count("buffer-byte-cap")) {
        opts.buffer_byte_cap = result["buffer-byte-cap"].as<uint64_t>();
    }
    
    // === Parse config overrides (Phase 1: ACTIVE) ===
    if (result.count("set")) {
        auto overrides = result["set"].as<std::vector<std::string>>();
        for (const auto& override : overrides) {
            size_t eq_pos = override.find('=');
            if (eq_pos == std::string::npos) {
                // Note: Logger not yet initialized at this point, use stderr
                std::cerr << "Error: Invalid --set format (expected KEY=VALUE): " << override << "\n";
                std::exit(1);
            }
            std::string key = override.substr(0, eq_pos);
            std::string value = override.substr(eq_pos + 1);
            opts.overrides[key] = value;
        }
    }
    
    // === Parse performance/profiling options ===
    opts.benchmark = result.count("benchmark") > 0;
    opts.profile = result.count("profile") > 0;
    
    if (result.count("profile-output")) {
        opts.profile_output = result["profile-output"].as<std::string>();
    } else if (opts.profile || opts.benchmark) {
        // Default output file if profiling is enabled
        opts.profile_output = "profile.json";
    }
    
    if (result.count("threads")) {
        int num_threads = result["threads"].as<int>();
        if (num_threads <= 0) {
            std::cerr << "Error: Invalid thread count " << num_threads << " (must be > 0)\n";
            std::exit(1);
        }
        opts.num_threads = num_threads;
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
ICARION HDF5 Output Schema v1.0 (FullConfig-based)
===================================================

Hierarchical Structure:
-----------------------

/metadata/
  /config/
    format_version      string      Config format version
    config_json         string      Full FullConfig as JSON (TODO(v1.1): serializer)
    dt_s                float64     Timestep [s]
    total_time_s        float64     Total simulation time [s]
    total_steps         int32       Number of integration steps
    write_interval      int32       Output write interval
    integrator          string      Solver type (RK4, RK45, Leapfrog, etc.)
    collision_model     string      Collision model (EHSS, HSS, etc.)
    enable_reactions        bool        Chemical reactions enabled?
    enable_space_charge     bool        Space charge enabled?
    enable_space_charge_gpu bool        Prefer GPU space charge (if available)?
    enable_gpu          bool        GPU acceleration enabled?
  
  /reproducibility/
    global_seed         uint32      RNG seed for reproducibility
    rng_algorithm       string      RNG type (std::mt19937_64)
    seed_scheme         string      Seeding strategy description
    git_hash            string      Git commit hash
    git_dirty           bool        Uncommitted changes in repo?
    code_version        string      ICARION version
    build_type          string      Release or Debug
    compiler_cxx        string      Compiler version
    build_info          string      Build information string
    openmp_threads      int32       OpenMP thread count (if enabled)
    cuda_version        string      CUDA version (if GPU enabled)
    /input_hash/
      config_sha256     string      SHA256 hash of config file
      species_db_sha256 string      SHA256 hash of species DB (N/A if embedded)
      reaction_db_sha256 string     SHA256 hash of reaction DB (N/A if embedded)
  
  /system/
    hostname            string      Machine hostname
    username            string      Username
    os                  string      Operating system + release
    kernel              string      Kernel version
    cpu_model           string      CPU model name
    cpu_cores           int32       Number of CPU cores
    memory_gb           float64     Total system memory [GB]
    gpu_model           string      GPU model (if CUDA enabled)
    gpu_memory_gb       float64     GPU memory [GB] (if CUDA enabled)
    driver_version      string      GPU driver version (if CUDA enabled)
    timestamp           string      ISO8601 simulation start time (UTC)
  
  /species/                         # Tabular format (pandas-compatible)
    names               string[N]   Species names (variable-length strings)
    mass_kg             float64[N]  Masses [kg]
    charge_C            float64[N]  Charges [C]
    mobility_m2Vs       float64[N]  Reduced mobilities [m²/(V·s)]
    ccs_m2              float64[N]  Collision cross sections [m²]
    
    Note: Only species referenced by ions are written (filters large databases)
  
  /reactions/                       # Tabular format (if reactions enabled)
    id                  string[R]   Reaction identifiers
    reactant_1          string[R]   First reactant species name
    reactant_2          string[R]   Second reactant (empty if unimolecular)
    product_1           string[R]   First product species name
    rate_constant   float64[R]  Rate constant [m³/s]
    type                int32[R]    Reaction type enum (2=two-body)
    
    Note: Only reactions involving present species are written (filters large networks)
  
  /completion/                      # Written at simulation end
    success             bool        Simulation completed successfully?
    final_time_s        float64     Final time reached [s]
    active_ions         int32       Number of active ions at completion
    completion_timestamp string     ISO8601 completion time (UTC)

/trajectory/                        # Time-series data (chunked, compressed)
  time                  float64[T]      Time snapshots [s]
  positions             float64[T×N×3]  Ion positions [m] (x,y,z)
  velocities            float64[T×N×3]  Ion velocities [m/s] (vx,vy,vz)
  species_ids           string[T×N]     Species name per ion per timestep
  domain_indices        int32[T×N]      Domain index per ion per timestep

/ions/                              # Per-ion initial conditions
  initial_species_id    string[N]       Initial species name
  initial_pos_x         float64[N]      Initial x position [m]
  initial_pos_y         float64[N]      Initial y position [m]
  initial_pos_z         float64[N]      Initial z position [m]
  initial_vel_x         float64[N]      Initial x velocity [m/s]
  initial_vel_y         float64[N]      Initial y velocity [m/s]
  initial_vel_z         float64[N]      Initial z velocity [m/s]
  birth_time_s          float64[N]      Ion birth time [s] (0 for initial ions)
  death_time_s          float64[N]      Ion death time [s] (-1 if still alive)
  charge_C              float64[N]      Ion charge [C]

/domains/                           # Per-domain configuration hierarchy
  /domain_0/
    name                string      Domain name
    instrument          string      Instrument type (IMS, TOF, LQIT, Orbitrap, etc.)
    solver              string      Solver type (RK4, RK45, Leapfrog, etc.)
    domain_index        int32       Domain index
    /geometry/
      length_m          float64     Domain length [m]
      radius_m          float64     Inner radius [m]
      radius_in_m       float64     Inner radius (for LQIT) [m]
      radius_out_m      float64     Outer radius (for LQIT) [m]
      origin_m          float64[3]  Origin position [m] (x,y,z)
    /environment/
      pressure_Pa       float64     Gas pressure [Pa]
      temperature_K     float64     Temperature [K]
      gas_species       string      Neutral gas species (N2, He, etc.)
      particle_density_m3 float64   Number density [m⁻³]
      mean_thermal_velocity_ms float64 Mean thermal velocity [m/s]
      gas_velocity_ms   float64[3]  Gas flow velocity [m/s] (x,y,z)
    /fields/
      /waveforms/                       # v1.1: Waveform library (if present)
        /<waveform_name>/
          type          string      Waveform type (constant, linear, quadratic, sinusoidal, pulsed, arbitrary)
          ...           varies      Type-specific parameters (start_value, end_value, frequency_Hz, etc.)
      /dc/
        axial_V         float64     DC drift voltage [V] (t=0 if waveform)
        EN_Td           float64     Reduced electric field [Td] (t=0 if waveform)
        quad_V          float64     Quadrupole DC voltage [V] (LQIT) (t=0 if waveform)
      /rf/
        voltage_V       float64     RF amplitude [V] (t=0 if waveform)
        frequency_Hz    float64     RF frequency [Hz] (t=0 if waveform)
        phase_rad       float64     RF phase [rad]
      /ac/
        voltage_V       float64     AC amplitude [V] (LQIT) (t=0 if waveform)
        frequency_Hz    float64     AC frequency [Hz] (t=0 if waveform)
  /domain_1/
    ...
  /domain_N/
    ...

Storage Properties:
-------------------
- Chunked storage: Trajectory datasets use [100 × N × 3] chunks for efficient appending
- GZIP compression: Level 6 applied to trajectory data (~3-5× size reduction)
- Unlimited dimensions: Time dimension is extendable (H5S_UNLIMITED)
- Variable-length strings: For species names, reaction IDs, and text fields
- SI units: All physical quantities use SI base units
- ISO8601 timestamps: UTC timezone (YYYY-MM-DDTHH:MM:SSZ format)
- Row-major order: C/Python convention (not Fortran)

Python Access Examples:
------------------------
See docs/HDF5_OUTPUT_STRUCTURE.md for complete Python examples with h5py and pandas.

Reproducibility Verification:
------------------------------
See docs/HDF5_OUTPUT_STRUCTURE.md for Python code to verify SHA256 hashes.

MATLAB Access:
--------------
MATLAB R2020a+ has native HDF5 support via h5read() and h5info().
See docs/HDF5_OUTPUT_STRUCTURE.md for MATLAB examples.

Full Documentation:
-------------------
- Detailed format specification: docs/HDF5_OUTPUT_STRUCTURE.md
- C++ API reference: src/core/io/hdf5Writer_v2.h
- Test examples and usage: tests/io/test_hdf5_writer_v2.cpp
- Migration plan: tmp/HDF5_REFACTORING_PLAN.md

Compatibility:
--------------
- HDF5 C library 1.10.0+
- Python h5py 2.10.0+, pandas 1.0.0+
- MATLAB R2020a+
- Julia HDF5.jl 0.15.0+

Command-Line Tools:
-------------------
View file structure:     h5ls -r icarion_output.h5
Dump structure:          h5dump -H icarion_output.h5
View specific dataset:   h5dump -d /trajectory/positions icarion_output.h5
Get compression info:    h5stat icarion_output.h5

Notes:
------
- SHA256 hashing is implemented for config files, species/reaction DBs are embedded
- All physical quantities use SI base units (meters, seconds, kilograms, etc.)
- Timestamps are in UTC timezone for international collaboration

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

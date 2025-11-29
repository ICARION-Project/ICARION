/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        main.cpp
 *   @brief       Entry point and orchestration layer of ICARION.
 *
 *   @details
 *   ICARION serves as the main execution driver for an ion trajectory simulation.
 *   It performs setup, data import, and solver execution in a defined pipeline:
 *
 *   1. Parse command-line arguments and apply overrides.
 *   2. Load configuration from JSON file (SSOT: FullConfig).
 *   3. Initialize ions from configuration.
 *   4. Create SimulationEngine with FullConfig (dependency injection).
 *   5. Run simulation via SimulationEngine::run().
 *   6. Report results and completion status.
 *
 *   The output consists of time-resolved trajectories, arrival time distributions,
 *   and optionally reaction histories (HDF5 format).
 *
 *   @date        2025-11-23
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include <chrono>
#include <exception>
#include <iostream>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

// Core systems (Phase 5: SSOT-compliant)
#include "core/config/loader/ConfigLoader.h"
#include "core/config/utils/ConfigOverride.h"
#include "core/integrator/SimulationEngine.h"
#include "core/log/Logger.h"
#include "core/utils/startupBanner.h"
#include "core/utils/Profiler.h"
#include "utils/cli_parser.h"
#include "main/setup/PhysicsSetup.h"

/**
 * @file main.cpp
 * @brief Entry point for the ICARION simulation.
 *
 * Initializes and executes a complete ICARION run using SimulationEngine.
 *
 * @param[in] argc Number of command-line arguments.
 * @param[in] argv Command-line arguments (expects JSON configuration path).
 * @return 0 on success, 1 on any runtime error.
 *
 * @throws std::runtime_error If configuration loading or simulation fails.
 *
 * @note The JSON file must define a complete FullConfig structure.
 *
 * @see config::ConfigLoader
 * @see integrator::SimulationEngine
 */

int main(int argc, char* argv[]) {
    using namespace ICARION;
    
    try {
        PROFILE_SCOPE("Total Execution");
        
        // === 1. Parse command-line arguments ===
        cli::CLIOptions opts = cli::parse_arguments(argc, argv);
        
        // Handle --help and --version (already printed by parser)
        if (opts.show_help || opts.show_version) {
            return 0;
        }
        
        // Enable profiler if requested
        if (opts.benchmark || opts.profile) {
            profiling::Profiler::getInstance().enable(true);
        }
        
        // === 2. Initialize logging ===
        {
            PROFILE_SCOPE("Logging Initialization");
            log::Logger::init(
                opts.log_level,
                opts.log_file.value_or(""),
                opts.log_format);
        }
        
        // Print startup banner (text format only)
        if (opts.log_format == "text") {
            utils::print_startup_banner(
                ICARION_VERSION,
                GIT_HASH,
                opts.config_file,
                opts.log_level,
                opts.log_file.value_or("")
            );
        } else {
            log::Logger::main()->info("ICARION v{} starting", ICARION_VERSION);
            log::Logger::main()->info("Git commit: {}", GIT_HASH);
            log::Logger::main()->info("Config file: {}", opts.config_file);
        }
        
        // === Handle information flags ===
        if (opts.dump_build_info) {
            cli::print_build_info();
            return 0;
        }
        
        if (opts.dump_hdf5_schema) {
            cli::print_hdf5_schema();
            return 0;
        }
        
        if (opts.dump_config_schema) {
            cli::print_config_schema();
            return 0;
        }
        
        if (opts.list_collision_models) {
            cli::list_collision_models();
            return 0;
        }
        
        if (opts.list_integrators) {
            cli::list_integrators();
            return 0;
        }
        
        if (opts.check_deps) {
            cli::check_dependencies();
            return 0;
        }
        
        if (opts.validate_schema) {
            cli::validate_schema(opts.config_file);
            return 0;
        }
    
    // === Apply logging options (Phase 1) ===
    // Note: Logger system handles file output via spdlog, no freopen needed
    // Old file redirection code removed (replaced by Logger::init with --log-file)

        // === Handle --validate-config ===
        if (opts.validate_config) {
            log::Logger::main()->info("=== ICARION Configuration Validation ===");
            log::Logger::main()->info("Config file: {}", opts.config_file);
            
            try {
                auto config = config::ConfigLoader::load(opts.config_file);
                auto validation = config.validate();
                
                log::Logger::main()->info("--- Validation Results ---");
                if (validation.valid && validation.warnings.empty()) {
                    log::Logger::main()->info("✓ Configuration is valid (no warnings)");
                    return 0;
                }
                
                if (validation.valid) {
                    log::Logger::main()->info("✓ Configuration is valid (with warnings)");
                } else {
                    log::Logger::main()->error("✗ Configuration has errors");
                }
                
                if (!validation.warnings.empty()) {
                    log::Logger::main()->warn("Warnings:");
                    for (const auto& warn : validation.warnings) {
                        log::Logger::main()->warn("  ⚠  {}", warn);
                    }
                }
                
                if (!validation.errors.empty()) {
                    log::Logger::main()->error("Errors:");
                    for (const auto& err : validation.errors) {
                        log::Logger::main()->error("  ✗  {}", err);
                    }
                }
                
                return validation.valid ? 0 : 1;
                
            } catch (const std::exception& e) {
                log::Logger::main()->error("✗ Configuration validation failed: {}", e.what());
                return 1;
            }
        }

        // === 3. Load configuration (SSOT: FullConfig) ===
        log::Logger::main()->info("Loading configuration: {}", opts.config_file);
        
        config::FullConfig config;
        {
            PROFILE_SCOPE("Config Loading");
            config = config::ConfigLoader::load(opts.config_file);
        }
        
        // === Apply CLI overrides ===
        if (!opts.overrides.empty()) {
            log::Logger::main()->info("Applying {} CLI overrides", opts.overrides.size());
            config::ConfigOverride::apply(config, opts.overrides);
            config.finalize_all();
        }
        
        if (opts.output_file.has_value()) {
            log::Logger::main()->info("Output file override: {}", opts.output_file.value());
            config.output.trajectory_file = opts.output_file.value();
        }
        
        if (opts.output_dir.has_value()) {
            log::Logger::main()->info("Output directory override: {}", opts.output_dir.value());
            config.output.folder = opts.output_dir.value();
        }
        
        if (opts.seed.has_value()) {
            log::Logger::main()->info("RNG seed override: {}", opts.seed.value());
            config.simulation.rng_seed = opts.seed.value();
        }
        
        if (opts.no_reactions) {
            log::Logger::main()->info("Disabling reactions (--no-reactions)");
            config.physics.enable_reactions = false;
        }
        
        // === Validate configuration ===
        auto validation = config.validate();
        if (!validation.valid) {
            log::Logger::main()->error("Configuration validation failed:");
            for (const auto& err : validation.errors) {
                log::Logger::main()->error("  - {}", err);
            }
            return 1;
        }
        
        if (!validation.warnings.empty()) {
            for (const auto& warn : validation.warnings) {
                log::Logger::main()->warn("  ⚠  {}", warn);
            }
        }
        
        // === Print simulation summary ===
        log::Logger::main()->info("");
        log::Logger::main()->info("=== Simulation Parameters ===");
        log::Logger::main()->info("Timestep:     {:.2e} s ({:.2f} ns)", 
                                  config.simulation.dt_s, config.simulation.dt_s * 1e9);
        log::Logger::main()->info("Total time:   {:.2e} s ({:.2f} µs)", 
                                  config.simulation.total_time_s, config.simulation.total_time_s * 1e6);
        log::Logger::main()->info("Domains:      {}", config.domains.size());
        log::Logger::main()->info("Species:      {}", config.species_db.size());
        log::Logger::main()->info("Reactions:    {}", config.reaction_db.size());
        log::Logger::main()->info("Output file:  {}", config.output.trajectory_file);
        log::Logger::main()->info("RNG seed:     {}", config.simulation.rng_seed);
        log::Logger::main()->info("=============================");
        log::Logger::main()->info("");
        
        // === Dry-run mode ===
        if (opts.dry_run) {
            log::Logger::main()->info("Dry-run mode: Configuration valid. Exiting.");
            return 0;
        }
        
        // === 4. Initialize ions ===
        size_t total_ion_count = 0;
        for (const auto& spec : config.ions.species) {
            total_ion_count += spec.count;
        }
        log::Logger::main()->info("Generating {} ions", total_ion_count);
        
        std::vector<IonState> ions;
        {
            PROFILE_SCOPE("Ion Generation");
            std::mt19937 rng(config.simulation.rng_seed);
            ions = config.generate_ions(rng);
        }
        log::Logger::main()->info("✓ {} ions generated", ions.size());
        
        // === 5. Create physics dependencies ===
        setup::PhysicsModules physics;
        {
            PROFILE_SCOPE("Physics Module Setup");
            physics = setup::PhysicsSetup::initialize(config, ions);
        }
        
        // === 6. Create SimulationEngine ===
        log::Logger::main()->info("Initializing SimulationEngine");
        
        integrator::SimulationEngine engine = [&]() {
            PROFILE_SCOPE("Engine Initialization");
            return integrator::SimulationEngine(
                config,
                physics.force_registries,
                physics.integrator,
                physics.collision_handler,
                physics.reaction_handler
            );
        }();
        
        // === 7. Run simulation ===
        log::Logger::main()->info("Starting simulation (t_max = {:.2e} s)", 
                                  config.simulation.total_time_s);
        
        auto start = std::chrono::high_resolution_clock::now();
        
        std::vector<IonState> final_ions;
        {
            PROFILE_SCOPE("Simulation Run");
            final_ions = engine.run(ions);
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_s = std::chrono::duration<double>(end - start).count();
        
        // === 8. Report results ===
        size_t active_count = std::count_if(
            final_ions.begin(), final_ions.end(),
            [](const auto& ion) { return ion.active; }
        );
        
        log::Logger::main()->info("");
        log::Logger::main()->info("=== Simulation Complete ===");
        log::Logger::main()->info("CPU time:     {:.3f} s", elapsed_s);
        log::Logger::main()->info("Active ions:  {}/{}", active_count, final_ions.size());
        log::Logger::main()->info("Output file:  {}", config.output.trajectory_file);
        log::Logger::main()->info("===========================");
        
        // Print profiling summary if enabled
        if (opts.benchmark || opts.profile) {
            profiling::Profiler::getInstance().printSummary();
            
            if (opts.profile_output.has_value()) {
                std::string filename = opts.profile_output.value();
                try {
                    if (filename.find(".csv") != std::string::npos) {
                        profiling::Profiler::getInstance().exportCSV(filename);
                    } else {
                        // Default to JSON
                        if (filename.find(".json") == std::string::npos) {
                            filename += ".json";
                        }
                        profiling::Profiler::getInstance().exportJSON(filename);
                    }
                    log::Logger::main()->info("Profile data written to: {}", filename);
                } catch (const std::exception& e) {
                    log::Logger::main()->error("Failed to write profile data: {}", e.what());
                }
            }
        }
        
    } catch (const std::exception& e) {
        log::Logger::main()->error("Fatal error: {}", e.what());
        log::Logger::shutdown();
        return 1;
    }
    
    log::Logger::shutdown();
    return 0;
}

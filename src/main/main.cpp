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
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/BorisStrategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/forces/SpaceChargeDirect.h"
#include "core/physics/forces/SpaceChargeGrid.h"
#include "core/physics/spacecharge/spaceChargeSolver.h"
#include "core/physics/collisions/CollisionHandlerFactory.h"
#include "core/physics/reactions/ReactionHandlerFactory.h"
#include "core/log/Logger.h"
#include "utils/cli_parser.h"
#include "core/utils/startupBanner.h"
#include "core/utils/Profiler.h"

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
        log::Logger::main()->info("Initializing physics modules");
        
        PROFILE_SCOPE("Physics Module Setup");
        
        // Create ForceRegistry for each domain (Phase 12 enhancement)
        std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries;
        for (const auto& domain : config.domains) {
            force_registries.push_back(std::make_shared<physics::ForceRegistry>(domain));
        }
        log::Logger::main()->info("Created {} ForceRegistry instances (one per domain)", 
                                  force_registries.size());
        
        // Auto-select space charge method based on ion count
        // N < 1000: Direct N-body (SpaceChargeForce, O(N²), exact)
        // N ≥ 1000: Grid-based Poisson solver (SpaceChargeSolver, O(N log N), fast)
        if (config.physics.enable_space_charge) {
            const size_t N = ions.size();
            constexpr size_t SPACE_CHARGE_THRESHOLD = 1000;
            
            if (N < SPACE_CHARGE_THRESHOLD) {
                // Use direct N-body Coulomb (exact, but O(N²))
                log::Logger::main()->info("Space charge: Using SpaceChargeDirect (N={} < {})", 
                                          N, SPACE_CHARGE_THRESHOLD);
                log::Logger::main()->info("  → Direct N-body Coulomb (exact, O(N²))");
                
                // Add SpaceChargeDirect to all domain registries
                constexpr double SOFTENING_LENGTH = 1e-10;  // 0.1 nm (prevents 1/r² divergence)
                for (auto& registry : force_registries) {
                    registry->add_force(std::make_unique<physics::SpaceChargeDirect>(SOFTENING_LENGTH));
                }
                log::Logger::main()->info("  ✓ SpaceChargeDirect added to {} registries (ε={:.2e} m)",
                                          force_registries.size(), SOFTENING_LENGTH);
            } else {
                // Use grid-based Poisson solver (fast, but approximate)
                log::Logger::main()->info("Space charge: Using SpaceChargeGrid (N={} >= {})", 
                                          N, SPACE_CHARGE_THRESHOLD);
                log::Logger::main()->info("  → Grid-based Poisson solver (fast, O(N log N))");
                
                // Estimate domain size from ion initial positions
                Vec3 min_pos = ions[0].pos;
                Vec3 max_pos = ions[0].pos;
                for (const auto& ion : ions) {
                    min_pos.x = std::min(min_pos.x, ion.pos.x);
                    min_pos.y = std::min(min_pos.y, ion.pos.y);
                    min_pos.z = std::min(min_pos.z, ion.pos.z);
                    max_pos.x = std::max(max_pos.x, ion.pos.x);
                    max_pos.y = std::max(max_pos.y, ion.pos.y);
                    max_pos.z = std::max(max_pos.z, ion.pos.z);
                }
                
                Vec3 domain_size = {max_pos.x - min_pos.x, max_pos.y - min_pos.y, max_pos.z - min_pos.z};
                Vec3 domain_center = {(min_pos.x + max_pos.x) / 2, (min_pos.y + max_pos.y) / 2, (min_pos.z + max_pos.z) / 2};
                
                // Add 50% margin to domain size (ions will move)
                domain_size = domain_size * 1.5;
                
                // Grid resolution: Aim for ~1mm cells (adjust based on domain)
                constexpr int TARGET_GRID_SIZE = 64;  // 64³ = 262k cells (good balance)
                double cell_size_x = domain_size.x / TARGET_GRID_SIZE;
                double cell_size_y = domain_size.y / TARGET_GRID_SIZE;
                double cell_size_z = domain_size.z / TARGET_GRID_SIZE;
                
                // Use uniform cell size (max of xyz)
                double cell_size = std::max({cell_size_x, cell_size_y, cell_size_z, 1e-4});  // Min 0.1mm
                
                Vec3 grid_origin = {
                    domain_center.x - (TARGET_GRID_SIZE * cell_size) / 2,
                    domain_center.y - (TARGET_GRID_SIZE * cell_size) / 2,
                    domain_center.z - (TARGET_GRID_SIZE * cell_size) / 2
                };
                
                log::Logger::main()->info("  Grid: {}³ cells, {:.2e} m cell size", TARGET_GRID_SIZE, cell_size);
                log::Logger::main()->info("  Domain: [{:.3f}, {:.3f}] x [{:.3f}, {:.3f}] x [{:.3f}, {:.3f}] mm",
                                          grid_origin.x * 1e3, (grid_origin.x + TARGET_GRID_SIZE * cell_size) * 1e3,
                                          grid_origin.y * 1e3, (grid_origin.y + TARGET_GRID_SIZE * cell_size) * 1e3,
                                          grid_origin.z * 1e3, (grid_origin.z + TARGET_GRID_SIZE * cell_size) * 1e3);
                
                // Create solver
                auto sc_solver = std::make_shared<SpaceChargeSolver>(
                    TARGET_GRID_SIZE, TARGET_GRID_SIZE, TARGET_GRID_SIZE,
                    cell_size, cell_size, cell_size,
                    grid_origin
                );
                
                // Wrap solver in IForce interface and add to registries
                for (auto& registry : force_registries) {
                    registry->add_force(std::make_unique<physics::SpaceChargeGrid>(sc_solver));
                }
                log::Logger::main()->info("  ✓ SpaceChargeGrid added to {} registries",
                                          force_registries.size());
            }
        }
        
        // Create integration strategy (from config.simulation.integrator)
        std::shared_ptr<integrator::IIntegrationStrategy> integration_strategy;
        if (config.simulation.integrator == "RK4" || config.simulation.integrator == "rk4") {
            integration_strategy = std::make_shared<integrator::RK4Strategy>();
            log::Logger::main()->info("Using RK4 integrator");
        } else if (config.simulation.integrator == "RK45" || config.simulation.integrator == "rk45") {
            integration_strategy = std::make_shared<integrator::RK45Strategy>();
            log::Logger::main()->info("Using RK45 integrator");
        } else if (config.simulation.integrator == "Boris" || config.simulation.integrator == "boris") {
            integration_strategy = std::make_shared<integrator::BorisStrategy>();
            log::Logger::main()->info("Using Boris integrator");
        } else {
            // Default fallback
            log::Logger::main()->warn("Unknown integrator '{}', defaulting to RK45", 
                                      config.simulation.integrator);
            integration_strategy = std::make_shared<integrator::RK45Strategy>();
        }
        
        // Create collision handler (from config.physics.collision_model)
        const double gamma_for_ou = 0.0;  // OU damping coefficient not used for stochastic models
        std::shared_ptr<physics::ICollisionHandler> collision_handler = 
            physics::CollisionHandlerFactory::create(
                config.physics,
                nullptr,        // geometry map (EHSS only; provided in SimulationEngine paths)
                gamma_for_ou,
                false,          // enable_logging
                &config.species_db
            );
        
        // Create reaction handler (from config.physics.enable_reactions)
        std::shared_ptr<physics::IReactionHandler> reaction_handler = 
            physics::ReactionHandlerFactory::create(config.physics);
        
        // === 6. Create SimulationEngine ===
        log::Logger::main()->info("Initializing SimulationEngine");
        
        integrator::SimulationEngine engine = [&]() {
            PROFILE_SCOPE("Engine Initialization");
            return integrator::SimulationEngine(
                config,
                force_registries,  // Vector of registries (one per domain)
                integration_strategy,
                collision_handler,
                reaction_handler
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

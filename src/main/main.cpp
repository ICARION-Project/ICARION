/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        ICARION.cpp
 *   @brief       Entry point and orchestration layer of ICARION.
 *
 *   @details
 *   ICARION serves as the main execution driver for an ion trajectory simulation.
 *   It performs setup, data import, and solver executation in a defined pipeline:
 *
 *   1. Parse global parameters from the provided JSON configuration file.
 *   2. Load instrument domains and their geometry, environment and field setting.
 *   3. Initialize species and reaction database (if enabled).
 *   4. Initialize the ion ensemble from the cloud definition.
 *   5. Execute the time integration using the selected solver (RK4 or RK45), 
 *      depending on the instrument definition.
 *   6. Export results to HDF5 format.
 *
 *   The output consists of time-resolved trajectories, arrival time distributions,
 *   and optionalky reaction histories.
 *
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include <chrono>
#include <exception>
#include <iostream>
#include <fstream>
#include <sstream>

#include "integrator/integrator.h"
#include "core/io/hdf5Writer.h"
#include "core/param/paramUtils.h"
#include "utils/cli_parser.h"  // CLI argument parser
#include "core/log/Logger.h"  // Structured logging
#include "core/utils/simulationsUtils.h"  // print_domain_summary, print_results
#include "core/physics/reactions/reactionUtils.h"  // ReactionEntry struct only
#include "core/io/fieldArrayLoader.h"
#include "core/io/speciesLoader.h"  // io::SpeciesDatabase (temporary until Phase 5)

// New config system (Phase 4)
#include "core/config/loader/ConfigLoader.h"
#include "core/config/adapter/LegacyAdapter.h"
#include "core/config/utils/ConfigOverride.h"

/**
 * @file ICARION.cpp
 * @brief Entry point for the ICARION simulation.
 *
 * Initializes and executes a complete ICARION run.
 *
 * @param[in] argc Number of command-line arguments.
 * @param[in] argv Command-line arguments (expects JSON configuration path).
 * @return 0 on success, 1 on any runtime error.
 *
 * @throws std::runtime_error If parameter parsing, data loading, or solver execution fails.
 *
 * @note The JSON file must define:
 *   - Global parameters (`global` section)
 *   - One or more instrument domains (`domains` array)
 *   - Optionally, species and reaction files.
 *
 * @see load_global_params()
 * @see integrate_trajectory()
 * @see load_single_domain()
 * @see load_field_array()
 * @see load_speciesDB()
 * @see load_reactions()
 */

int main(int argc, char* argv[]) {
    // Parse command-line arguments
    ICARION::cli::CLIOptions opts = ICARION::cli::parse_arguments(argc, argv);
    
    // Handle --help and --version (already printed by parser)
    if (opts.show_help || opts.show_version) {
        return 0;
    }
    
    // === Initialize logging FIRST (before any other operations) ===
    bool json_logs = false;  // TODO: Add --log-format flag in future
    ICARION::log::Logger::init(
        opts.log_level,
        opts.log_file.value_or(""),
        json_logs
    );
    
    // Log startup information
    ICARION::log::Logger::main()->info("ICARION v{} starting", ICARION_VERSION);
    ICARION::log::Logger::main()->info("Git commit: {}", GIT_HASH);
    ICARION::log::Logger::main()->info("Config file: {}", opts.config_file);
    
    // === Handle information flags ===
    if (opts.dump_build_info) {
        ICARION::cli::print_build_info();
        return 0;
    }
    
    if (opts.dump_hdf5_schema) {
        ICARION::cli::print_hdf5_schema();
        return 0;
    }
    
    if (opts.dump_config_schema) {
        ICARION::cli::print_config_schema();
        return 0;
    }
    
    if (opts.list_collision_models) {
        ICARION::cli::list_collision_models();
        return 0;
    }
    
    if (opts.list_integrators) {
        ICARION::cli::list_integrators();
        return 0;
    }
    
    if (opts.check_deps) {
        ICARION::cli::check_dependencies();
        return 0;
    }
    
    if (opts.validate_schema) {
        ICARION::cli::validate_schema(opts.config_file);
        return 0;
    }
    
    // === Apply logging options (Phase 1) ===
    // Note: Logger system handles file output via spdlog, no freopen needed
    // Old file redirection code removed (replaced by Logger::init with --log-file)

    // === Handle --validate-config (Phase 4) ===
    if (opts.validate_config) {
        using ICARION::log::Logger;
        
        Logger::main()->info("=== ICARION Configuration Validation ===");
        Logger::main()->info("Config file: {}", opts.config_file);
        
        try {
            // Load config using new system
            auto config = ICARION::config::ConfigLoader::load(opts.config_file);
            
            // Validate (already done in load(), but we want to show results)
            auto validation = config.validate();
            
            Logger::main()->info("--- Validation Results ---");
            if (validation.valid && validation.warnings.empty()) {
                Logger::main()->info("✓ Configuration is valid (no warnings)");
                return 0;
            }
            
            if (validation.valid) {
                Logger::main()->info("✓ Configuration is valid (with warnings)");
            } else {
                Logger::main()->error("✗ Configuration has errors");
            }
            
            if (!validation.warnings.empty()) {
                Logger::main()->warn("Warnings:");
                for (const auto& warn : validation.warnings) {
                    Logger::main()->warn("  ⚠  {}", warn);
                }
            }
            
            if (!validation.errors.empty()) {
                Logger::main()->error("Errors:");
                for (const auto& err : validation.errors) {
                    Logger::main()->error("  ✗  {}", err);
                }
            }
            
            return validation.valid ? 0 : 1;
            
        } catch (const std::exception& e) {
            Logger::main()->error("✗ Configuration validation failed: {}", e.what());
            return 1;
        }
    }

    try {
        // === 1. Load configuration (new system with legacy adapter) ===
        // TODO: Remove LegacyAdapter once integrator/physics/IO are refactored
        auto full_config = ICARION::config::ConfigLoader::load(opts.config_file);
        
        // === Apply CLI overrides (Phase 1: --set support) ===
        if (!opts.overrides.empty()) {
            ICARION::config::ConfigOverride::apply(full_config, opts.overrides);
            // Recompute derived values after overrides
            full_config.finalize_all();
        }
        
        // === Apply output overrides (Phase 1) ===
        if (opts.output_file.has_value()) {
            std::cout << "[CLI] Output file override: " << opts.output_file.value() << "\n";
            full_config.output.trajectory_file = opts.output_file.value();
        }
        
        if (opts.output_dir.has_value()) {
            std::cout << "[CLI] Output directory override: " << opts.output_dir.value() << "\n";
            full_config.output.folder = opts.output_dir.value();
        }
        
        // Convert to legacy GlobalParams
        GlobalParams gParams = ICARION::config::LegacyAdapter::to_global_params(full_config);
        
        // Override RNG seed if --seed was provided
        if (opts.seed.has_value()) {
            std::cout << "[CLI] Overriding RNG seed: " << opts.seed.value() << "\n";
            gParams.rng_seed = opts.seed.value();
        }
        
        // Override reaction flag if --no-reactions was provided
        if (opts.no_reactions) {
            std::cout << "[CLI] Disabling reactions (--no-reactions)\n";
            gParams.enable_reactions = false;
        }
        
        // Read full configuration JSON string for metadata/logging
        std::string config_json;
        std::ifstream config_file(opts.config_file);
        if (config_file.is_open()) {
            std::stringstream buffer;
            buffer << config_file.rdbuf();
            config_json = buffer.str();
            config_file.close();
        }

        // === 3. Load instrument domains (via legacy adapter) ===
        std::vector<InstrumentDomain> domains = ICARION::config::LegacyAdapter::to_instrument_domains(full_config);
        
        // Load field arrays (if specified)
        for (auto& dom : domains) {
            if (!dom.FA_file.empty()) {
                dom.fieldArray = load_field_array(dom.FA_file);
                dom.fieldArrayLoaded = dom.fieldArray.is_valid();
            }
            if (!dom.FA_terms.empty()) {
                size_t ok = 0;
                for (auto& t : dom.FA_terms) {
                    try {
                        t.field = load_field_array(t.file);
                        t.loaded = t.field.is_valid();
                        if (t.loaded) ok++;
                    } catch (...) {
                        t.loaded = false;
                    }
                }
                dom.fieldArrayLoaded = (ok == dom.FA_terms.size() && ok > 0) || dom.fieldArrayLoaded;
            }
        }
        ICARION::log::Logger::config()->info("Loaded {} instrument domains", domains.size());
        ICARION::utils::print_domain_summary(domains);

        // === Save input config for reproducibility ===
        // Log the input configuration to output folder so simulation can be exactly reproduced
        std::string config_log_path = gParams.output_file;
        // Replace .h5 extension with _config.json
        size_t ext_pos = config_log_path.find_last_of('.');
        if (ext_pos != std::string::npos) {
            config_log_path = config_log_path.substr(0, ext_pos) + "_config.json";
        } else {
            config_log_path += "_config.json";
        }
        
        std::ofstream config_log(config_log_path);
        if (config_log.is_open()) {
            config_log << config_json;
            config_log.close();
            ICARION::log::Logger::main()->info("✓ Input configuration saved to: {}", config_log_path);
        } else {
            ICARION::log::Logger::main()->warn("Could not save input configuration to {}", config_log_path);
        }
        
        // === Print simulation parameters summary ===
        using ICARION::log::Logger;
        
        Logger::main()->info("");
        Logger::main()->info("=== Simulation Parameters Summary ===");
        Logger::main()->info("Timestep:        {} ns", gParams.dt_s * 1e9);
        Logger::main()->info("Total steps:     {}", gParams.sim_time_steps);
        Logger::main()->info("Max time:        {} µs", gParams.dt_s * gParams.sim_time_steps * 1e6);
        Logger::main()->info("Write interval:  {}", gParams.write_interval);
        
        std::string collision_model;
        switch (gParams.collisionModel) {
            case ICARION::core::CollisionModel::NoCollisions: collision_model = "NoCollisions"; break;
            case ICARION::core::CollisionModel::HardSphere: collision_model = "HardSphere"; break;
            case ICARION::core::CollisionModel::Langevin: collision_model = "Langevin"; break;
            case ICARION::core::CollisionModel::Friction: collision_model = "Friction"; break;
            case ICARION::core::CollisionModel::EHSS: collision_model = "EHSS"; break;
            case ICARION::core::CollisionModel::HSMC: collision_model = "HSMC"; break;
            default: collision_model = "Unknown"; break;
        }
        Logger::main()->info("Collision model: {}", collision_model);
        Logger::main()->info("Reactions:       {}", gParams.enable_reactions ? "enabled" : "disabled");
        Logger::main()->info("Space charge:    {}", gParams.enable_space_charge ? "enabled" : "disabled");
        Logger::main()->info("GPU:             {}", gParams.enable_gpu ? "enabled" : "disabled");
        Logger::main()->info("OpenMP:          {}", gParams.parallelization ? "enabled" : "disabled");
        Logger::main()->info("RNG seed:        {}", gParams.rng_seed);
        Logger::main()->info("Output file:     {}", gParams.output_file);
        Logger::main()->info("=====================================");
        Logger::main()->info("");

        run_guard_check_global(gParams);

        // === 4. Load physical models ===
        // Species and reactions are ALREADY loaded in full_config by ConfigLoader!
        // (ConfigLoader calls full_config.load_databases() automatically)
        Logger::main()->info("");
        Logger::main()->info("=== Physical Models ===");
        Logger::config()->info("Species loaded: {}", full_config.species_db.size());
        Logger::config()->info("Reactions loaded: {}", full_config.reaction_db.size());
        
        // --- Convert species for integrator (temporary until Phase 5) ---
        ICARION::io::SpeciesDatabase speciesDB;
        
        for (const auto& [id, props] : full_config.species_db.species) {
            ICARION::io::Species species;
            species.id = id;
            species.name = props.name.value_or(id);
            species.mass_u = props.mass_amu;
            species.mass_kg = props.mass_kg;
            species.charge_e = props.charge;
            species.charge_C = props.charge_C;
            species.mobility_m2Vs = props.mobility_m2Vs;
            species.CCS_m2 = props.CCS_m2;
            species.geometry_file = props.geometry_file;
            
            // Calculate reduced mass (if domain available)
            if (!domains.empty()) {
                double neutral_mass_kg = domains[0].env.neutral_mass_kg;
                species.reduced_mass_kg = (species.mass_kg * neutral_mass_kg) / 
                                         (species.mass_kg + neutral_mass_kg);
            } else {
                species.reduced_mass_kg = species.mass_kg;  // Fallback
            }
            
            speciesDB.add(species);
        }
        
        // --- Convert reactions for integrator (inline, no adapter class) ---
        std::vector<ReactionEntry> reaction_list;
        
        if (gParams.enable_reactions && full_config.reaction_db.size() > 0) {
            reaction_list.reserve(full_config.reaction_db.reactions.size());
            
            size_t skipped = 0;
            for (const auto& rxn : full_config.reaction_db.reactions) {
                // Validate: single-reactant → single-product only
                if (rxn.reactant.empty() || rxn.product.empty()) {
                    Logger::config()->warn("Skipping reaction '{}' - missing reactant or product", rxn.id);
                    skipped++;
                    continue;
                }
                
                // Create legacy entry
                ReactionEntry entry;
                entry.reactant = rxn.reactant;
                entry.product = rxn.product;
                entry.rate_constant = rxn.rate_constant_m3s;
                entry.neutral_concentration = 0.0;  // Computed dynamically in integrator
                
                // Convert order terms
                for (const auto& term : rxn.order_terms) {
                    ReactionOrderTerm legacy_term;
                    legacy_term.species = term.species;
                    legacy_term.exponent = term.exponent;
                    entry.order.push_back(legacy_term);
                }
                
                reaction_list.push_back(entry);
            }
            
            if (skipped > 0) {
                Logger::config()->info("✓ {} reactions converted for integrator ({} skipped)", reaction_list.size(), skipped);
            } else {
                Logger::config()->info("✓ {} reactions converted for integrator", reaction_list.size());
            }
            
        } else if (!gParams.enable_reactions) {
            Logger::config()->info("ℹ Reactions disabled (enable_reactions=false)");
        } else {
            ICARION::log::Logger::config()->info("No reactions loaded");
        }

        // === Dry-run mode: Validate configuration and exit ===
        if (opts.dry_run) {
            Logger::main()->info("");
            Logger::main()->info("=== Dry-run mode: Configuration validation ===");
            Logger::config()->info("JSON configuration loaded successfully");
            Logger::config()->info("Species database loaded: {} species", speciesDB.size());
            Logger::config()->info("Reactions loaded: {} reactions", reaction_list.size());
            Logger::main()->info("✓ Domains configured: {} domain(s)", domains.size());
            Logger::main()->info("✓ RNG seed: {}", gParams.rng_seed);
            Logger::main()->info("✓ Output file: {}", gParams.output_file);
            Logger::main()->info("");
            Logger::main()->info("Configuration valid. Exiting without running simulation.");
            return 0;
        }

        // === 5. Initialize ions (or load from checkpoint) ===
        std::vector<IonState> ions;
        double t_start = gParams.t_eval.front();
        double t_end = gParams.t_eval.back();
        
        if (!gParams.continue_from.empty()) {
            std::cout << "Continue mode: Loading state from " << gParams.continue_from << "\n";
            
            // Load ions from HDF5 checkpoint
            ions = ICARION::io::load_final_state_from_HDF5(gParams.continue_from);
            
            // Read final time from HDF5 metadata
            try {
                H5::H5File file(gParams.continue_from, H5F_ACC_RDONLY);
                H5::Attribute attr = file.openAttribute("final_time_s");
                double checkpoint_time;
                attr.read(H5::PredType::NATIVE_DOUBLE, &checkpoint_time);
                file.close();
                
                t_start = checkpoint_time;
                t_end = checkpoint_time + gParams.continue_time_s;
                
                std::cout << "  Checkpoint time: " << checkpoint_time << " s\n";
                std::cout << "  Continuing for:  " << gParams.continue_time_s << " s\n";
                std::cout << "  New end time:    " << t_end << " s\n";
                
            } catch (const std::exception& e) {
                throw std::runtime_error("Failed to read checkpoint metadata: " + std::string(e.what()));
            }
            
        } else {
            // Normal initialization from ion configuration
            std::mt19937 rng(gParams.rng_seed);
            ions = full_config.generate_ions(rng);
            ICARION::log::Logger::main()->info("Generated {} ions from configuration", ions.size());
        }

        // === 6. Run simulation ===
        LOG_TIMER("simulation_total");
        
        ICARION::log::Logger::main()->info("Starting integration ({} ions, {:.6f} s total time, dt={:.2e} s)",
                                            ions.size(), t_end - t_start, gParams.dt_s);
        
        auto start = std::chrono::high_resolution_clock::now();
        std::vector<IonState> final_ions = integrate_trajectory(
            ions, 
            t_start, 
            t_end, 
            gParams.dt_s, 
            gParams,
            speciesDB, 
            reaction_list, 
            domains,
            RK45Settings(),
            nullptr  // No more RunLogger
        );
        
        auto             end    = std::chrono::high_resolution_clock::now();
        
        // === 8. Write completion metadata ===
        int active_count = 0;
        for (const auto& ion : final_ions) {
            if (ion.active) active_count++;
        }
        
        try {
            H5::H5File hdf5_file(gParams.output_file + ".h5", H5F_ACC_RDWR);
            bool simulation_complete = (active_count == 0);  // All ions exited/deactivated
            
            // Get git hash at compile time if available
            std::string git_hash;
            #ifdef GIT_COMMIT_HASH
                git_hash = GIT_COMMIT_HASH;
            #endif
            
            ICARION::io::write_simulation_metadata(hdf5_file, simulation_complete, active_count, t_end, 
                                     gParams.rng_seed, git_hash, config_json);
            hdf5_file.close();
            
            if (!simulation_complete && gParams.auto_continue_if_active) {
                ICARION::log::Logger::main()->warn("Simulation incomplete: {} ions still active at t={:.6f} s",
                                                    active_count, t_end);
                ICARION::log::Logger::main()->info("Consider using continue mode: \"continue_from\": \"{}.h5\", \"continue_time_s\": {}",
                                                    gParams.output_file, gParams.continue_time_s);
            }
        } catch (const std::exception& e) {
            ICARION::log::Logger::hdf5()->warn("Failed to write metadata: {}", e.what());
        }
        
        // === 10. Optional result printout ===
        if (gParams.print_results) {
            ICARION::utils::print_results(final_ions, 100);
        }
        
        ICARION::log::Logger::hdf5()->info("HDF5 file written successfully");
        
        double elapsed_s = std::chrono::duration<double>(end - start).count();
        ICARION::log::Logger::main()->info("Simulation completed in {:.3f} s CPU time", elapsed_s);
        ICARION::log::Logger::main()->info("Output file: {}.h5", gParams.output_file);
    } catch (const std::exception& e) {
        ICARION::log::Logger::main()->error("Fatal error: {}", e.what());
        ICARION::log::Logger::shutdown();
        return 1;
    }

    ICARION::log::Logger::shutdown();
    return 0;
}

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
#include "core/io/logger.h"
#include "core/io/hdf5Writer.h"
#include "core/param/paramUtils.h"
#include "utils/cli_parser.h"  // CLI argument parser
#include "core/utils/simulationsUtils.h"  // print_domain_summary, print_results
#include "core/physics/reactions/reactionUtils.h"  // ReactionEntry struct only
#include "core/io/fieldArrayLoader.h"
#include "core/io/speciesLoader.h"  // io::SpeciesDatabase (temporary until Phase 5)

// New config system (Phase 4)
#include "core/config/loader/ConfigLoader.h"
#include "core/config/adapter/LegacyAdapter.h"

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
    
    // Handle --help and --version
    if (opts.show_help) {
        ICARION::cli::print_help(argv[0]);
        return 0;
    }
    if (opts.show_version) {
        ICARION::cli::print_version();
        return 0;
    }

    // === Handle --validate-config (Phase 4) ===
    if (opts.validate_config) {
        std::cout << "=== ICARION Configuration Validation ===" << std::endl;
        std::cout << "Config file: " << opts.config_file << std::endl << std::endl;
        
        try {
            // Load config using new system
            auto config = ICARION::config::ConfigLoader::load(opts.config_file);
            
            // Validate (already done in load(), but we want to show results)
            auto validation = config.validate();
            
            std::cout << "--- Validation Results ---" << std::endl;
            if (validation.valid && validation.warnings.empty()) {
                std::cout << "✓ Configuration is valid (no warnings)" << std::endl;
                return 0;
            }
            
            if (validation.valid) {
                std::cout << "✓ Configuration is valid (with warnings)" << std::endl;
            } else {
                std::cout << "✗ Configuration has errors" << std::endl;
            }
            
            if (!validation.warnings.empty()) {
                std::cout << std::endl << "Warnings:" << std::endl;
                for (const auto& warn : validation.warnings) {
                    std::cout << "  ⚠  " << warn << std::endl;
                }
            }
            
            if (!validation.errors.empty()) {
                std::cout << std::endl << "Errors:" << std::endl;
                for (const auto& err : validation.errors) {
                    std::cout << "  ✗  " << err << std::endl;
                }
            }
            
            return validation.valid ? 0 : 1;
            
        } catch (const std::exception& e) {
            std::cerr << "✗ Configuration validation failed: " << e.what() << std::endl;
            return 1;
        }
    }

    try {
        // === 1. Load configuration (new system with legacy adapter) ===
        // TODO: Remove LegacyAdapter once integrator/physics/IO are refactored
        auto full_config = ICARION::config::ConfigLoader::load(opts.config_file);
        GlobalParams gParams = ICARION::config::LegacyAdapter::to_global_params(full_config);
        
        // Override RNG seed if --seed was provided
        if (opts.seed.has_value()) {
            std::cout << "Overriding RNG seed from command line: " << opts.seed.value() << "\n";
            gParams.rng_seed = opts.seed.value();
        }
        
        // Override reaction flag if --no-reactions was provided
        if (opts.no_reactions) {
            std::cout << "Disabling reactions (--no-reactions)\n";
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
        std::cout << "Loaded " << domains.size() << " instrument domains.\n";
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
            std::cout << "✓ Input configuration saved to: " << config_log_path << "\n";
        } else {
            std::cerr << "Warning: Could not save input configuration to " << config_log_path << "\n";
        }
        
        // === Print simulation parameters summary ===
        std::cout << "\n=== Simulation Parameters Summary ===\n";
        std::cout << "Timestep:        " << gParams.dt_s * 1e9 << " ns\n";
        std::cout << "Total steps:     " << gParams.sim_time_steps << "\n";
        std::cout << "Max time:        " << (gParams.dt_s * gParams.sim_time_steps * 1e6) << " µs\n";
        std::cout << "Write interval:  " << gParams.write_interval << "\n";
        std::cout << "Collision model: ";
        switch (gParams.collisionModel) {
            case ICARION::core::CollisionModel::NoCollisions: std::cout << "NoCollisions"; break;
            case ICARION::core::CollisionModel::HardSphere: std::cout << "HardSphere"; break;
            case ICARION::core::CollisionModel::Langevin: std::cout << "Langevin"; break;
            case ICARION::core::CollisionModel::Friction: std::cout << "Friction"; break;
            case ICARION::core::CollisionModel::EHSS: std::cout << "EHSS"; break;
            case ICARION::core::CollisionModel::HSMC: std::cout << "HSMC"; break;
            default: std::cout << "Unknown"; break;
        }
        std::cout << "\n";
        std::cout << "Reactions:       " << (gParams.enable_reactions ? "enabled" : "disabled") << "\n";
        std::cout << "Space charge:    " << (gParams.enable_space_charge ? "enabled" : "disabled") << "\n";
        std::cout << "GPU:             " << (gParams.enable_gpu ? "enabled" : "disabled") << "\n";
        std::cout << "OpenMP:          " << (gParams.parallelization ? "enabled" : "disabled") << "\n";
        std::cout << "RNG seed:        " << gParams.rng_seed << "\n";
        std::cout << "Output file:     " << gParams.output_file << "\n";
        std::cout << "=====================================\n\n";

        run_guard_check_global(gParams);

        // === 4. Load physical models ===
        // Species and reactions are ALREADY loaded in full_config by ConfigLoader!
        // (ConfigLoader calls full_config.load_databases() automatically)
        std::cout << "\n=== Physical Models ===\n";
        std::cout << "Species loaded:  " << full_config.species_db.size() << "\n";
        std::cout << "Reactions loaded: " << full_config.reaction_db.size() << "\n";
        
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
                    std::cerr << "Warning: Skipping reaction '" << rxn.id 
                              << "' - missing reactant or product\n";
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
            
            std::cout << "✓ " << reaction_list.size() << " reactions converted for integrator";
            if (skipped > 0) {
                std::cout << " (" << skipped << " skipped)";
            }
            std::cout << "\n";
            
        } else if (!gParams.enable_reactions) {
            std::cout << "ℹ  Reactions disabled (enable_reactions=false)\n";
        } else {
            std::cout << "ℹ  No reactions loaded\n";
        }

        // === Dry-run mode: Validate configuration and exit ===
        if (opts.dry_run) {
            std::cout << "\n=== Dry-run mode: Configuration validation ===\n";
            std::cout << "✓ JSON configuration loaded successfully\n";
            std::cout << "✓ Species database loaded: " << speciesDB.size() << " species\n";
            std::cout << "✓ Reactions loaded: " << reaction_list.size() << " reactions\n";
            std::cout << "✓ Domains configured: " << domains.size() << " domain(s)\n";
            std::cout << "✓ RNG seed: " << gParams.rng_seed << "\n";
            std::cout << "✓ Output file: " << gParams.output_file << "\n";
            std::cout << "\nConfiguration valid. Exiting without running simulation.\n";
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
            std::cout << "Generated " << ions.size() << " ions from configuration.\n";
        }

        // === 6. Initialize logger ===
        ICARION::io::RunLogger logger(gParams.output_file);
        logger.writeHeader();
        logger.writeGlobalParams(gParams, ions, speciesDB, reaction_list);
        logger.writeInstrumentDomains(domains);

        // === 7. Run simulation ===
        auto             start  = std::chrono::high_resolution_clock::now();
        logger.log("Starting integration...");
        SimulationResult result = integrate_trajectory(
            ions, 
            t_start, 
            t_end, 
            gParams.dt_s, 
            gParams,
            speciesDB, 
            reaction_list, 
            domains,
            RK45Settings(),
            &logger
        );
        
        auto             end    = std::chrono::high_resolution_clock::now();
        
        // === 8. Write completion metadata ===
        int active_count = 0;
        for (const auto& ion : result.ions) {
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
                std::cout << "\n[WARNING] Simulation incomplete: " << active_count 
                          << " ions still active at t=" << t_end << " s\n";
                std::cout << "Consider using continue mode:\n";
                std::cout << "  \"continue_from\": \"" << gParams.output_file << ".h5\",\n";
                std::cout << "  \"continue_time_s\": " << gParams.continue_time_s << "\n\n";
            }
        } catch (const std::exception& e) {
            std::cerr << "Warning: Failed to write metadata: " << e.what() << "\n";
        }
        
        // === 10. Optional result printout ===
        if (gParams.print_results) {
            ICARION::utils::print_results(result, 100);
        }
        logger.log("HDF5 file written successfully.");
        double elapsed_s = std::chrono::duration<double>(end - start).count();
        logger.log("Simulation completed in " + std::to_string(elapsed_s) + " s CPU time.");

        logger.finalize(result, gParams.output_file);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}

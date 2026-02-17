// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include <chrono>
#include <exception>
#include <iostream>
#include <algorithm>
#include <filesystem>
#include <fstream>
#include <cctype>

#ifdef _OPENMP
#include <omp.h>
#endif

// Core systems (Phase 5: SSOT-compliant)
#include "core/config/loader/ConfigLoader.h"
#include "core/config/utils/ConfigOverride.h"
#include "core/integrator/SimulationEngine.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/log/Logger.h"
#include "core/utils/startupBanner.h"
#include "core/utils/Profiler.h"
#include "utils/cli_parser.h"
#include "main/setup/PhysicsSetup.h"
#include <nlohmann/json.hpp>

/**
 * @file main.cpp
 * @brief Entry point for the ICARION simulation.
 *
 * Initializes and executes a complete ICARION run using SimulationEngine.
 * Also supports info/validation modes that exit before running a simulation.
 *
 * @param[in] argc Number of command-line arguments.
 * @param[in] argv Command-line arguments (expects JSON configuration path).
 * @return 0 on success (including info/validation-only modes), non-zero on
 *         validation failures or runtime errors.
 *
 * @note The JSON file must define a complete FullConfig structure; CLI overrides
 *       are applied on top of it before validation.
 *
 * @see config::ConfigLoader
 * @see integrator::SimulationEngine
 */

// ------------------------------------------------------------------------
// Helpers
// ------------------------------------------------------------------------
static bool parse_bool_cli(const std::string& value) {
    std::string lower = value;
    std::transform(lower.begin(), lower.end(), lower.begin(), ::tolower);
    return (lower == "true" || lower == "1" || lower == "yes" || lower == "on");
}

static void apply_override_to_json(nlohmann::json& j, const std::string& key, const std::string& value) {
    try {
        if (key == "simulation.dt_s" || key == "simulation.timestep") {
            j["simulation"]["dt_s"] = std::stod(value);
        } else if (key == "simulation.total_time_s" || key == "simulation.total_time") {
            j["simulation"]["total_time_s"] = std::stod(value);
        } else if (key == "simulation.write_interval") {
            j["simulation"]["write_interval"] = std::stoi(value);
        } else if (key == "simulation.rng_seed" || key == "simulation.seed") {
            j["simulation"]["rng_seed"] = std::stoul(value);
        } else if (key == "simulation.integrator") {
            j["simulation"]["integrator"] = value;
        } else if (key == "simulation.enable_gpu") {
            j["simulation"]["enable_gpu"] = parse_bool_cli(value);
        } else if (key == "simulation.enable_openmp") {
            j["simulation"]["enable_openmp"] = parse_bool_cli(value);
        }
        else if (key == "physics.collision_model") {
            j["physics"]["collision_model"] = value;
        } else if (key == "physics.enable_reactions") {
            j["physics"]["enable_reactions"] = parse_bool_cli(value);
        } else if (key == "physics.enable_space_charge") {
            j["physics"]["enable_space_charge"] = parse_bool_cli(value);
        } else if (key == "physics.enable_space_charge_gpu") {
            j["physics"]["enable_space_charge_gpu"] = parse_bool_cli(value);
        } else if (key == "physics.enable_ou_thermalization") {
            j["physics"]["enable_ou_thermalization"] = parse_bool_cli(value);
        }
        else if (key == "output.folder") {
            j["output"]["folder"] = value;
        } else if (key == "output.trajectory_file" || key == "output.file") {
            j["output"]["trajectory_file"] = value;
        } else if (key == "output.print_progress") {
            j["output"]["print_progress"] = parse_bool_cli(value);
        } else if (key == "output.buffer_byte_cap") {
            j["output"]["buffer_byte_cap"] = std::stoull(value);
        }
        else if (key == "species_database" || key == "database.species") {
            j["species_database"] = value;
        } else if (key == "reaction_database" || key == "database.reactions") {
            j["reaction_database"] = value;
        } else {
            // Unknown override key -> ignore for snapshot
        }
    } catch (const std::exception& e) {
        ICARION::log::Logger::main()->warn("Config snapshot override skipped for {}: {}", key, e.what());
    }
}

static std::string write_config_snapshot(
    const std::string& config_path,
    const ICARION::config::FullConfig& config,
    const ICARION::cli::CLIOptions& opts
) {
    // Snapshot is based on the original JSON plus CLI overrides (no derived/finalized fields).
    // Read original JSON
    nlohmann::json j;
    try {
        std::ifstream in(config_path);
        if (!in) {
            ICARION::log::Logger::main()->warn("Config snapshot skipped: cannot read {}", config_path);
            return {};
        }
        in >> j;
    } catch (const std::exception& e) {
        ICARION::log::Logger::main()->warn("Config snapshot skipped: failed to parse {} ({})", config_path, e.what());
        return {};
    }
    
    // Apply CLI overrides (same keys as ConfigOverride)
    for (const auto& [key, value] : opts.overrides) {
        apply_override_to_json(j, key, value);
    }
    
    // Apply direct CLI flags (not in overrides map)
    if (opts.output_file.has_value()) {
        j["output"]["trajectory_file"] = opts.output_file.value();
    }
    if (opts.output_dir.has_value()) {
        j["output"]["folder"] = opts.output_dir.value();
    }
    if (opts.buffer_byte_cap.has_value()) {
        j["output"]["buffer_byte_cap"] = opts.buffer_byte_cap.value();
    }
    if (opts.seed.has_value()) {
        j["simulation"]["rng_seed"] = opts.seed.value();
    }
    if (opts.no_reactions) {
        j["physics"]["enable_reactions"] = false;
    }
    
    // Snapshot path: same folder as output, base of trajectory file + ".config.json"
    std::filesystem::path out_dir = config.output.folder;
    std::filesystem::path traj_file = config.output.trajectory_file;
    std::string base = traj_file.stem().string();
    if (base.empty()) {
        base = "config_snapshot";
    }
    std::filesystem::create_directories(out_dir);
    std::filesystem::path snapshot_path = out_dir / (base + ".config.json");
    
    try {
        std::string content = j.dump(2);
        std::ofstream out(snapshot_path);
        out << content;
        out.close();
        ICARION::log::Logger::main()->info("Wrote config snapshot: {}", snapshot_path.string());
        return content;
    } catch (const std::exception& e) {
        ICARION::log::Logger::main()->warn("Failed to write config snapshot to {} ({})", snapshot_path.string(), e.what());
        throw;
    }
}

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
        
        // Set OpenMP thread count if specified (before logger init)
        if (opts.num_threads.has_value()) {
#ifdef _OPENMP
            omp_set_num_threads(opts.num_threads.value());
#endif
        }
        
        // === 2. Initialize logging ===
        {
            PROFILE_SCOPE("Logging Initialization");
            log::Logger::init(
                opts.log_level,
                opts.log_file.value_or(""),
                opts.log_format);
        }
        
        // Log thread count after logger is initialized
        if (opts.num_threads.has_value()) {
#ifdef _OPENMP
            log::Logger::main()->info("OpenMP threads set to: {}", opts.num_threads.value());
#else
            log::Logger::main()->warn("--threads specified but OpenMP is not enabled");
#endif
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
        
        // === Persist effective config snapshot (after overrides) ===
        config.resolved_config_json = write_config_snapshot(opts.config_file, config, opts);
        
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
        core::IonEnsemble ensemble = core::IonEnsemble::from_legacy(ions);

        // Prepopulate species pool/index to make reaction processing thread-safe.
        // Reactions can create product ions of species that are not present initially.
        {
            std::vector<std::string> all_species_ids;
            all_species_ids.reserve(config.species_db.species.size());
            for (const auto& kv : config.species_db.species) {
                all_species_ids.push_back(kv.first);
            }
            ensemble.prepopulate_species_pool(all_species_ids);
        }
        
        // === 5. Create physics dependencies ===
        setup::PhysicsModules physics;
        {
            PROFILE_SCOPE("Physics Module Setup");
            physics = setup::PhysicsSetup::initialize(config, ensemble);
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
        
        core::IonEnsemble final_ensemble;
        {
            PROFILE_SCOPE("Simulation Run");
            final_ensemble = engine.run(ensemble);
        }
        
        auto end = std::chrono::high_resolution_clock::now();
        double elapsed_s = std::chrono::duration<double>(end - start).count();
        
        // === 8. Report results ===
        size_t active_count = 0;
        const auto* active_ptr = final_ensemble.active_data();
        const auto* born_ptr = final_ensemble.born_data();
        for (size_t i = 0; i < final_ensemble.size(); ++i) {
            if (active_ptr[i] && born_ptr[i]) {
                ++active_count;
            }
        }
        
        log::Logger::main()->info("");
        log::Logger::main()->info("=== Simulation Complete ===");
        log::Logger::main()->info("CPU time:     {:.3f} s", elapsed_s);
        log::Logger::main()->info("Active ions:  {}/{}", active_count, final_ensemble.size());
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

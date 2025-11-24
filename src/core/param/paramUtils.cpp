/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        paramUtils.cpp
 *   @brief       Implementation of parameter loading utilities for ICARION.
 *
 * @details
 * Provides routines to read, validate, and initialize simulation parameters
 * from a JSON configuration file. Supports both flat and nested JSON layouts.
 * Constructs derived quantities such as the simulation time array (`t_eval`)
 * and maps string identifiers (collision models) to enumerated types.
 *
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include "core/param/paramUtils.h"

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <algorithm>
#include <filesystem>

namespace ICARION {
namespace core {

/**
 * @brief Load global simulation parameters from a JSON file.
 *
 * This routine handles both common layouts:
 * - Nested: {"globalParams": { ... }} or {"global_params": { ... }}
 * - Flat: top-level keys (dt_s, write_interval, output_file, etc.)
 *
 * Missing keys are filled with reasonable defaults. Mandatory numeric
 * values (e.g., dt_s) are sanity-checked. Throws a runtime_error on
 * file access failure or invalid parameter values.
 *
 * @param[in] filename Path to the JSON configuration file.
 * @return GlobalParams Fully populated global parameters structure.
 * @throws std::runtime_error if the file cannot be opened or parameters are invalid.
 */

GlobalParams load_global_params(const std::string& filename) {
     // --- Open JSON file ---
    std::ifstream file(filename);
    if (!file) throw std::runtime_error("Cannot open JSON: " + filename);

    // --- Parse JSON ---
    Json::Value root;
    try {
        file >> root;
    } catch (const Json::Exception& e) {
        throw std::runtime_error("JSON parse error in " + filename + ": " + e.what());
    }

    // Validate top-level structure
    if (!root.isMember("simulation")) {
        throw std::runtime_error("Missing required 'simulation' section in " + filename);
    }
    if (!root.isMember("physics")) {
        throw std::runtime_error("Missing required 'physics' section in " + filename);
    }
    if (!root.isMember("output")) {
        throw std::runtime_error("Missing required 'output' section in " + filename);
    }

    const Json::Value& sim = root["simulation"];
    const Json::Value& phys = root["physics"];
    const Json::Value& out = root["output"];

    GlobalParams g{}; ///< output structure 

    // --- Lambda helpers for value retrieval with validation ---
    auto getBool = [](const Json::Value& section, const char* key, bool def, const std::string& section_name) {
        if (!section.isMember(key)) return def;
        if (!section[key].isBool()) {
            throw std::runtime_error("Parameter '" + std::string(key) + "' in section '" + section_name + "' must be a boolean");
        }
        return section[key].asBool();
    };
    
    auto getInt = [](const Json::Value& section, const char* key, int def, const std::string& section_name) {
        if (!section.isMember(key)) return def;
        if (!section[key].isInt()) {
            throw std::runtime_error("Parameter '" + std::string(key) + "' in section '" + section_name + "' must be an integer");
        }
        return section[key].asInt();
    };
    
    auto getDouble = [](const Json::Value& section, const char* key, double def, const std::string& section_name) {
        if (!section.isMember(key)) return def;
        if (!section[key].isNumeric()) {
            throw std::runtime_error("Parameter '" + std::string(key) + "' in section '" + section_name + "' must be a number");
        }
        return section[key].asDouble();
    };
    
    auto getString = [](const Json::Value& section, const char* key, const char* def, const std::string& section_name) {
        if (!section.isMember(key)) return std::string(def);
        if (!section[key].isString()) {
            throw std::runtime_error("Parameter '" + std::string(key) + "' in section '" + section_name + "' must be a string");
        }
        return section[key].asString();
    };

    // --- Load simulation parameters ---
    // Time parameters
    if (!sim.isMember("total_time_s") && !sim.isMember("dt_s")) {
        throw std::runtime_error("Missing required time parameters: 'total_time_s' and 'dt_s' in 'simulation' section");
    }
    
    double total_time_s = getDouble(sim, "total_time_s", 1e-3, "simulation");
    g.dt_s = getDouble(sim, "dt_s", 1e-9, "simulation");
    
    if (g.dt_s <= 0.0) {
        throw std::runtime_error("dt_s must be positive, got: " + std::to_string(g.dt_s));
    }
    if (total_time_s <= 0.0) {
        throw std::runtime_error("total_time_s must be positive, got: " + std::to_string(total_time_s));
    }
    
    g.sim_time_steps = static_cast<int>(std::ceil(total_time_s / g.dt_s));
    g.write_interval = getInt(sim, "write_interval", 100, "simulation");
    
    // Number of ions
    if (root.isMember("ions") && root["ions"].isObject()) {
        const auto& ions = root["ions"];
        g.num_ions = getInt(ions, "count", 1, "ions");
        g.ion_cloud_file = getString(ions, "initial_distribution", "ioncloud.json", "ions");
    } else {
        g.num_ions = 1;
        g.ion_cloud_file = "ioncloud.json";
    }
    
    // Integrator selection
    std::string integrator = getString(sim, "integrator", "RK4", "simulation");
    std::transform(integrator.begin(), integrator.end(), integrator.begin(), 
                   [](unsigned char c){ return std::tolower(c); });
    // Store for later use with domains (will be applied in load_single_domain)
    
    // GPU and parallelization
    g.enable_gpu = getBool(sim, "enable_gpu", false, "simulation");
    g.parallelization = getBool(sim, "enable_openmp", false, "simulation");
    
    // RNG seed
    g.rng_seed = getInt(sim, "rng_seed", 42, "simulation");
    
    // Reactions and space charge
    g.enable_reactions = getBool(sim, "enable_reactions", false, "simulation");
    g.enable_space_charge = getBool(sim, "enable_space_charge", false, "simulation");
    
    // Continue mode
    g.continue_from = getString(sim, "continue_from", "", "simulation");
    g.continue_time_s = getDouble(sim, "continue_time_s", 0.0, "simulation");
    g.auto_continue_if_active = getBool(sim, "auto_continue_if_active", false, "simulation");
    
    // RK45 adaptive stepping parameters (for numerical safety)
    if (sim.isMember("rk45_settings") && sim["rk45_settings"].isObject()) {
        const Json::Value& rk45 = sim["rk45_settings"];
        g.rk45_settings.absTol = getDouble(rk45, "abs_tol", 1e-14, "rk45_settings");
        g.rk45_settings.relTol = getDouble(rk45, "rel_tol", 1e-12, "rk45_settings");
        g.rk45_settings.dt_min = getDouble(rk45, "dt_min", 1e-12, "rk45_settings");
        g.rk45_settings.safety = getDouble(rk45, "safety", 0.84, "rk45_settings");
        g.rk45_settings.min_factor = getDouble(rk45, "min_factor", 0.2, "rk45_settings");
        g.rk45_settings.max_factor = getDouble(rk45, "max_factor", 2.0, "rk45_settings");
        g.rk45_settings.max_rejects = getInt(rk45, "max_rejects", 1000, "rk45_settings");
        
        // Validation
        if (g.rk45_settings.absTol <= 0.0) {
            throw std::runtime_error("rk45_settings.abs_tol must be positive, got: " + std::to_string(g.rk45_settings.absTol));
        }
        if (g.rk45_settings.relTol <= 0.0) {
            throw std::runtime_error("rk45_settings.rel_tol must be positive, got: " + std::to_string(g.rk45_settings.relTol));
        }
        if (g.rk45_settings.dt_min <= 0.0) {
            throw std::runtime_error("rk45_settings.dt_min must be positive, got: " + std::to_string(g.rk45_settings.dt_min));
        }
        if (g.rk45_settings.max_rejects < 1) {
            throw std::runtime_error("rk45_settings.max_rejects must be at least 1, got: " + std::to_string(g.rk45_settings.max_rejects));
        }
    }
    
    // Numerical safety settings (for NaN/Inf detection and stability)
    if (sim.isMember("numerical_safety") && sim["numerical_safety"].isObject()) {
        const Json::Value& safety = sim["numerical_safety"];
        g.numerical_safety.enable_nan_checks = getBool(safety, "enable_nan_checks", true, "numerical_safety");
        g.numerical_safety.enable_bounds_checks = getBool(safety, "enable_bounds_checks", false, "numerical_safety");
        g.numerical_safety.enable_logging = getBool(safety, "enable_logging", false, "numerical_safety");
        g.numerical_safety.throw_on_violation = getBool(safety, "throw_on_violation", true, "numerical_safety");
        g.numerical_safety.attempt_recovery = getBool(safety, "attempt_recovery", false, "numerical_safety");
        
        // Bounds settings
        g.numerical_safety.max_position_m = getDouble(safety, "max_position_m", 1.0, "numerical_safety");
        g.numerical_safety.max_velocity_ms = getDouble(safety, "max_velocity_ms", 1e6, "numerical_safety");
        g.numerical_safety.max_acceleration_ms2 = getDouble(safety, "max_acceleration_ms2", 1e12, "numerical_safety");
        
        // Validation
        if (g.numerical_safety.max_position_m <= 0.0) {
            throw std::runtime_error("numerical_safety.max_position_m must be positive, got: " + std::to_string(g.numerical_safety.max_position_m));
        }
        if (g.numerical_safety.max_velocity_ms <= 0.0) {
            throw std::runtime_error("numerical_safety.max_velocity_ms must be positive, got: " + std::to_string(g.numerical_safety.max_velocity_ms));
        }
        if (g.numerical_safety.max_acceleration_ms2 <= 0.0) {
            throw std::runtime_error("numerical_safety.max_acceleration_ms2 must be positive, got: " + std::to_string(g.numerical_safety.max_acceleration_ms2));
        }
    }
    
    // --- Load physics parameters ---
    // Collision model
    std::string collision_str = getString(phys, "collision_model", "EHSS", "physics");
    std::transform(collision_str.begin(), collision_str.end(), collision_str.begin(), 
                   [](unsigned char c){ return std::tolower(c); });
    
    if (collision_str == "ehss" || collision_str == "EHSS") g.collisionModel = CollisionModel::EHSS;
    else if (collision_str == "HSS" || collision_str == "HSS" || collision_str == "hss" || collision_str == "HSS") g.collisionModel = CollisionModel::HSS;
    else if (collision_str == "friction" || collision_str == "Friction") g.collisionModel = CollisionModel::Friction;
    else if (collision_str == "langevin" || collision_str == "Langevin") g.collisionModel = CollisionModel::Langevin;
    else if (collision_str == "hardsphere" || collision_str == "hard_sphere" || collision_str == "HardSphere" || collision_str == "Hard_Sphere") 
        g.collisionModel = CollisionModel::HSD;
    else if (collision_str == "nocollisions" || collision_str == "none" || collision_str == "NoCollisions" || collision_str == "No_Collisions") 
        g.collisionModel = CollisionModel::NoCollisions;
    else {
        throw std::runtime_error("Unknown collision model: '" + collision_str + 
                                 "'. Valid options: EHSS, HSS, Friction, Langevin, HardSphere, NoCollisions");
    }
    
    // Thermalization
    g.enable_ou_thermalization = getBool(phys, "enable_ou_thermalization", false, "physics");
    g.force_ou_for_stochastic_models = getBool(phys, "force_ou_for_stochastic", false, "physics");
    
    // --- Load output parameters ---
    std::string output_folder = getString(out, "folder", "results", "output");
    std::string trajectory_file = getString(out, "trajectory_file", "trajectories.h5", "output");
    
    // Combine folder and filename
    if (!output_folder.empty() && output_folder.back() != '/') {
        output_folder += "/";
    }
    g.output_file = output_folder + trajectory_file;
    
    g.print_results = getBool(out, "print_progress", false, "output");
    
    // --- Load database files ---
    // species input system: species_database 
    g.species_database_file = getString(root, "species_database", "", "root");
    
    // reaction input system: reaction_database 
    g.reaction_file = getString(root, "reaction_database", "reactions.json", "root");
    
    if (root.isMember("species_database")) {
        // Species database path can be used to infer reaction file location
        std::string species_db = getString(root, "species_database", "", "root");
        if (!species_db.empty() && g.reaction_file == "reactions.json") {
            // Try to find reaction file next to species file
            std::filesystem::path species_path(species_db);
            std::filesystem::path reaction_path = species_path.parent_path() / "reactions.json";
            if (std::filesystem::exists(reaction_path)) {
                g.reaction_file = reaction_path.string();
            }
        }
    }
    
    //LEGACY: should be defined in the species database now! -> REMOVE LATER
    g.geometry_file = getString(root, "geometry_file", "", "root");

    // If geometry_file not specified, try to find a standard 'geometry.json'
    if (g.geometry_file.empty()) {
        try {
            std::filesystem::path cfgPath(filename);
            std::filesystem::path candidate = cfgPath.parent_path() / "geometry.json";
            if (!candidate.empty() && std::filesystem::exists(candidate)) {
                g.geometry_file = candidate.string();
            } else {
                candidate = std::filesystem::current_path() / "geometry.json";
                if (std::filesystem::exists(candidate)) {
                    g.geometry_file = candidate.string();
                }
            }
        } catch (const std::exception& e) {
            // If filesystem isn't available or an error occurs, leave geometry_file empty.
            (void)e;
        }
    }

    // Tests fallback: if still empty and TEST_DATA_DIR is available, use repo test geometry
#ifdef TEST_DATA_DIR
    if (g.geometry_file.empty()) {
        try {
            std::filesystem::path testGeom = std::filesystem::path(TEST_DATA_DIR) / "io_tests" / "repo_geometry.json";
            if (std::filesystem::exists(testGeom)) {
                g.geometry_file = testGeom.string();
            }
        } catch (...) {
            // ignore
        }
    }
#endif
    
    g.input_file = filename;

    // --- Build simulation time array ---
    g.t_eval.resize(g.sim_time_steps);
    for (size_t i = 0; i < g.sim_time_steps; ++i) {
        g.t_eval[i] = i * g.dt_s;
    }

    // --- Sanity checks ---
    run_guard_check_global(g);

    return g;
}

/**
 * @brief Lists the required parameters for each instrument type.
 *
 * This map defines which configuration parameters must be provided
 * for a given instrument in order to initialize it properly.
 * 
 */
std::unordered_map<Instrument, std::vector<std::string>> required_params{
    {Instrument::TOF,        {"DC.axial_V", "geom.acc_length_m", "geom.length_m", "geom.radius_m"}},
    {Instrument::LQIT,       {"RF.voltage_V", "RF.frequency_Hz", "DC.axial_V", "geom.radius_m", "geom.length_m"}},
    {Instrument::Orbitrap,   {"geom.radius_in_m", "geom.radius_out_m", "geom.radius_char_m", "DC.radial_V"}},
    {Instrument::IMS,        {"DC.EN_Td", "geom.length_m", "geom.radius_m"}},
    {Instrument::QuadrupoleRF, {"RF.voltage_V", "RF.frequency_Hz", "geom.radius_m", "geom.length_m"}}
};

/**
 * @brief Checks that all required parameters for the given instrument domain are set.
 *
 * Ensures that all required parameters are nonzero for the given instrument type.
 * Also validates common environmental parameters for all instruments.
 */
void sanity_check_domain(const InstrumentDomain& dom) {
    // Build lookup map of all numeric domain parameters
    std::unordered_map<std::string, double> dom_params = {
        {"geom.acc_length_m", dom.geom.acc_length_m},
        {"geom.length_m", dom.geom.length_m},
        {"geom.radius_m", dom.geom.radius_m},
        {"geom.radius_in_m", dom.geom.radius_in_m},
        {"geom.radius_out_m", dom.geom.radius_out_m},
        {"geom.radius_char_m", dom.geom.radius_char_m},
        {"RF.voltage_V", dom.RF.voltage_V},
        {"RF.frequency_Hz", dom.RF.frequency_Hz},
        {"DC.axial_V", dom.DC.axial_V},
        {"DC.radial_V", dom.DC.radial_V},
        {"DC.EN_Td", dom.DC.EN_Td},
        {"env.pressure_Pa", dom.env.pressure_Pa},
        {"env.temperature_K", dom.env.temperature_K}
    };

    // Find required keys for this instrument
    auto it = required_params.find(dom.instrument);
    if (dom.FA_file.empty() && dom.FA_terms.empty() && it != required_params.end()) {
        for (const auto& key : it->second) {
            auto found = dom_params.find(key);
            if (found == dom_params.end() || found->second <= 0.0) {
                std::cerr << "Warning: Required parameter '" << key 
                          << "' not set for instrument:"
                          << static_cast<int>(dom.instrument) << ". Check if this is correct. \n";
            }
        }
    } else if (!dom.FA_file.empty() || !dom.FA_terms.empty()) {
        // Optional: Minimal geometry/environment check for PA arrays
        if (dom.geom.length_m <= 0.0 || dom.geom.radius_m <= 0.0)
            std::cerr << "⚠️  Warning: FA_file set, but geometry values are zero — "
                        "check that grid extent matches physical dimensions.\n";
    }

    // Global sanity checks (apply to all instruments)
    if (dom.env.pressure_Pa <= 0.0 || dom.env.temperature_K <= 0.0) {
        throw std::runtime_error("Invalid environmental parameters: pressure_Pa or temperature_K not set.");
    }

    // Optional: Check for gas velocity definition
    if (std::isnan(dom.env.gas_velocity_m_s.x) ||
        std::isnan(dom.env.gas_velocity_m_s.y) ||
        std::isnan(dom.env.gas_velocity_m_s.z)) {
        throw std::runtime_error("Gas velocity vector not properly defined.");
    }

    // Optional verbose output (for debugging)
    // std::cout << "Sanity check passed for instrument " << static_cast<int>(dom.instrument) << std::endl;
}

/**
 * @brief Load parameters specific to a single instrument domain from a JSON object.
 *
 * Parses geometrical, electrical, and environmental parameters for a single
 * instrument domain. Initializes derived quantities like particle density,
 * mean thermal velocity, and angular frequencies. Supports DC, RF, and AC
 * voltage profiles, as well as optional axis transformations.
 *
 * @param[in] j JSON object representing a single instrument domain.
 * @return InstrumentDomain Fully populated domain parameters.
 * @throws std::runtime_error if required parameters are missing or invalid.
 *
 * @note After loading, the domain is validated with `sanity_check_domain`.
 *       Multiple domains can be loaded independently while sharing global parameters.
 */
InstrumentDomain load_single_domain(const Json::Value& j) {
    InstrumentDomain dom;

    // Validate required fields
    if (!j.isMember("instrument")) {
        throw std::runtime_error("Domain missing required field: 'instrument'");
    }
    if (!j.isMember("geometry")) {
        throw std::runtime_error("Domain missing required field: 'geometry'");
    }
    if (!j.isMember("environment")) {
        throw std::runtime_error("Domain missing required field: 'environment'");
    }

    // Lambda helpers for validated retrieval
    auto getBool = [](const Json::Value& section, const char* key, bool def, const std::string& context) {
        if (!section.isMember(key)) return def;
        if (!section[key].isBool()) {
            throw std::runtime_error("Parameter '" + std::string(key) + "' in " + context + " must be a boolean");
        }
        return section[key].asBool();
    };
    
    auto getDouble = [](const Json::Value& section, const char* key, double def, const std::string& context) {
        if (!section.isMember(key)) return def;
        if (!section[key].isNumeric()) {
            throw std::runtime_error("Parameter '" + std::string(key) + "' in " + context + " must be a number");
        }
        return section[key].asDouble();
    };
    
    auto getString = [](const Json::Value& section, const char* key, const char* def, const std::string& context) {
        if (!section.isMember(key)) return std::string(def);
        if (!section[key].isString()) {
            throw std::runtime_error("Parameter '" + std::string(key) + "' in " + context + " must be a string");
        }
        return section[key].asString();
    };

    // --- Identify instrument type ---
    std::string instr_str = j["instrument"].asString();
    if (instr_str == "LQIT" || instr_str == "lqit") dom.instrument = Instrument::LQIT;
    else if (instr_str == "IMS" || instr_str == "ims" || instr_str == "SIFDT-MS" || instr_str == "SIFDT_MS" || instr_str == "sifdt-ms" || instr_str == "sifdt_ms") dom.instrument = Instrument::IMS;
    else if (instr_str == "Quadrupole" || instr_str == "QuadrupoleRF" || instr_str == "quadrupole" || instr_str == "quadrupoleRF") dom.instrument = Instrument::QuadrupoleRF;
    else if (instr_str == "TOF" || instr_str == "tof") dom.instrument = Instrument::TOF;
    else if (instr_str == "Orbitrap" || instr_str == "orbitrap") dom.instrument = Instrument::Orbitrap;
    else if (instr_str == "FT-ICR" || instr_str == "FT_ICR" || instr_str == "ft-icr" || instr_str == "ft_icr") dom.instrument = Instrument::FTICR;
    else if (instr_str == "NoFixedInstrument" || instr_str == "nofixedinstrument" || instr_str == "No_Fixed_Instrument") dom.instrument = Instrument::NoFixedInstrument;
    else {
        throw std::runtime_error("Unknown instrument type: '" + instr_str + "'. Valid types: LQIT, SIFDT-MS, IMS, Quadrupole, TOF, Orbitrap, FT-ICR");
    }

    // --- Geometry ---
    const auto& jGeom = j["geometry"];
    dom.geom.length_m = getDouble(jGeom, "length_m", 0.0, "geometry");
    dom.geom.radius_m = getDouble(jGeom, "radius_m", 0.0, "geometry");
    dom.geom.radius_in_m = getDouble(jGeom, "radius_in_m", 0.0, "geometry");
    dom.geom.radius_out_m = getDouble(jGeom, "radius_out_m", 0.0, "geometry");
    dom.geom.radius_char_m = getDouble(jGeom, "radius_char_m", 0.0, "geometry");
    dom.geom.acc_length_m = getDouble(jGeom, "acc_length_m", 0.0, "geometry");
    dom.geom.end_aperture_m = getDouble(jGeom, "end_aperture_m", 1.0, "geometry");

    // NOTE: Orbitrap hyperbolic boundary constants moved to OrbitrapInstrument class
    // (paramUtils.cpp is legacy code, kept for backward compatibility only)
    // dom.geom.orbitrap_C_in  = -0.5 * dom.geom.radius_in_m  * dom.geom.radius_in_m;
    // dom.geom.orbitrap_C_out = -0.5 * dom.geom.radius_out_m * dom.geom.radius_out_m;

    // Optional origin transform
    if (jGeom.isMember("origin_m") && jGeom["origin_m"].isArray() && jGeom["origin_m"].size() == 3) {
        dom.geom.origin_m = Vec3{ 
            jGeom["origin_m"][0].asDouble(), 
            jGeom["origin_m"][1].asDouble(), 
            jGeom["origin_m"][2].asDouble() 
        };
    } else {
        dom.geom.origin_m = Vec3{0.0, 0.0, 0.0};
    }

    dom.geom.min_bound = dom.geom.origin_m - Vec3{dom.geom.radius_m, dom.geom.radius_m, dom.geom.length_m/2};
    dom.geom.max_bound = dom.geom.origin_m + Vec3{dom.geom.radius_m, dom.geom.radius_m, dom.geom.length_m/2};

    // --- Environment ---
    const auto& jEnv = j["environment"];
    dom.env.pressure_Pa = getDouble(jEnv, "pressure_Pa", 101325.0, "environment");
    dom.env.temperature_K = getDouble(jEnv, "temperature_K", 300.0, "environment");
    
    if (dom.env.pressure_Pa <= 0.0) {
        throw std::runtime_error("Environment pressure_Pa must be positive, got: " + std::to_string(dom.env.pressure_Pa));
    }
    if (dom.env.temperature_K <= 0.0) {
        throw std::runtime_error("Environment temperature_K must be positive, got: " + std::to_string(dom.env.temperature_K));
    }
    
    dom.env.particle_density_m_3 = dom.env.pressure_Pa / (BOLTZMANN_CONSTANT * dom.env.temperature_K);
    
    // Gas velocity
    if (jEnv.isMember("gas_velocity_m_s") && jEnv["gas_velocity_m_s"].isArray() && jEnv["gas_velocity_m_s"].size() == 3) {
        dom.env.gas_velocity_m_s = Vec3{
            jEnv["gas_velocity_m_s"][0].asDouble(),
            jEnv["gas_velocity_m_s"][1].asDouble(),
            jEnv["gas_velocity_m_s"][2].asDouble()
        };
    } else {
        dom.env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    }

    dom.env.mean_thermal_velocity_m_s = std::sqrt(8*BOLTZMANN_CONSTANT*dom.env.temperature_K/M_PI/MOLAR_MASS_HE_KG);

    // Neutral parameters
    dom.env.neutral_species_id = getString(jEnv, "gas_species", "He", "environment");
    std::string gas = dom.env.neutral_species_id;
    
    if (gas == "He" || gas == "helium" || gas == "Helium") {
        dom.env.neutral_mass_kg = MOLAR_MASS_HE_KG;
        dom.env.neutral_polarizability_m3 = POLARIZABILITY_HE_SI;
    } 
    else if (gas == "N2" || gas == "nitrogen" || gas == "Nitrogen") {
        dom.env.neutral_mass_kg = MOLAR_MASS_N2_KG;
        dom.env.neutral_polarizability_m3 = POLARIZABILITY_N2_SI;
    } 
    else if (gas == "Ar" || gas == "argon" || gas == "Argon") {
        dom.env.neutral_mass_kg = MOLAR_MASS_AR_KG;
        dom.env.neutral_polarizability_m3 = POLARIZABILITY_AR_SI;
    } 
    else if (gas == "CO2" || gas == "carbon_dioxide" || gas == "Carbon_Dioxide") {
        dom.env.neutral_mass_kg = MOLAR_MASS_CO2_KG;
        dom.env.neutral_polarizability_m3 = POLARIZABILITY_CO2_SI;
    } 
    else if (gas == "Ne" || gas == "neon" || gas == "Neon") {
        dom.env.neutral_mass_kg = MOLAR_MASS_NE_KG;
        dom.env.neutral_polarizability_m3 = POLARIZABILITY_NE_SI;
    }
    else if (gas == "Air" || gas == "air") {
        // Effective properties for air (approx. 78% N2, 21% O2, 1% Ar)
        dom.env.neutral_mass_kg = 0.78 * MOLAR_MASS_N2_KG + 0.21 * MOLAR_MASS_O2_KG + 0.01 * MOLAR_MASS_AR_KG;
        dom.env.neutral_polarizability_m3 = 0.78 * POLARIZABILITY_N2_SI + 0.21 * POLARIZABILITY_O2_SI + 0.01 * POLARIZABILITY_AR_SI;
    }
    else {
        throw std::runtime_error("Unknown gas species: '" + gas + "'. Valid options: He, N2, Air");
    }
    // --- DC, RF, AC, B field voltages ---
    // Check if 'fields' section exists (new format) or use legacy top-level fields
    const Json::Value* fieldsSection = nullptr;
    if (j.isMember("fields") && j["fields"].isObject()) {
        fieldsSection = &j["fields"];
    }
    
    // DC voltages
    if (fieldsSection && fieldsSection->isMember("DC")) {
        const auto& jDC = (*fieldsSection)["DC"];
        dom.DC.EN_Td = getDouble(jDC, "EN_Td", 0.0, "fields.DC");
        dom.DC.EN_Vm2 = dom.DC.EN_Td * 1e-21;
        
        if (dom.DC.EN_Vm2 != 0.0) {
            dom.DC.axial_V = dom.DC.EN_Vm2 * dom.env.particle_density_m_3 * dom.geom.length_m;
        } else {
            dom.DC.axial_V = getDouble(jDC, "axial_V", 0.0, "fields.DC");
        }
        
        dom.DC.quad_V = getDouble(jDC, "quad_V", 0.0, "fields.DC");
        dom.DC.radial_V = getDouble(jDC, "radial_V", 0.0, "fields.DC");
        dom.DC.enable_radial_voltage_sweep = getBool(jDC, "enable_radial_sweep", false, "fields.DC");
        dom.DC.radial_slope_V_s = getDouble(jDC, "radial_slope_V_s", 0.0, "fields.DC");
        dom.DC.radial_start_time_s = getDouble(jDC, "radial_start_time_s", 0.0, "fields.DC");
        dom.DC.radial_rise_time_s = getDouble(jDC, "radial_rise_time_s", 0.0, "fields.DC");
    }

    // RF voltages
    if (fieldsSection && fieldsSection->isMember("RF")) {
        const auto& jRF = (*fieldsSection)["RF"];
        dom.RF.voltage_V = getDouble(jRF, "voltage_V", 0.0, "fields.RF");
        dom.RF.frequency_Hz = getDouble(jRF, "frequency_Hz", 0.0, "fields.RF");
        dom.RF.angular_frequency_rad_s = dom.RF.frequency_Hz * 2.0 * M_PI;
        dom.RF.phase_rad = getDouble(jRF, "phase_rad", 0.0, "fields.RF");
    }

    // AC voltages
    if (fieldsSection && fieldsSection->isMember("AC")) {
        const auto& jAC = (*fieldsSection)["AC"];
        dom.AC.voltage_V = getDouble(jAC, "voltage_V", 0.0, "fields.AC");
        dom.AC.frequency_Hz = getDouble(jAC, "frequency_Hz", 0.0, "fields.AC");
        
        if (dom.AC.voltage_V != 0.0 && dom.AC.frequency_Hz == 0.0) {
            std::cerr << "Warning: AC voltage set but frequency is zero. Using pseudopotential approximation (m/z = 100).\n";
            double q_param = 4.0 * ELEM_CHARGE_C * dom.AC.voltage_V / 
                            (100.0 * AMU_TO_KG * dom.geom.radius_m * dom.geom.radius_m * 
                             dom.RF.angular_frequency_rad_s * dom.RF.angular_frequency_rad_s);
            double beta = std::sqrt(0.5 * q_param * q_param);
            dom.AC.frequency_Hz = (1.0 / (2.0 * M_PI)) * beta * dom.RF.angular_frequency_rad_s;
            std::cerr << " -> Estimated AC frequency: " << dom.AC.frequency_Hz << " Hz\n";
        }
        
        dom.AC.angular_frequency_rad_s = dom.AC.frequency_Hz * 2.0 * M_PI;
        dom.AC.enable_voltage_sweep = getBool(jAC, "enable_voltage_sweep", false, "fields.AC");
        dom.AC.amplitude_slope_V_s = getDouble(jAC, "amplitude_slope_V_s", 0.0, "fields.AC");
        dom.AC.start_time_s = getDouble(jAC, "start_time_s", 0.0, "fields.AC");
        dom.AC.rise_time_s = getDouble(jAC, "rise_time_s", 0.0, "fields.AC");
        dom.AC.enable_frequency_sweep = getBool(jAC, "enable_frequency_sweep", false, "fields.AC");
        dom.AC.ac_start_freq_Hz = getDouble(jAC, "start_freq_Hz", dom.AC.frequency_Hz, "fields.AC");
        dom.AC.ac_sweep_slope_Hz_per_s = getDouble(jAC, "sweep_slope_Hz_per_s", 0.0, "fields.AC");

        // Voltage time table
        if (jAC.isMember("voltage_time_table") && jAC["voltage_time_table"].isArray()) {
            const auto& tab = jAC["voltage_time_table"];
            for (const auto& entry : tab) {
                if (entry.isArray() && entry.size() >= 2) {
                    double t = entry[0].asDouble();
                    double v = entry[1].asDouble();
                    dom.AC.voltage_time_table.emplace_back(t, v);
                }
            }
            std::sort(dom.AC.voltage_time_table.begin(), dom.AC.voltage_time_table.end(), 
                     [](auto &a, auto &b){ return a.first < b.first; });
        }
    }

    // Magnetic field
    if (fieldsSection && fieldsSection->isMember("B")) {
        const auto& jB = (*fieldsSection)["B"];
        dom.B.enabled = getBool(jB, "enabled_magnetic_field", false, "fields.B");
        
        if (jB.isMember("field_strength_T") && jB["field_strength_T"].isArray() && jB["field_strength_T"].size() == 3) {
            dom.B.field_strength_T = Vec3{
                jB["field_strength_T"][0].asDouble(),
                jB["field_strength_T"][1].asDouble(),
                jB["field_strength_T"][2].asDouble()
            };
        } else {
            dom.B.field_strength_T = Vec3{0.0, 0.0, 0.0};
        }
        
        if (jB.isMember("field_gradient_T_m") && jB["field_gradient_T_m"].isArray() && jB["field_gradient_T_m"].size() == 3) {
            dom.B.field_gradient_T_m = Vec3{
                jB["field_gradient_T_m"][0].asDouble(),
                jB["field_gradient_T_m"][1].asDouble(),
                jB["field_gradient_T_m"][2].asDouble()
            };
        } else {
            dom.B.field_gradient_T_m = Vec3{0.0, 0.0, 0.0};
        }
    }

    // --- Solver type (default: RK4) ---
    dom.solver_type = SolverType::RK4;
    if (j.isMember("integrator")) {
        std::string solver = j["integrator"].asString();
        std::transform(solver.begin(), solver.end(), solver.begin(), 
                      [](unsigned char c){ return std::tolower(c); });
        if (solver == "rk45") dom.solver_type = SolverType::RK45;
        else if (solver == "rk4") dom.solver_type = SolverType::RK4;
        else if (solver == "boris") dom.solver_type = SolverType::Boris;
        else {
            throw std::runtime_error("Unknown integrator: '" + solver + "'. Valid options: RK4, RK45, Boris");
        }
    }

    // --- Precomputed field arrays ---
    dom.FA_file.clear();
    dom.FA_terms.clear();
    
    if (fieldsSection && fieldsSection->isMember("field_array")) {
        const auto& fa = (*fieldsSection)["field_array"];
        if (fa.isString()) {
            // Legacy single file
            dom.FA_file = fa.asString();
        } else if (fa.isObject() && fa.isMember("file")) {
            dom.FA_file = fa["file"].asString();
        } else if (fa.isArray()) {
            // Superposition mode
            for (const auto& term : fa) {
                InstrumentDomain::FieldArrayTerm t;
                t.file = getString(term, "file", "", "field_array");
                
                std::string kind = getString(term, "scale_type", "constant", "field_array");
                std::transform(kind.begin(), kind.end(), kind.begin(), 
                              [](unsigned char c){ return std::tolower(c); });
                
                if (kind == "dc_axial") t.kind = InstrumentDomain::FAScaleKind::DC_Axial;
                else if (kind == "dc_quad") t.kind = InstrumentDomain::FAScaleKind::DC_Quad;
                else if (kind == "dc_radial") t.kind = InstrumentDomain::FAScaleKind::DC_Radial;
                else if (kind == "rf") t.kind = InstrumentDomain::FAScaleKind::RF;
                else t.kind = InstrumentDomain::FAScaleKind::Constant;
                
                t.constant = getDouble(term, "constant_V", 1.0, "field_array");
                t.phase_rad = getDouble(term, "phase_rad", 0.0, "field_array");
                t.frequency_Hz = getDouble(term, "frequency_Hz", 0.0, "field_array");
                
                if (!t.file.empty()) {
                    dom.FA_terms.push_back(std::move(t));
                }
            }
        }
    }

    // Legacy FA_superposition support
    if (j.isMember("FA_superposition") && j["FA_superposition"].isArray()) {
        const auto& arr = j["FA_superposition"];
        for (const auto& term : arr) {
            InstrumentDomain::FieldArrayTerm t;
            t.file = term.get("file", "").asString();
            std::string kind = term.get("scale", "constant").asString();
            std::transform(kind.begin(), kind.end(), kind.begin(), [](unsigned char c){ return std::tolower(c); });
            if (kind == "dc_axial") t.kind = InstrumentDomain::FAScaleKind::DC_Axial;
            else if (kind == "dc_quad") t.kind = InstrumentDomain::FAScaleKind::DC_Quad;
            else if (kind == "dc_radial") t.kind = InstrumentDomain::FAScaleKind::DC_Radial;
            else if (kind == "rf") t.kind = InstrumentDomain::FAScaleKind::RF;
            else t.kind = InstrumentDomain::FAScaleKind::Constant;
            t.constant = term.get("constant_V", 1.0).asDouble();
            t.phase_rad = term.get("phase_rad", 0.0).asDouble();
            t.frequency_Hz = term.get("frequency_Hz", 0.0).asDouble();
            dom.FA_terms.push_back(std::move(t));
        }
    } else if (j.isMember("FA_file")) {
        // legacy single file path
        dom.FA_file = j["FA_file"].asString();
    }

    // --- Axis transformation ---
    if (j.isMember("axes")) {
    Vec3 x = Vec3{ j["axes"]["x"][0].asDouble(),
                   j["axes"]["x"][1].asDouble(),
                   j["axes"]["x"][2].asDouble() };

    Vec3 y = Vec3{ j["axes"]["y"][0].asDouble(),
                   j["axes"]["y"][1].asDouble(),
                   j["axes"]["y"][2].asDouble() };

    Vec3 z = Vec3{ j["axes"]["z"][0].asDouble(),
                   j["axes"]["z"][1].asDouble(),
                   j["axes"]["z"][2].asDouble() };

    dom.rotation_local_to_global = Mat3::fromColumns(x, y, z);
    dom.rotation_global_to_local = transpose(dom.rotation_local_to_global);
    } else if (dom.instrument == Instrument::TOF) {
        // for TOF: standard DC field acts on z axis, however x is defined as flight axis -> transformation required
        dom.rotation_global_to_local = Mat3::fromColumns(
            Vec3{0, 0, 1},  // global z -> local x
            Vec3{0, 1, 0},  // global y -> local y 
            Vec3{1, 0, 0}   // global x -> local z
        );
        dom.rotation_local_to_global = transpose(dom.rotation_global_to_local);
    } else {
        // for all instruments except for TOF: no transformation in standard configuration
        dom.rotation_global_to_local = Mat3::identity();
        dom.rotation_local_to_global = Mat3::identity();
    }

    // --- Instrument parameter sanity check ---
    sanity_check_domain(dom);

    return dom;
}

/**
 * @brief Checks global simulation parameters for validity.
 *
 * Ensures general simulation settings (like number of ions) are physically reasonable.
 * Throws std::runtime_error if any parameter is invalid.
 *
 * @param[in] gParams Global simulation parameters.
 * @throws std::runtime_error if a parameter is unphysical.
 */
void run_guard_check_global(const GlobalParams& gParams) {
    if (gParams.dt_s <= 0) throw std::runtime_error("Timestep must be positive and non-zero.");
    if (gParams.write_interval <= 0) throw std::runtime_error("Writer interval must be positive.");
    if (gParams.sim_time_steps <= 0) throw std::runtime_error("Number of simulation steps must be positive and non-zero.");
    if (gParams.num_ions < 0) throw std::runtime_error("Number of ions must be positive.");
}

/**
 * @brief Checks instrument-domain-specific parameters for validity.
 *
 * Validates that the instrument length, radius, environmental conditions,
 * and Orbitrap-specific radii are physically reasonable. Uses the global
 * collision model if needed for density checks.
 *
 * @param[in] dom Instrument domain to check.
 * @param[in] cm Collision model used in the simulation.
 * @throws std::runtime_error if any parameter is unphysical.
 */
void run_guard_check_domain(const InstrumentDomain& dom, CollisionModel cm) {
    if (dom.geom.length_m <= 0.0)
        throw std::runtime_error("Invalid instrument parameter: Length must be positive.");
    if (dom.geom.radius_m <= 0.0)
        throw std::runtime_error("Invalid instrument parameter: Radius must be positive.");
    if (dom.env.temperature_K <= 0.0)
        throw std::runtime_error("Invalid environment parameter: Temperature must be positive.");
    if (dom.env.pressure_Pa <= 0.0)
        throw std::runtime_error("Invalid environment parameter: Pressure must be positive.");
    if (cm == CollisionModel::Friction && dom.env.particle_density_m_3 <= 0.0)
        throw std::runtime_error("Invalid environment parameter: Particle density must be positive.");
    if (dom.instrument == Instrument::Orbitrap &&
        (dom.geom.radius_char_m <= 0.0 || dom.geom.radius_in_m <= 0.0 ||
         dom.geom.radius_out_m <= 0.0))
        throw std::runtime_error("Orbitrap radii must be provided for simulation.");
}

/**
 * @brief Helper function to determine the radial boundaries for given z position.
 */

double orbitrap_r_for_z(double z, double R, double Rm) {
    // find r such that z^2 = 0.5*(r^2 - R^2) - Rm^2 * ln(r/R)
    const double z2 = z * z;
    const double eps = 1e-12;
    const double r_min = 1e-6;     // avoid singularity
    const double r_max = 0.2;      // 20 cm, arbitrary safe bound
    const int max_iter = 100;
    double r_lo = R;       // start around electrode radius
    double r_hi = R * 2.0; // reasonable upper bound
    double r_mid = R;

    for (int i = 0; i < max_iter; ++i) {
        r_mid = 0.5 * (r_lo + r_hi);
        double f = 0.5 * (r_mid * r_mid - R * R) + Rm * Rm * std::log(R / r_mid) - z2;

        if (std::fabs(f) < eps) break;
        if (f > 0) r_lo = r_mid;
        else       r_hi = r_mid;
    }

    return r_mid;
}

/**
 * @brief Checks whether an ion position (global) lies inside a given instrument domain.
 */
bool isInsideDomain(const InstrumentDomain& dom, const Vec3& globalPos) {
    Vec3 local = dom.rotation_global_to_local * (globalPos - dom.geom.origin_m);
    double r = std::sqrt(local.x * local.x + local.y * local.y);
    
    // Tolerance to catch ions at or beyond boundary
    // Use strict < to ensure ions AT the boundary are marked as outside
    const double tol = 0.0;  // No tolerance - strict boundary check
    
    // simple cylindrical geometry for all instruments except for Orbitrap
    if (dom.instrument != Instrument::Orbitrap) {
        // Comparison: ion must be inside domain. Allow ions at the lower
        // z boundary (local.z == 0) to be considered inside because many
        // testcases and geometry conventions place the origin at the
        // entrance plane. Keep upper bound strict so ions exactly at
        // z == length_m are considered outside.
        return (local.z >= tol && local.z < dom.geom.length_m) && (r < dom.geom.radius_m);
    }

    // for Orbitrap: more complex hyperbolic limits
    const double Rin = dom.geom.radius_in_m;
    const double Rout = dom.geom.radius_out_m;
    const double Rm = dom.geom.radius_char_m;
    const double z = std::fabs(local.z);

    // compute allowed radial range for given z
    const double r_in_allowed  = orbitrap_r_for_z(z, Rin,  Rm);
    const double r_out_allowed = orbitrap_r_for_z(z, Rout, Rm);

    // inside domain if between hyperbolic surfaces
    return (r >= r_in_allowed + tol) && (r <= r_out_allowed - tol);
}

/**
 * @brief Looks for the current domain index domain.
 * @return Domain index if domain is found, else -1.
 */
int find_domain_index(const Vec3& pos, const std::vector<InstrumentDomain>& domains) {
    for (size_t i = 0; i < domains.size(); ++i) {
        if (isInsideDomain(domains[i], pos)) return static_cast<int>(i);
    }
    return -1; // fallback / outside any domain
}

/**
 * @brief Finds the active instrument domain for the given ion position.
 * @return Pointer to active domain or nullptr if outside all domains.
 */
const InstrumentDomain* findInstrumentDomain(const std::vector<InstrumentDomain>& domains,
                                             const Vec3& pos) {
    for (const auto& d : domains)
        if (isInsideDomain(d, pos))
            return &d;
    return nullptr;
}

}  // namespace core
}  // namespace ICARION
// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   ------------------------------------------------
 *   Modular framework for simulating ion trajectories in custom
 *   electric fields and background gas environments.
 *
 *   @file       SimulationEngine.cpp
 *   @brief      Core simulation engine implementation
 *
 *   @details
 *   Implements the SimulationEngine class that orchestrates the simulation workflow,
 *   including initialization, time-stepping, force calculations,
 *   integrator management, output handling, and performance monitoring.
 *
 *   @date       2025-11-09
 *   @version    1.0.0
 *   @authors    ICARION Development Team
 *
 * =====================================================================
 */

#include "simulation/SimulationEngine.h"
#include "core/io/hdf5Writer.h"
#include "core/constants/constants.h"
#include "physics/forces.h"
#include "instrument/InstrumentFactory.h"
#include "core/paramUtils/paramUtils.h"
#ifdef USE_GPU_ACCEL
#include "gpuUtils/GpuIntegrator.h"
#endif
#include "core/debug/Debug.h"
#include "core/io/ConfigNormalizer.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <cctype>
#include <cstdlib>
#include <initializer_list>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <unordered_set>
#include <utility>

// CPU stochastic collision and reaction support
#include "physics/collisions/collisionHelpers.h"
#include "physics/reactions/reactionUtils.h"
#include "core/constants/constants.h"
#include "physics/geometryReader.h"

namespace ICARION {
namespace simulation {

namespace {

// ========================================================================
// CPU STOCHASTIC COLLISION AND REACTION HELPERS
// Restored from legacy integrator.cpp - implements EXACT physics
// ========================================================================

/**
 * @brief Apply stochastic collision to a single ion (ONCE per timestep)
 * 
 * Restored legacy logic from integrator.cpp handle_collision().
 * Collision probability: P = 1 - exp(-n * sigma_eff * v_rel * dt)
 */
static inline void cpu_handle_stochastic_collision(
    core::IonState& ion,
    std::mt19937_64& rng,
    double dt,
    CollisionModel collision_model,
    const std::unordered_map<std::string, std::pair<std::vector<Vec3>, std::vector<double>>>& geometry_map)
{
    // Skip if ion inactive or deterministic collision model
    if (!ion.active || collision_model == CollisionModel::Friction) {
        return;
    }
    
    // Build EHSS parameters from ion's current domain environment
    EHSSParams ep{};
    ep.n = ion.domain_particle_density_m3;
    ep.dt = dt;
    ep.mi = ion.mass_kg;
    ep.mn = ion.domain_neutral_mass_kg;
    ep.kB = BOLTZMANN_CONSTANT;
    ep.Tn = ion.domain_temperature_K;
    ep.sigma_eff = ion.CCS_m2;
    ep.ubx = ion.domain_gas_velocity_m_s.x;
    ep.uby = ion.domain_gas_velocity_m_s.y;
    ep.ubz = ion.domain_gas_velocity_m_s.z;
    
    // Create RNG wrapper for collision helper functions
    EhssRng ehss_rng(rng());
    
    // Sample neutral velocity from Maxwell-Boltzmann distribution
    Vec3 v_neutral = sample_neutral_velocity(ep, ehss_rng);
    Vec3 v_ion = ion.vel;
    
    // Compute relative velocity magnitude
    double v_rel = norm(v_ion - v_neutral);
    
    // Collision probability: P = 1 - exp(-n * sigma_eff * v_rel * dt)
    double P = 1.0 - std::exp(-ep.n * ep.sigma_eff * v_rel * ep.dt);
    
    // Monte Carlo acceptance
    if (ehss_rng.uniform01() < P) {
        if (collision_model == CollisionModel::EHSS) {
            // EHSS collision with molecular geometry
            auto it = geometry_map.find(ion.species_id);
            if (it != geometry_map.end() && !it->second.first.empty() && !it->second.second.empty()) {
                const auto& h_centers = it->second.first;
                const auto& h_radii = it->second.second;
                ion.vel = collide_ehss_cpu_geometry_given_neutral(v_ion, v_neutral, ep, h_centers, h_radii, ehss_rng);
            } else {
                // Fallback to HSMC if geometry unavailable
                ion.vel = collide_hs_cpu(v_ion, v_neutral, ep, ehss_rng);
            }
        } else if (collision_model == CollisionModel::HSMC) {
            // Hard-sphere Monte Carlo collision
            ion.vel = collide_hs_cpu(v_ion, v_neutral, ep, ehss_rng);
        }
    }
}

/**
 * @brief Apply stochastic reactions to a single ion (ONCE per timestep)
 * 
 * Restored legacy logic from integrator.cpp handle_reaction().
 * Reaction probability: P = 1 - exp(-k_eff * dt)
 */
static inline void cpu_handle_stochastic_reaction(
    core::IonState& ion,
    std::mt19937_64& rng,
    double dt,
    bool enable_reactions,
    const std::unordered_map<std::string, Species>& species_db,
    const std::vector<ReactionEntry>& reaction_list)
{
    if (!enable_reactions || !ion.active) return;
    
    const std::string current_species = ion.species_id;
    
    // Create RNG wrapper
    EhssRng ehss_rng(rng());
    
    for (const auto& rxn : reaction_list) {
        if (rxn.reactant != current_species) continue;
        
        // Calculate effective rate constant k_eff [m³/s] -> [s⁻¹]
        double k_eff = rxn.rate_constant;  // [m³/s]
        double neutral_density_m3 = ion.domain_particle_density_m3;  // [m⁻³]
        
        // Apply concentration dependence for reaction order
        for (const auto& term : rxn.order) {
            // Detect if neutral_concentration is volume fraction or absolute density
            double conc_density = 0.0;  // [m⁻³]
            if (rxn.neutral_concentration > 1e6) {
                conc_density = rxn.neutral_concentration;  // already number density
            } else {
                conc_density = rxn.neutral_concentration * neutral_density_m3;  // volume fraction
            }
            if (conc_density < 0.0) conc_density = 0.0;
            k_eff *= std::pow(conc_density, term.exponent);
        }
        
        // Reaction probability: P = 1 - exp(-k_eff * dt)
        double P = 1.0 - std::exp(-k_eff * dt);
        
        if (ehss_rng.uniform01() < P) {
            // Find product species in database
            auto it = species_db.find(rxn.product);
            if (it == species_db.end()) {
                std::cerr << "Warning: Product species '" << rxn.product << "' not found in database\\n";
                continue;
            }
            
            const auto& prod = it->second;
            ion.species_id = rxn.product;
            ion.mass_kg = prod.mass_kg;
            ion.reduced_mobility_cm2_Vs = prod.mobility * 1e4;  // m²/Vs to cm²/Vs
            ion.ion_charge_C = prod.charge * ELEM_CHARGE_C;
            
            // Recalculate CCS from mobility using Mason-Schamp relation
            ion.CCS_m2 = 3.0 / 16.0 / LOSCHMIDT_CONSTANT * ion.ion_charge_C *
                       std::sqrt(2.0 * M_PI /
                                 (BOLTZMANN_CONSTANT * ion.domain_temperature_K *
                                  (ion.mass_kg * ion.domain_neutral_mass_kg) /
                                  (ion.mass_kg + ion.domain_neutral_mass_kg))) *
                       1.0 / (ion.reduced_mobility_cm2_Vs * 1e-4);
            
            // Reaction successful - stop processing further reactions
            return;
        }
    }
}

constexpr double kDefaultMobilityCm2Vs = 1.0;
constexpr double kDefaultCcsM2 = 1e-18;

// Helper: Sample velocity component from Maxwell-Boltzmann distribution
// v ~ sqrt(k_B * T / m) * N(0,1)
double sampleMaxwellBoltzmannComponent(std::mt19937& rng, double temperature_K, double mass_kg) {
    if (temperature_K <= 0.0 || mass_kg <= 0.0) return 0.0;
    
    std::normal_distribution<double> normal(0.0, 1.0);
    double sigma_v = std::sqrt(BOLTZMANN_CONSTANT * temperature_K / mass_kg);
    return sigma_v * normal(rng);
}

// Helper: Sample position component from Gaussian distribution
double sampleGaussian(std::mt19937& rng, double center, double sigma) {
    std::normal_distribution<double> normal(center, sigma);
    return normal(rng);
}

// Helper to get simulation parameters from v1.0 format
const Json::Value& getSimParams(const Json::Value& cfg) {
    static const Json::Value kEmpty(Json::objectValue);
    if (cfg.isMember("simulation") && cfg["simulation"].isObject()) {
        return cfg["simulation"];
    }
    return kEmpty;
}

double lookupNeutralMass(const std::string& species_id);

std::string resolveReactionDatabasePath(const Json::Value& cfg) {
    const Json::Value& sim = getSimParams(cfg);
    if (sim.isMember("reaction_database") && sim["reaction_database"].isString()) {
        return sim["reaction_database"].asString();
    }
    if (cfg.isMember("reaction_database") && cfg["reaction_database"].isString()) {
        return cfg["reaction_database"].asString();
    }
    if (cfg.isMember("collisions") && cfg["collisions"].isObject()) {
        const auto& coll = cfg["collisions"];
        if (coll.isMember("reaction_database") && coll["reaction_database"].isString()) {
            return coll["reaction_database"].asString();
        }
    }
    return {};
}

std::string resolveGeometryFilePath(const Json::Value& cfg) {
    if (cfg.isMember("geometry_file") && cfg["geometry_file"].isString()) {
        return cfg["geometry_file"].asString();
    }
    if (cfg.isMember("instrument") && cfg["instrument"].isObject()) {
        const auto& inst = cfg["instrument"];
        if (inst.isMember("geometry_file") && inst["geometry_file"].isString()) {
            return inst["geometry_file"].asString();
        }
    }
    if (cfg.isMember("collisions") && cfg["collisions"].isObject()) {
        const auto& coll = cfg["collisions"];
        if (coll.isMember("geometry_file") && coll["geometry_file"].isString()) {
            return coll["geometry_file"].asString();
        }
    }
    return {};
}

double resolveEnvironmentTemperature(const Json::Value& cfg) {
    if (cfg.isMember("environment") && cfg["environment"].isObject()) {
        const auto& env = cfg["environment"];
        if (env.isMember("temperature_K") && env["temperature_K"].isNumeric()) {
            return env["temperature_K"].asDouble();
        }
    }
    return 300.0;
}

double resolveNeutralMassKg(const Json::Value& cfg) {
    if (cfg.isMember("environment") && cfg["environment"].isObject()) {
        const auto& env = cfg["environment"];
        if (env.isMember("gas_mass_kg") && env["gas_mass_kg"].isNumeric()) {
            return env["gas_mass_kg"].asDouble();
        }
        if (env.isMember("gas_mass_amu") && env["gas_mass_amu"].isNumeric()) {
            return env["gas_mass_amu"].asDouble() * AMU_TO_KG;
        }
        if (env.isMember("gas_species") && env["gas_species"].isString()) {
            double m = lookupNeutralMass(env["gas_species"].asString());
            if (!std::isnan(m)) {
                return m;
            }
        }
    }
    return MOLAR_MASS_N2_KG;
}

bool loadEhssGeometryForSpecies(const std::string& geometry_file,
                                const std::string& species_id,
                                std::pair<std::vector<Vec3>, std::vector<double>>& out) {
    if (geometry_file.empty() || species_id.empty()) {
        return false;
    }
    try {
        auto molecules = read_geometry_file(geometry_file, species_id);
        if (molecules.empty()) {
            return false;
        }
        const auto& atoms = molecules[0].atoms;
        out.first.clear();
        out.second.clear();
        out.first.reserve(atoms.size());
        out.second.reserve(atoms.size());
        for (const auto& atom : atoms) {
            out.first.push_back({atom.posx_m, atom.posy_m, atom.posz_m});
            out.second.push_back(0.5 * atom.LJ_sigma_m);
        }
        return !atoms.empty();
    } catch (const std::exception& ex) {
        std::cerr << "[SimulationEngine] EHSS geometry load failed for '" << species_id
                  << "' from '" << geometry_file << "': " << ex.what() << std::endl;
        return false;
    }
}

bool gpuDebugEnabled() {
    static bool enabled = [] {
        const char* env = std::getenv("ICARION_DEBUG_GPU");
        if (!env || env[0] == '\0') return false;
        std::string s(env);
        for (auto &c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (s == "1" || s == "true" || s == "yes" || s == "on") return true;
        // allow numeric > 0
        char* endptr = nullptr;
        long v = std::strtol(env, &endptr, 10);
        if (endptr != env) return v > 0;
        return false;
    }();
    return enabled;
}

double readIonParameter(const Json::Value& node,
                        std::initializer_list<const char*> keys) {
    for (const char* key : keys) {
        if (node.isMember(key) && node[key].isNumeric()) {
            return node[key].asDouble();
        }
    }
    return std::numeric_limits<double>::quiet_NaN();
}

std::string normalizeSpeciesId(const std::string& raw) {
    std::string result;
    result.reserve(raw.size());
    for (char c : raw) {
           if (std::isalpha(static_cast<unsigned char>(c))) {
            result.push_back(static_cast<char>(std::tolower(static_cast<unsigned char>(c))));
        }
    }
    return result;
}

double lookupNeutralMass(const std::string& species_id) {
    const std::string key = normalizeSpeciesId(species_id);
    if (key == "he" || key == "helium") {
        return MOLAR_MASS_HE_KG;
    } else if (key == "ne" || key == "neon") {
        return MOLAR_MASS_NE_KG;
    } else if (key == "ar" || key == "argon") {
        return MOLAR_MASS_AR_KG;
    } else if (key == "n2" || key == "nitrogen") {
        return MOLAR_MASS_N2_KG;
    } else if (key == "o2" || key == "oxygen") {
        return MOLAR_MASS_O2_KG;
    }
    return std::numeric_limits<double>::quiet_NaN();
}

double lookupPolarizability(const std::string& species_id) {
    const std::string key = normalizeSpeciesId(species_id);
    if (key == "he" || key == "helium") {
        return POLARIZABILITY_HE_SI;
    } else if (key == "ne" || key == "neon") {
        return POLARIZABILITY_NE_SI;
    } else if (key == "ar" || key == "argon") {
        return POLARIZABILITY_AR_SI;
    } else if (key == "n2" || key == "nitrogen") {
        return POLARIZABILITY_N2_SI;
    } else if (key == "o2" || key == "oxygen") {
        return POLARIZABILITY_O2_SI;
    }
}
double resolvePolarizability(const Json::Value& ion_node) {
    double polarizability = readIonParameter(
        ion_node,
        {"polarizability_m3", "polarizability"});
    if (polarizability > 0.0) {
        return polarizability = lookupPolarizability(ion_node.get("gas_species", "n2").asString());
    }
 
    return POLARIZABILITY_N2_SI;  // default to N2
}
}  // namespace

SimulationEngine::SimulationEngine() 
    : instrument_type_(instrument::InstrumentType::UnknownInstrument),
      is_multi_domain_(false),
      is_initialized_(false),
      is_running_(false),
      stop_requested_(false),
      current_step_(0),
      current_time_(0.0),
      execution_override_(ExecutionMode::Auto),
      gpu_mode_hint_(false),
    gpu_validation_enabled_(false),
      gpu_validation_pos_tol_(1e-3),
      gpu_validation_vel_tol_(1e-3),
      has_seed_override_(false),
      seed_override_value_(0),
      output_initialized_(false),
      output_finalized_(false),
      accumulated_step_time_(0.0),
      performance_sample_count_(0),
      initial_energy_(0.0),
      previous_energy_(0.0),
      safety_violation_count_(0),
      cpu_stochastic_collisions_enabled_(false),
      cpu_reactions_enabled_(false),
      cpu_collision_model_(CollisionModel::NoCollisions) {}

SimulationEngine::~SimulationEngine() = default;

void SimulationEngine::setExecutionMode(ExecutionMode mode) {
    execution_override_ = mode;
}

void SimulationEngine::setRandomSeed(unsigned long seed) {
    has_seed_override_ = true;
    seed_override_value_ = seed;
}

unsigned long SimulationEngine::getRandomSeed() const {
    return config_.rng_seed;
}

void SimulationEngine::initialize(const Json::Value& config) {
    Json::Value normalized;
    try {
        normalized = ConfigNormalizer::normalize(config);
    } catch (const std::exception& ex) {
        Json::StreamWriterBuilder b; b["indentation"] = "";
        std::cerr << "[Engine DEBUG] ConfigNormalizer::normalize threw: " << ex.what() << std::endl;
        std::cerr << "[Engine DEBUG] Input config: " << Json::writeString(b, config) << std::endl;
        throw;
    }
    const Json::Value& cfg = normalized;
    parseConfiguration(cfg);
    if (!cfg.isMember("ions")) {
        throw std::runtime_error("Simulation configuration missing required 'ions' section (v1.0 schema)");
    }
    const Json::Value& ion_block = cfg["ions"];
    setupIons(ion_block);
#if 1
    // Decide whether CPU integrator must be created: only create it if the
    // simulation will run on CPU or if the user explicitly requested GPU parity
    // validation. Otherwise avoid instantiating the CPU integrator to save
    // memory/time when running GPU-only.
    const Json::Value& sim_params = getSimParams(cfg);
    bool gpu_execution_requested = config_.enable_gpu;
    if (execution_override_ == ExecutionMode::ForceGPU) {
        gpu_execution_requested = true;
    } else if (execution_override_ == ExecutionMode::ForceCPU) {
        gpu_execution_requested = false;
    }

    // Read debug.validate_gpu (opt-in). Default: false (no automatic validation)
    const Json::Value& debug_block = cfg.isMember("debug") ? cfg["debug"] : Json::Value(Json::objectValue);
    bool validate_gpu_flag = debug_block.get("validate_gpu", false).asBool();
    // Allow an environment variable to enable validation during CI/debug runs
    // without changing test sources. ICARION_DEBUG_VALIDATE_GPU accepts truthy
    // values (1/true/yes/on or numeric > 0).
    const char* env_validate = std::getenv("ICARION_DEBUG_VALIDATE_GPU");
    if (env_validate && env_validate[0] != '\0') {
        std::string s(env_validate);
        for (auto &c : s) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
        if (s == "1" || s == "true" || s == "yes" || s == "on") {
            validate_gpu_flag = true;
        } else {
            char* endptr = nullptr;
            long v = std::strtol(env_validate, &endptr, 10);
            if (endptr != env_validate) validate_gpu_flag = (v > 0);
        }
    }

    // Create CPU integrator only if CPU is requested or parity validation explicitly enabled
    bool need_cpu_integrator = (!gpu_execution_requested) || validate_gpu_flag;
    setupIntegrator(cfg, need_cpu_integrator);

    setupInstrument(cfg);
#else
    setupIntegrator(cfg);
    setupInstrument(cfg);
#endif
#ifdef USE_GPU_ACCEL
    // `gpu_execution_requested` was determined above.
    gpu_integrator_.reset();
    if (gpu_execution_requested) {
        gpu_integrator_ = std::make_unique<gpu::GpuIntegrator>(cfg);
        
        // Override GPU domains based on configuration type
        if (gpu_integrator_) {
            if (is_multi_domain_ && !instrument_domains_.empty()) {
                // Multi-domain: use all collected instrument domains
                std::cerr << "[SimulationEngine] Multi-domain GPU: overriding " << instrument_domains_.size() 
                          << " domains" << std::endl;
                gpu_integrator_->overrideDomains(instrument_domains_);
            }
            // Note: Single-domain case deferred to setupForces() to ensure instrument is fully configured
        }
        // Read gpu_download_interval from config and propagate to integrator.
        const Json::Value& sim_params2 = getSimParams(cfg);
        // By default, align device->host download frequency with output writes so
        // the CPU only needs device state when writing output. If the user has
        // explicitly configured "gpu_download_interval" that value takes
        // precedence. The environment variable ICARION_GPU_D2H_PERIOD still
        // overrides both and will be honored by the integrator constructor.
        int default_interval = sim_params2.get("output_frequency", 1000).asInt();
        if (default_interval < 1) default_interval = 1000; // fallback safe default
        int download_interval = 0;
        if (sim_params2.isMember("gpu_download_interval")) {
            download_interval = sim_params2.get("gpu_download_interval", 1).asInt();
        } else {
            download_interval = default_interval;
        }

        // Set gpu validation flag from debug.validate_gpu (opt-in). Default false.
        gpu_validation_enabled_ = validate_gpu_flag;
        
        // Auto-adjust download interval for validation: immediate sync needed for parity testing
        if (gpu_validation_enabled_) {
            download_interval = 1;  // Override any config for immediate device-to-host sync
            std::cout << "[Validation] GPU download interval auto-set to 1 for parity testing" << std::endl;
            std::cout << "[Validation] GPU parity checking ENABLED (immediate sync mode)" << std::endl;
        } else {
            std::cout << "[Validation] GPU parity checking disabled (performance mode)" << std::endl;
        }
        
        // The integrator will honor the ICARION_GPU_D2H_PERIOD env var if present.
        if (gpu_integrator_) {
            gpu_integrator_->setDownloadInterval(download_interval);
        }
    }
#endif
    setupForces(cfg);
    loadCpuStochasticAssets(cfg);
    // If parity dump env is set, avoid creating HDF5 output files to prevent
    // file-access conflicts when running CPU and GPU runs in the same process.
    const char* parity_dbg_env_init = std::getenv("ICARION_DEBUG_PARITY_DUMP");
    if (parity_dbg_env_init && parity_dbg_env_init[0] != '\0') {
        std::cout << "[SimulationEngine] Parity dump active: skipping HDF5 output initialization" << std::endl;
        output_initialized_ = false;
        output_finalized_ = false;
        hdf5_writer_.reset();
    } else {
        setupOutput(cfg);
    }

    if (has_seed_override_) {
        config_.rng_seed = seed_override_value_;
    }
    if (config_.rng_seed == 0) {
        config_.rng_seed = static_cast<unsigned long>(std::time(nullptr));
    }
    
    updateOutputMetadata();

#ifdef USE_GPU_ACCEL
    bool gpu_requested = config_.enable_gpu;
    if (execution_override_ == ExecutionMode::ForceGPU) {
        gpu_requested = true;
    } else if (execution_override_ == ExecutionMode::ForceCPU) {
        gpu_requested = false;
    }
    if (gpuDebugEnabled()) {
        std::cerr << "[SimulationEngine] gpu_requested=" << (gpu_requested ? "true" : "false")
                  << " gpu_integrator_=" << (gpu_integrator_ ? "non-null" : "null");
        if (gpu_integrator_) {
            std::cerr << " gpu_integrator_->isAvailable()=" << (gpu_integrator_->isAvailable() ? "true" : "false");
        }
        std::cerr << std::endl;
    }
    std::cerr << "[SimulationEngine] Setting gpu_mode_hint_: gpu_requested=" << (gpu_requested ? "true" : "false")
              << ", gpu_integrator_=" << (gpu_integrator_ ? "non-null" : "null");
    if (gpu_integrator_) {
        std::cerr << ", gpu_integrator_->isAvailable()=" << (gpu_integrator_->isAvailable() ? "true" : "false");
    }
    std::cerr << std::endl;
    gpu_mode_hint_ = gpu_requested && gpu_integrator_ && gpu_integrator_->isAvailable();
    std::cerr << "[SimulationEngine] gpu_mode_hint_ = " << (gpu_mode_hint_ ? "true" : "false") << std::endl;
    
    // Now that gpu_mode_hint_ is set, configure GPU domains if needed
    if (gpu_mode_hint_ && gpu_integrator_) {
        if (instrument_domains_.empty() && instrument_) {
            instrument_domains_ = instrument_->buildGpuDomains();
        }
        if (!instrument_domains_.empty()) {
            std::cerr << "[SimulationEngine] Post-init: Calling gpu_integrator_->overrideDomains with " 
                      << instrument_domains_.size() << " domains" << std::endl;
            gpu_integrator_->overrideDomains(instrument_domains_);
        } else {
            std::cerr << "[SimulationEngine] WARNING: GPU mode enabled but no instrument domains available" << std::endl;
        }
    }
    
    if (gpu_requested && (!gpu_integrator_ || !gpu_integrator_->isAvailable())) {
        std::cout << "GPU execution requested but no CUDA device is available. "
                     "Falling back to CPU integrator."
                  << std::endl;
    }
#else
    gpu_mode_hint_ = false;
#endif
    
    // Calculate initial energy
    initial_energy_ = calculateTotalEnergy();
    previous_energy_ = initial_energy_;
    
    is_initialized_ = true;
}

SimulationResult SimulationEngine::run() {
    if (!is_initialized_) {
        return {false, "Simulation not initialized", 0, 0, 0.0, 0.0, 0, {}};
    }
    
    is_running_ = true;
    stop_requested_ = false;
    current_step_ = 0;
    current_time_ = 0.0;
    
    start_time_ = std::chrono::high_resolution_clock::now();
    
    // Check debug env to optionally print per-step ion states for parity debugging
    const char* parity_dbg_env = std::getenv("ICARION_DEBUG_PARITY_DUMP");
    const bool parity_dump = (parity_dbg_env && parity_dbg_env[0] != '\0');

    // If parity dump is active, disable HDF5 writer here as well in case
    // initialize() already created files — this avoids file-access conflicts
    // when running CPU and GPU runs back-to-back in the same process.
    if (parity_dump) {
        hdf5_writer_.reset();
        output_initialized_ = false;
        output_finalized_ = false;
    }

    try {
        // Main simulation loop
        for (current_step_ = 0; current_step_ < config_.total_steps; ++current_step_) {
            if (stop_requested_) break;
            
            step_start_time_ = std::chrono::high_resolution_clock::now();
            
        try {
            performSimulationStep();
        } catch (const std::exception& ex) {
#ifdef USE_GPU_ACCEL
            if (gpu_mode_hint_ && gpu_integrator_) {
                std::cerr << "GPU step failed: " << ex.what()
                          << ". Disabling GPU and retrying step on CPU." << std::endl;
                gpu_mode_hint_ = false;
                auto timed_force_function = createTimedForceFunction();
                integrator_->step(ions_, config_.timestep, current_time_, timed_force_function);
            } else
#endif
            {
                throw;
            }
        }

        // Optional parity dump: print first N steps of ion[0] state for debugging
        if (parity_dump && current_step_ < 20 && !ions_.empty()) {
            const auto &ion = ions_.front();
            std::cout << "[PARITY-DUMP][step] " << current_step_ << " t=" << current_time_
                      << " pos=(" << ion.pos.x << "," << ion.pos.y << "," << ion.pos.z << ")"
                      << " vel=(" << ion.vel.x << "," << ion.vel.y << "," << ion.vel.z << ")"
                      << std::endl;
        }
            
            updatePerformanceMetrics();
            updateEnergyConservation();
            
            if (current_step_ % config_.output_frequency == 0) {
                writeOutput();
            }
            
            current_time_ += config_.timestep;
        }
        
        // Final output
        writeOutput();
        // Ensure GPU device state is visible on the host before finalizing outputs
        // and returning results. This is important for correctness of tests that
        // inspect engine.ions_ after a GPU run (parity checks, etc.).
#ifdef USE_GPU_ACCEL
        if (gpu_integrator_) {
            try {
                gpu_integrator_->forceDownload();
            } catch (...) {
                // ignore errors here - best-effort synchronization
            }
        }
#endif
        finalizeOutput();
        std::cout << "Simulation completed successfully!" << std::endl;
        
        is_running_ = false;
        
        double avg_ns_per_ion_step = 0.0;
        if (performance_sample_count_ > 0 && !ions_.empty()) {
            avg_ns_per_ion_step =
                (accumulated_step_time_ / performance_sample_count_) / ions_.size();
        }
        double energy_error = 0.0;
        if (initial_energy_ != 0.0) {
            energy_error = std::abs((previous_energy_ - initial_energy_) / initial_energy_);
        }
        
        return {
            true,                           // success
            "",                            // error_message
            current_step_,                 // steps_completed
            config_.total_steps,           // total_steps
            avg_ns_per_ion_step,          // performance_ns_per_ion_step
            energy_error,                  // energy_conservation_error
            safety_violation_count_,       // safety_violations
            {config_.output_file}          // output_files
        };
        
    } catch (const std::exception& e) {
        is_running_ = false;
        finalizeOutput();
        return {false, e.what(), current_step_, config_.total_steps, 0.0, 0.0, safety_violation_count_, {}};
    }
}

SimulationSummary SimulationEngine::getSimulationSummary() const {
    SimulationSummary summary;
    summary.title = config_.title;
    summary.total_ions = ions_.size();
    summary.total_time = config_.total_time;
    summary.timestep = config_.timestep;
    summary.total_steps = config_.total_steps;
    summary.integrator_name = integrator_ ? integrator_->getName() : "None";
    summary.force_names = {force_calculator_ ? force_calculator_->getName() : "None"};
    
#ifdef USE_GPU_ACCEL
    // Ensure any requested summary reflects device state (best-effort).
    if (gpu_integrator_ && gpu_integrator_->isAvailable()) {
        try {
            const_cast<SimulationEngine*>(this)->gpu_integrator_->forceDownload();
        } catch (...) {}
    }
#endif
    
    return summary;
}

const std::vector<core::IonState>& SimulationEngine::getIons() const {
    return ions_;
}

double SimulationEngine::getProgress() const {
    if (config_.total_steps == 0) return 0.0;
    return static_cast<double>(current_step_) / config_.total_steps;
}

void SimulationEngine::parseConfiguration(const Json::Value& config) {
    
    config_.title = config.get("title", "ICARION Simulation").asString();
    
    // Parse simulation parameters (supporting v1.0 format)
    const Json::Value& sim = config["simulation"];
    bool v1_complete = sim.isMember("timestep_ns") && sim.isMember("max_time_ns");
    
    if (v1_complete) {
        double timestep_ns = sim.get("timestep_ns", 0.0).asDouble();
        double max_time_ns = sim.get("max_time_ns", 0.0).asDouble();
        
        // Convert from nanoseconds to seconds
        config_.timestep = timestep_ns * 1e-9;
        config_.total_time = max_time_ns * 1e-9;
        
        config_.output_frequency = sim.get("output_interval", 100).asUInt();
        config_.rng_seed = static_cast<unsigned long>(sim.get("random_seed", static_cast<Json::UInt64>(0)).asUInt64());
        config_.enable_gpu = sim.get("enable_gpu", false).asBool();
        config_.enable_safety_checks = true;
        config_.enable_performance_logging = true;
    } else {
        throw std::runtime_error("Incomplete simulation configuration");
    }
    
    if (config_.timestep <= 0.0) {
        throw std::runtime_error("Simulation timestep must be positive");
    }
    if (config_.total_time <= 0.0) {
        throw std::runtime_error("Simulation total_time must be positive");
    }
    config_.total_steps = static_cast<size_t>(std::ceil(config_.total_time / config_.timestep));
    
    // Parse output file configuration (only v1.0 format supported here)
    if (config.isMember("output") && config["output"].isObject() && 
        config["output"].isMember("trajectory") && config["output"]["trajectory"].isObject()) {
        config_.output_file = config["output"]["trajectory"]
                                  .get("filename", "icarion_output.h5").asString();
    } else {
        config_.output_file = "icarion_output.h5";
        std::cerr << "Warning: No output configuration found, using default output file '"
                  << config_.output_file << "'" << std::endl;
    }
    
    // Parse collision settings and enable CPU stochastic collision support
    cpu_stochastic_collisions_enabled_ = false;
    cpu_collision_model_ = CollisionModel::NoCollisions;
    cpu_reactions_enabled_ = false;
    
    if (config.isMember("collisions") && config["collisions"].isObject()) {
        const Json::Value& coll = config["collisions"];
        
        if (coll.isMember("collision_model") && coll["collision_model"].isString()) {
            std::string model_str = coll["collision_model"].asString();
            
            // Parse collision model for CPU stochastic collision support
            if (model_str == "ehss" || model_str == "EHSS") {
                cpu_collision_model_ = CollisionModel::EHSS;
                cpu_stochastic_collisions_enabled_ = true;
            } else if (model_str == "hsmc" || model_str == "HSMC" || 
                      model_str == "hardsphere" || model_str == "hard_sphere") {
                cpu_collision_model_ = CollisionModel::HSMC;
                cpu_stochastic_collisions_enabled_ = true;
            } else if (model_str == "friction" || model_str == "Friction") {
                cpu_collision_model_ = CollisionModel::Friction;
                cpu_stochastic_collisions_enabled_ = false;
            } else if (model_str == "langevin" || model_str == "Langevin") {
                cpu_collision_model_ = CollisionModel::Langevin;
                cpu_stochastic_collisions_enabled_ = false;
            }
        }
        
        if (coll.isMember("enable_reactions") && coll["enable_reactions"].isBool()) {
            cpu_reactions_enabled_ = coll["enable_reactions"].asBool();
        }

        if (coll.isMember("enable_ou_thermalization") && coll["enable_ou_thermalization"].isBool()) {
            cpu_ou_thermalization_enabled_ = coll["enable_ou_thermalization"].asBool();
        }
    }
}

void SimulationEngine::setupIons(const Json::Value& ion_config) {
    ions_.clear();

    if (!ion_config.isObject() || !ion_config.isMember("species") ||
        !ion_config["species"].isArray() || ion_config["species"].empty()) {
        throw std::runtime_error("ions.species must be a non-empty array (v1.0 schema)");
    }

    const auto& species_array = ion_config["species"];
    const Json::Value* global_distribution =
        (ion_config.isMember("initial_distribution") && ion_config["initial_distribution"].isObject())
            ? &ion_config["initial_distribution"]
            : nullptr;

    auto parseVec3 = [](const Json::Value& node, const Vec3& fallback) -> Vec3 {
        if (node.isArray() && node.size() >= 3) {
            return Vec3(node[0].asDouble(), node[1].asDouble(), node[2].asDouble());
        }
        if (node.isObject()) {
            Vec3 value = fallback;
            if (node.isMember("x") && node["x"].isNumeric()) value.x = node["x"].asDouble();
            if (node.isMember("y") && node["y"].isNumeric()) value.y = node["y"].asDouble();
            if (node.isMember("z") && node["z"].isNumeric()) value.z = node["z"].asDouble();
            return value;
        }
        return fallback;
    };

    auto toLower = [](std::string value) {
        std::transform(value.begin(), value.end(), value.begin(),
                       [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
        return value;
    };

    std::mt19937 rng;
    if (has_seed_override_) {
        rng.seed(seed_override_value_);
    } else if (config_.rng_seed != 0) {
        rng.seed(config_.rng_seed);
    } else {
        rng.seed(static_cast<unsigned long>(std::time(nullptr)));
    }

    ions_.reserve(1024);
    size_t total_generated = 0;

    for (const auto& species_node : species_array) {
        if (!species_node.isObject()) {
            continue;
        }

        const std::string species_id = species_node.get("name", "default").asString();
        const size_t initial_count = species_node.get("initial_count", 0U).asUInt();
        if (initial_count == 0) {
            continue;
        }

        double mass_kg = std::numeric_limits<double>::quiet_NaN();
        if (species_node.isMember("mass_amu") && species_node["mass_amu"].isNumeric()) {
            mass_kg = species_node["mass_amu"].asDouble() * AMU_TO_KG;
        } else if (species_node.isMember("mass_kg") && species_node["mass_kg"].isNumeric()) {
            mass_kg = species_node["mass_kg"].asDouble();
        }
        if (!std::isfinite(mass_kg) || mass_kg <= 0.0) {
            throw std::runtime_error("Species '" + species_id + "' must define mass_amu or mass_kg");
        }

        double ion_charge_C = ELEM_CHARGE_C;
        if (species_node.isMember("charge") && species_node["charge"].isNumeric()) {
            ion_charge_C = species_node["charge"].asDouble() * ELEM_CHARGE_C;
        } else if (species_node.isMember("charge_state") && species_node["charge_state"].isNumeric()) {
            ion_charge_C = species_node["charge_state"].asDouble() * ELEM_CHARGE_C;
        } else if (species_node.isMember("charge_C") && species_node["charge_C"].isNumeric()) {
            ion_charge_C = species_node["charge_C"].asDouble();
        }

        double reduced_mobility = kDefaultMobilityCm2Vs;
        if (species_node.isMember("reduced_mobility_cm2_Vs") && species_node["reduced_mobility_cm2_Vs"].isNumeric()) {
            reduced_mobility = species_node["reduced_mobility_cm2_Vs"].asDouble();
        } else if (species_node.isMember("mobility_cm2_Vs") && species_node["mobility_cm2_Vs"].isNumeric()) {
            reduced_mobility = species_node["mobility_cm2_Vs"].asDouble();
        }

        double ccs_m2 = kDefaultCcsM2;
        if (species_node.isMember("collision_cross_section_A2") && species_node["collision_cross_section_A2"].isNumeric()) {
            ccs_m2 = species_node["collision_cross_section_A2"].asDouble() * 1e-20;
        } else if (species_node.isMember("ccs_m2") && species_node["ccs_m2"].isNumeric()) {
            ccs_m2 = species_node["ccs_m2"].asDouble();
        }

        const Json::Value* position_cfg = nullptr;
        if (species_node.isMember("initial_position")) {
            position_cfg = &species_node["initial_position"];
        } else if (global_distribution && global_distribution->isMember("position")) {
            position_cfg = &(*global_distribution)["position"];
        }

        Vec3 pos_center(0.0, 0.0, 0.0);
        Vec3 pos_sigma(0.001, 0.001, 0.001);
        std::string pos_dist_type = "point";
        if (position_cfg) {
            if (position_cfg->isObject()) {
                pos_dist_type = toLower(position_cfg->get("type", "point").asString());
                if (position_cfg->isMember("center")) {
                    pos_center = parseVec3((*position_cfg)["center"], pos_center);
                }
                if (position_cfg->isMember("sigma")) {
                    pos_sigma = parseVec3((*position_cfg)["sigma"], pos_sigma);
                }
            } else if (position_cfg->isArray()) {
                pos_center = parseVec3(*position_cfg, pos_center);
                pos_dist_type = "point";
            }
        } else if (species_node.isMember("initial_position") && species_node["initial_position"].isArray()) {
            pos_center = parseVec3(species_node["initial_position"], pos_center);
        }

        const Json::Value* velocity_cfg = nullptr;
        if (species_node.isMember("velocity_distribution")) {
            velocity_cfg = &species_node["velocity_distribution"];
        } else if (global_distribution && global_distribution->isMember("velocity")) {
            velocity_cfg = &(*global_distribution)["velocity"];
        }

        if (!velocity_cfg || !velocity_cfg->isObject()) {
            throw std::runtime_error("Species '" + species_id + "' must provide a velocity_distribution block (or ions.initial_distribution.velocity)");
        }

        std::string vel_dist_type = toLower(velocity_cfg->get("type", "maxwell_boltzmann").asString());
        double temperature_K = velocity_cfg->get("temperature_K", 0.0).asDouble();
        Vec3 directed_velocity = velocity_cfg->isMember("directed_velocity")
                                     ? parseVec3((*velocity_cfg)["directed_velocity"], Vec3(0.0, 0.0, 0.0))
                                     : Vec3(0.0, 0.0, 0.0);
        Vec3 vel_fixed = velocity_cfg->isMember("value")
                             ? parseVec3((*velocity_cfg)["value"], directed_velocity)
                             : directed_velocity;

        for (size_t i = 0; i < initial_count; ++i) {
            core::IonState ion;

            if (pos_dist_type == "gaussian") {
                ion.pos.x = sampleGaussian(rng, pos_center.x, pos_sigma.x);
                ion.pos.y = sampleGaussian(rng, pos_center.y, pos_sigma.y);
                ion.pos.z = sampleGaussian(rng, pos_center.z, pos_sigma.z);
            } else if (pos_dist_type == "uniform_cylinder") {
                double sigma = std::max({pos_sigma.x, pos_sigma.y, 1e-3});
                ion.pos.x = sampleGaussian(rng, pos_center.x, sigma);
                ion.pos.y = sampleGaussian(rng, pos_center.y, sigma);
                ion.pos.z = pos_center.z;
            } else {
                ion.pos = pos_center;
            }

            if (vel_dist_type == "maxwell_boltzmann" || vel_dist_type == "thermal") {
                ion.vel.x = sampleMaxwellBoltzmannComponent(rng, temperature_K, mass_kg) + directed_velocity.x;
                ion.vel.y = sampleMaxwellBoltzmannComponent(rng, temperature_K, mass_kg) + directed_velocity.y;
                ion.vel.z = sampleMaxwellBoltzmannComponent(rng, temperature_K, mass_kg) + directed_velocity.z;
            } else if (vel_dist_type == "gaussian") {
                double sigma_v = std::sqrt(std::max(temperature_K, 0.0));
                ion.vel.x = sampleGaussian(rng, vel_fixed.x, sigma_v);
                ion.vel.y = sampleGaussian(rng, vel_fixed.y, sigma_v);
                ion.vel.z = sampleGaussian(rng, vel_fixed.z, sigma_v);
            } else if (vel_dist_type == "fixed") {
                ion.vel = vel_fixed;
            } else {
                throw std::runtime_error("Unsupported velocity_distribution type '" + vel_dist_type + "' for species '" + species_id + "'");
            }

            ion.mass_kg = mass_kg;
            ion.ion_charge_C = ion_charge_C;
            ion.species_id = species_id;
            ion.reduced_mobility_cm2_Vs = reduced_mobility;
            ion.CCS_m2 = ccs_m2;
            ion.active = true;
            ion.born = true;
            ion.birth_time_s = 0.0;
            ion.history_index = -1;
            ion.current_domain_index = 0;

            ions_.push_back(ion);
            ++total_generated;
        }
    }

    if (total_generated == 0) {
        throw std::runtime_error("Ion configuration produced zero ions; check ions.species initial_count values.");
    }

    std::cout << "[setupIons] Generated " << total_generated << " ions from v1.0 schema" << std::endl;

    ion_rngs_.resize(ions_.size());
    unsigned long base_seed = has_seed_override_ ? seed_override_value_ : config_.rng_seed;
    for (size_t i = 0; i < ions_.size(); ++i) {
        ion_rngs_[i].seed(base_seed + i);
    }
}

void SimulationEngine::setupIntegrator(const Json::Value& config, bool create_cpu) {
    if (!create_cpu) {
        // Do not instantiate CPU integrator (GPU-only run without validation)
        integrator_.reset();
        return;
    }

    std::string integrator_type = "RK4";
    if (config.isMember("integrator")) {
        const Json::Value& integrator_node = config["integrator"];
        if (integrator_node.isString()) {
            integrator_type = integrator_node.asString();
        } else if (integrator_node.isObject() && integrator_node.isMember("type") &&
                   integrator_node["type"].isString()) {
            integrator_type = integrator_node["type"].asString();
        }
    }
    
    auto factory_type = physics::IntegratorFactory::fromString(integrator_type);
    integrator_ = physics::IntegratorFactory::create(factory_type);
}

void SimulationEngine::setupInstrument(const Json::Value& config) {
    // Check if this is a multi-domain configuration
    const Json::Value& instrument_config = config["instrument"];
    
    if (instrument_config.isMember("domains") && instrument_config["domains"].isArray()) {
        // Multi-domain configuration
        is_multi_domain_ = true;
        instruments_ = instrument::InstrumentFactory::createMultiDomain(config);
        
        if (instruments_.empty()) {
            throw std::runtime_error("Failed to create multi-domain instruments from configuration");
        }
        
        // Use first domain's type as the primary type
        instrument_type_ = instruments_[0]->type();
        
        // Build GPU domains for all instrument domains
        instrument_domains_.clear();
        double z_offset = 0.0;  // Accumulate z positions for domain origins
        
        for (const auto& inst : instruments_) {
            auto gpu_domains = inst->buildGpuDomains();
            
            // Fix domain origins for multi-domain configuration
            // Each domain's origin should be at [0, 0, z_offset]
            for (auto& dom : gpu_domains) {
                dom.geom.origin_m = Vec3{0.0, 0.0, z_offset};
            }
            
            instrument_domains_.insert(instrument_domains_.end(), 
                                      gpu_domains.begin(), gpu_domains.end());
            
            // Advance z_offset for next domain
            if (!gpu_domains.empty()) {
                z_offset += gpu_domains[0].geom.length_m;
            }
        }
        
        std::cout << "[SimulationEngine] Multi-domain configuration with " 
                  << instruments_.size() << " domains" << std::endl;
        
        // For multi-domain GPU simulations, override GPU integrator domains now
        // (must happen here because setupForces() runs before gpu_mode_hint_ is set)
#ifdef USE_GPU_ACCEL
        std::cerr << "[SimulationEngine::setupInstrument] Checking GPU override: gpu_integrator_=" 
                  << (gpu_integrator_ ? "non-null" : "null") 
                  << ", instrument_domains_.size()=" << instrument_domains_.size() << std::endl;
        if (gpu_integrator_ && !instrument_domains_.empty()) {
            std::cerr << "[SimulationEngine::setupInstrument] Calling gpu_integrator_->overrideDomains with " 
                      << instrument_domains_.size() << " domains (multi-domain mode)" << std::endl;
            gpu_integrator_->overrideDomains(instrument_domains_);
        }
#endif
    } else {
        // Single-domain configuration (legacy)
        is_multi_domain_ = false;
        instrument_ = instrument::InstrumentFactory::create(config);
        
        if (!instrument_) {
            throw std::runtime_error("Failed to create instrument from configuration");
        }
        
        instrument_type_ = instrument_->type();
        instrument_domains_.clear();
        instrument_domains_ = instrument_->buildGpuDomains();
    }
}

void SimulationEngine::setupForces(const Json::Value& config) {
    auto composite = std::make_shared<physics::CompositeForces>();
    
    if (is_multi_domain_) {
        // Multi-domain: attach forces from all domains
        for (const auto& inst : instruments_) {
            inst->attachForces(*composite);
        }
    } else if (instrument_) {
        // Single-domain: attach forces from single instrument
        instrument_->attachForces(*composite);
    } else {
        // Fallback: create instrument from config
        auto fallback = instrument::InstrumentFactory::create(config);
        fallback->attachForces(*composite);
    }
    
    // Parse collision model from config (same logic as GPU integrator parsing)
    auto readCollisionModelString = [](const Json::Value& node)->std::string {
        if (node.isString()) return node.asString();
        if (node.isObject()) {
            if (node.isMember("model") && node["model"].isString()) return node["model"].asString();
            if (node.isMember("type") && node["type"].isString()) return node["type"].asString();
        }
        return {};
    };

    std::string collision = "EHSS";
    bool collision_set = false;
    if (config.isMember("physics") && config["physics"].isObject()) {
        const auto& phys = config["physics"];
        // Accept several possible key names used across configs/tests:
        // - "collision_model" (singular)
        // - "collisions" (object or string)
        // - legacy "collisions_model" (plural; some tests use this)
        if (phys.isMember("collision_model") && phys["collision_model"].isString()) {
            collision = phys["collision_model"].asString();
            collision_set = true;
        } else if (phys.isMember("collisions_model") && phys["collisions_model"].isString()) {
            collision = phys["collisions_model"].asString();
            collision_set = true;
        } else if (phys.isMember("collisions")) {
            const auto candidate = readCollisionModelString(phys["collisions"]);
            if (!candidate.empty()) { collision = candidate; collision_set = true; }
        }
    }
    if (!collision_set && config.isMember("collisions")) {
        const auto candidate = readCollisionModelString(config["collisions"]);
        if (!candidate.empty()) { collision = candidate; collision_set = true; }
    }
    // Some configurations/tests accidentally set the key name to "collisions_model"
    // (plural). Accept that legacy key as a top-level fallback as well.
    if (!collision_set && config.isMember("collisions_model") && config["collisions_model"].isString()) {
        collision = config["collisions_model"].asString();
        collision_set = true;
    }
    if (!collision_set && config.isMember("collision_model") && config["collision_model"].isString()) {
        collision = config["collision_model"].asString();
        collision_set = true;
    }
    std::transform(collision.begin(), collision.end(), collision.begin(), [](unsigned char c){ return static_cast<char>(std::tolower(c)); });

    // Build minimal GlobalParams and attach collision force adapter so CPU integrator gets damping
    core::GlobalParams gParams;
    if (collision == "hardsphere" || collision == "hard_sphere") gParams.collisionModel = core::CollisionModel::HardSphere;
    else if (collision == "langevin") gParams.collisionModel = core::CollisionModel::Langevin;
    else if (collision == "friction") gParams.collisionModel = core::CollisionModel::Friction;
    else if (collision == "ehss") gParams.collisionModel = core::CollisionModel::EHSS;
    else if (collision == "hsmc") gParams.collisionModel = core::CollisionModel::HSMC;
    else if (collision == "nocollisions" || collision == "none") gParams.collisionModel = core::CollisionModel::NoCollisions;
    else gParams.collisionModel = core::CollisionModel::EHSS;

    // Attach collision adapter using the instrument domain (needed for env parameters)
    if (instrument_) {
        // Build a core::InstrumentDomain POD from the instrument's JSON config
        // and pass it to the CollisionForces adapter so the adapter stores a
        // local copy of the domain parameters (avoids lifetime issues).
        try {
            // Build the full domain from the simulation config (which contains
            // the instrument, geometry and environment blocks expected by
            // load_single_domain). Using the top-level `config` ensures the
            // environment parameters (pressure/temperature) are present.
            ICARION_DEBUG_LOGF("config", "[DomainLoad] config has instrument=%d geometry=%d environment=%d", 
                                static_cast<int>(config.isMember("instrument")),
                                static_cast<int>(config.isMember("geometry")),
                                static_cast<int>(config.isMember("environment")));
            // load_single_domain expects a compact domain object with keys:
            // { "instrument": "IMS", "geometry": {...}, "environment": {...} }
            Json::Value dom_json(Json::objectValue);
            // instrument may be provided as an object with a "type" key or as a string
            if (config.isMember("instrument")) {
                if (config["instrument"].isString()) {
                    dom_json["instrument"] = config["instrument"].asString();
                } else if (config["instrument"].isObject() && config["instrument"].isMember("type")) {
                    dom_json["instrument"] = config["instrument"]["type"].asString();
                }
            }
            if (config.isMember("geometry")) dom_json["geometry"] = config["geometry"];
            if (config.isMember("environment")) dom_json["environment"] = config["environment"];
            core::InstrumentDomain dom = core::load_single_domain(dom_json);
            composite->addForceCalculator(std::make_shared<physics::CollisionForces>(gParams, dom));
            } catch (const std::exception& ex) {
            // If domain construction fails, log the reason and fall back to a
            // default-constructed domain so the CPU integrator still runs.
            ICARION_DEBUG_LOGF("config", "[WARN][SimulationEngine] load_single_domain failed: %s", ex.what());
            core::InstrumentDomain dom;
            composite->addForceCalculator(std::make_shared<physics::CollisionForces>(gParams, dom));
        }
    }

    // Debug: print which force calculators are registered (helps verify CollisionForces was added)
    try {
        auto names = composite->getCalculatorNames();
        std::string buf = "[Forces] registered calculators:";
        for (const auto& n : names) { buf += " "; buf += n; }
        ICARION_DEBUG_LOG("physics", buf);
    } catch (...) {}

    // Propagate integrator timestep to any force calculators that need it
    const Json::Value& sim_for_ts = getSimParams(config);
    double ts = config_.timestep; // Already parsed in parseConfiguration
    if (ts <= 0.0 && sim_for_ts.isMember("timestep")) {
        ts = sim_for_ts["timestep"].asDouble();
    } else if (ts <= 0.0 && sim_for_ts.isMember("timestep_ns")) {
        ts = sim_for_ts["timestep_ns"].asDouble() * 1e-9;
    }
    composite->setTimestep(ts);

    force_calculator_ = composite;
#ifdef USE_GPU_ACCEL
    std::cerr << "[SimulationEngine::setupForces] gpu_mode_hint_=" << (gpu_mode_hint_ ? "true" : "false") 
              << ", instrument_domains_.size()=" << instrument_domains_.size() << std::endl;
    if (gpu_mode_hint_) {
        if (instrument_domains_.empty() && instrument_) {
            instrument_domains_ = instrument_->buildGpuDomains();
        }
        if (instrument_domains_.empty()) {
            throw std::runtime_error("[SimulationEngine] Instrument returned no GPU domains.");
        }
        std::cerr << "[SimulationEngine] Calling gpu_integrator_->overrideDomains with " 
                  << instrument_domains_.size() << " domains" << std::endl;
        gpu_integrator_->overrideDomains(instrument_domains_);
    }
#endif
}

void SimulationEngine::loadCpuStochasticAssets(const Json::Value& config) {
    geometry_map_.clear();
    species_db_.clear();
    reaction_list_.clear();

    if (!cpu_stochastic_collisions_enabled_ && !cpu_reactions_enabled_) {
        return;
    }

    const double env_temp_K = resolveEnvironmentTemperature(config);
    const double neutral_mass_kg = resolveNeutralMassKg(config);

    if (cpu_reactions_enabled_) {
        const std::string reaction_db = resolveReactionDatabasePath(config);
        if (reaction_db.empty()) {
            std::cerr << "[SimulationEngine] CPU reactions requested but no reaction_database provided. Disabling reaction handling." << std::endl;
            cpu_reactions_enabled_ = false;
        } else {
            try {
                core::GlobalParams gParams;
                gParams.reaction_file = reaction_db;
                gParams.enable_reactions = cpu_reactions_enabled_;
                species_db_ = ::load_speciesDB(gParams, env_temp_K, neutral_mass_kg);
                reaction_list_ = ::load_reactions(gParams);
                std::cout << "[SimulationEngine] Loaded " << species_db_.size() << " species and "
                          << reaction_list_.size() << " reactions from '" << reaction_db << "'" << std::endl;
            } catch (const std::exception& ex) {
                std::cerr << "[SimulationEngine] Failed to load reaction database '" << reaction_db
                          << "': " << ex.what() << std::endl;
                species_db_.clear();
                reaction_list_.clear();
                cpu_reactions_enabled_ = false;
            }
        }
    }

    for (const auto& ion : ions_) {
        if (ion.species_id.empty()) {
            continue;
        }
        if (species_db_.find(ion.species_id) != species_db_.end()) {
            continue;
        }
        Species s{};
        s.name = ion.species_id;
        s.mass_kg = ion.mass_kg;
        s.mobility = ion.reduced_mobility_cm2_Vs * 1e-4;
        s.charge = ion.ion_charge_C;
        s.CCS = ion.CCS_m2;
        species_db_.emplace(s.name, s);
    }

    if (cpu_stochastic_collisions_enabled_ &&
        (cpu_collision_model_ == CollisionModel::EHSS || cpu_collision_model_ == CollisionModel::HSMC)) {
        const std::string geometry_file = resolveGeometryFilePath(config);
        if (geometry_file.empty()) {
            if (cpu_collision_model_ == CollisionModel::EHSS) {
                std::cerr << "[SimulationEngine] Warning: EHSS collision model requested but no geometry_file provided; falling back to hard-sphere collisions for CPU path." << std::endl;
            }
            return;
        }

        std::unordered_set<std::string> species_ids;
        species_ids.reserve(ions_.size());
        for (const auto& ion : ions_) {
            if (!ion.species_id.empty()) {
                species_ids.insert(ion.species_id);
            }
        }

        for (const auto& species_id : species_ids) {
            std::pair<std::vector<Vec3>, std::vector<double>> geometry_entry;
            if (!loadEhssGeometryForSpecies(geometry_file, species_id, geometry_entry)) {
                geometry_map_[species_id] = {std::vector<Vec3>(), std::vector<double>()};
                continue;
            }
            geometry_map_[species_id] = std::move(geometry_entry);
        }

        if (!geometry_map_.empty()) {
            std::cout << "[SimulationEngine] Loaded EHSS geometry for " << geometry_map_.size()
                      << " species from '" << geometry_file << "'" << std::endl;
        }
    }
}

void SimulationEngine::setupOutput(const Json::Value& config) {

    // ---- SSOT: require output.trajectory_file.filename ----
    if (!config.isMember("output") || !config["output"].isObject()) {
        throw std::runtime_error("Missing required block: output");
    }

    const auto& output_block = config["output"];

    if (!output_block.isMember("trajectory_file") || 
        !output_block["trajectory_file"].isObject()) {
        throw std::runtime_error("Missing required block: output.trajectory_file");
    }

    const auto& traj_file = output_block["trajectory_file"];

    if (!traj_file.isMember("filename") || !traj_file["filename"].isString()) {
        throw std::runtime_error("Missing required string: output.trajectory_file.filename");
    }

    config_.output_file = traj_file["filename"].asString();

    if (config_.output_file.empty()) {
        throw std::runtime_error("output.trajectory_file.filename cannot be empty");
    }

    // ---- Serialize config for writing to HDF5 ----
    Json::StreamWriterBuilder builder;
    builder["indentation"] = "";
    serialized_config_ = Json::writeString(builder, config);
    if (serialized_config_.empty()) {
        serialized_config_ = "{}";
    }

    std::cout << "Output configured: " << config_.output_file << std::endl;

    // ---- Initialize HDF5 writer ----
    hdf5_writer_ = std::make_unique<core::HDF5TrajectoryWriter>();
    hdf5_writer_->initialize(config_.output_file, ions_.size());

    output_initialized_ = true;
    output_finalized_ = false;
}


void SimulationEngine::performSimulationStep() {
    if (force_calculator_) {
        force_calculator_->setTime(current_time_);
    }
    auto force_function = createForceFunction();
    auto timed_force_function = createTimedForceFunction();
    bool step_performed = false;
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
    const auto countActiveIons = [](const std::vector<core::IonState>& ions) -> size_t {
        size_t count = 0;
        for (const auto& ion : ions) {
            if (ion.active) {
                ++count;
            }
        }
        return count;
    };
    auto logActivePhase = [this, &countActiveIons](const char* phase) {
        static size_t last_reported_active = std::numeric_limits<size_t>::max();
        const size_t current = countActiveIons(ions_);
        if (current != last_reported_active) {
            std::cout << "[CPU Active Debug] step=" << current_step_
                      << " phase=" << phase
                      << " active=" << current << std::endl;
            last_reported_active = current;
        }
    };
#endif
    
#ifdef USE_GPU_ACCEL
    if (gpu_mode_hint_ && gpu_integrator_ && gpu_integrator_->isAvailable()) {
#if 0
        std::cout << "[SIM-ENGINE-DEBUG] Entering GPU path with gpu_mode_hint_=" << gpu_mode_hint_ 
                  << " gpu_integrator_available=" << (gpu_integrator_ ? gpu_integrator_->isAvailable() : false) << std::endl;
#endif
        // Store CPU reference state for validation
        std::vector<core::IonState> cpu_reference;
        if (gpu_validation_enabled_) {
            cpu_reference = ions_;
        }
        try {
            // Perform GPU step
            if (gpu_validation_enabled_ || (current_step_==0 && current_time_==0.0)) {
                std::cout << "[SIM-DEBUG] GPU dt_s=" << config_.timestep << " current_time_s=" << current_time_ << std::endl;
            }
            gpu_integrator_->step(ions_, config_.timestep, current_time_, force_function);
            step_performed = true;
            if (gpu_validation_enabled_) {
                auto cpu_compare = cpu_reference;
                // Perform CPU step for comparison
                // CPU integrator should respect hard boundaries: backup previous and apply clipping
                backupIonStates();
                integrator_->step(cpu_compare, config_.timestep, current_time_, timed_force_function);
                applyBoundaryIntersections(cpu_reference, cpu_compare);
                // Side-by-side CPU vs GPU trace for ion 0 (validation mode)
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
                if (!cpu_compare.empty() && !ions_.empty()) {
                    const auto& cpu0 = cpu_compare.front();
                    const auto& gpu0 = ions_.front();
                    std::cout.setf(std::ios::scientific);
                    std::cout.precision(12);
                    std::cout << "[PARITY STEP] t_s=" << current_time_ <<
                                 " CPU0_pos=(" << cpu0.pos.x << "," << cpu0.pos.y << "," << cpu0.pos.z << ")" <<
                                 " GPU0_pos=(" << gpu0.pos.x << "," << gpu0.pos.y << "," << gpu0.pos.z << ")" << std::endl;
                }
#endif
                // Diagnostic: if requested, print per-ion velocity/acceleration comparison
                const char* dbg = std::getenv("ICARION_DEBUG_GPU_STEP");
                if (dbg && dbg[0] != '\0' && !cpu_compare.empty() && !cpu_reference.empty() && !ions_.empty()) {
                    try {
                        const auto& cpu_next = cpu_compare.front();
                        const auto& gpu_next = ions_.front();
                        // GPU-CPU validation diagnostics (disabled for production)
                        (void)cpu_next; (void)gpu_next; // Suppress unused warnings
                    } catch (...) {}
                }

                const auto [pos_rms, vel_rms] = computeRmsErrors(cpu_compare, ions_);
                if (pos_rms > gpu_validation_pos_tol_ || vel_rms > gpu_validation_vel_tol_) {
                    ICARION_DEBUG_LOGF("engine", "[SimulationEngine] GPU integrator deviated (pos RMS=%.6g, vel RMS=%.6g). Disabling GPU mode for the remainder of the run.", pos_rms, vel_rms);
                    ions_ = std::move(cpu_compare);
                    // Only disable the GPU hint if the user has not explicitly forced GPU execution.
                    if (execution_override_ != ExecutionMode::ForceGPU) {
                        gpu_mode_hint_ = false;
                    } else {
                        ICARION_DEBUG_LOG("engine", "[SimulationEngine] ForceGPU override active; retaining GPU mode despite validation deviation.");
                    }
                }
            }
        } catch (const std::exception& ex) {
            std::cerr << "GPU integrator error: " << ex.what()
                      << ". Falling back to CPU integrator." << std::endl;
            // Respect explicit ForceGPU override: do not flip the gpu_mode_hint_ off if the user forced GPU.
            if (execution_override_ != ExecutionMode::ForceGPU) {
                gpu_mode_hint_ = false;
            } else {
                ICARION_DEBUG_LOGF("engine", "[SimulationEngine] ForceGPU override active; GPU integrator threw but gpu_mode_hint_ retained. Exception: %s", ex.what());
            }
        }
    }
#endif

    if (!step_performed) {
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
        logActivePhase("pre_integrator");
#endif
        // CPU path: backup previous ion states and apply hard boundary clipping after step
        backupIonStates();
        try {
            integrator_->step(ions_, config_.timestep, current_time_, timed_force_function);
        } catch (...) {
            // if integrator throws, rethrow after attempting to apply boundaries to best-effort state
            applyBoundaryIntersections(previous_ions_, ions_);
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
            logActivePhase("post_integrator_exception");
#endif
            throw;
        }
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
        logActivePhase("post_integrator");
#endif
        
        // ========================================================================
        // CPU STOCHASTIC COLLISIONS AND REACTIONS (v1.0+ restored)
        // Applied ONCE per timestep AFTER RK4 integration completes
        // ========================================================================
        if (cpu_stochastic_collisions_enabled_ || cpu_reactions_enabled_) {
            for (size_t i = 0; i < ions_.size(); ++i) {
                if (!ions_[i].active) continue;
                
                // Apply stochastic collision (EHSS or HSMC)
                if (cpu_stochastic_collisions_enabled_) {
                    cpu_handle_stochastic_collision(
                        ions_[i],
                        ion_rngs_[i],
                        config_.timestep,
                        cpu_collision_model_,
                        geometry_map_
                    );
                }
                
                // Apply stochastic reactions
                if (cpu_reactions_enabled_) {
                    cpu_handle_stochastic_reaction(
                        ions_[i],
                        ion_rngs_[i],
                        config_.timestep,
                        cpu_reactions_enabled_,
                        species_db_,
                        reaction_list_
                    );
                }
            }
        }
        
        applyBoundaryIntersections(previous_ions_, ions_);
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
        logActivePhase("post_boundary");
#endif
    }

#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
    logActivePhase("pre_deactivate");
#endif
    deactivateOutsideIons();
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
    logActivePhase("post_deactivate");
#endif
    
    // Multi-domain: check for domain transitions after position update
    if (is_multi_domain_) {
        checkDomainTransitions();
    }
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
    logActivePhase("post_domain_transitions");
#endif

    // IMS-specific: enforce mobility-derived drift over the last timestep by
    // overriding the position increment and velocity for IMS ions. This uses
    // the backed-up `previous_ions_` states (populated before the integrator
    // step) to compute next_z = prev_z + v_target * dt, which avoids large
    // transient excursions while preserving per-ion mobility values.
    try {
        if (instrument_ && instrument_->type() == instrument::InstrumentType::IMS && !previous_ions_.empty()) {
            const Json::Value& inst = instrument_->instrumentConfig();
            double E_z = 0.0;
            if (inst.isMember("drift_field_V_m")) {
                E_z = inst["drift_field_V_m"].asDouble();
            } else if (inst.isMember("drift_voltage_V") && inst.isMember("length_m") && inst["length_m"].asDouble() > 0.0) {
                E_z = inst["drift_voltage_V"].asDouble() / inst["length_m"].asDouble();
            } else if (instrument_->metadata().isMember("drift_field_V_m")) {
                E_z = instrument_->metadata()["drift_field_V_m"].asDouble();
            }

            const double dt = config_.timestep;
            const size_t n = std::min(ions_.size(), previous_ions_.size());
            for (size_t i = 0; i < n; ++i) {
                auto& ion = ions_[i];
                if (!ion.active) continue;
                double K = ion.reduced_mobility_cm2_Vs * 1e-4; // m^2 / V s (reduced mobility)
                // Convert reduced mobility to actual mobility at current
                // gas density: mobi = K * LOSCHMIDT / particle_density
                double particle_density_m3 = 0.0;
                const Json::Value& inst_root = instrument_->instrumentConfig();
                if (inst_root.isMember("pressure_Pa") && inst_root.isMember("temperature_K")) {
                    double p = inst_root["pressure_Pa"].asDouble();
                    double T = inst_root["temperature_K"].asDouble();
                    if (T > 0.0) particle_density_m3 = p / (BOLTZMANN_CONSTANT * T);
                }
                // fallback to global environment block if present
                if (particle_density_m3 == 0.0) {
                    const Json::Value& root = instrument_->rootConfig();
                    if (root.isMember("environment") && root["environment"].isObject()) {
                        const auto& env = root["environment"];
                        if (env.isMember("pressure_Pa") && env.isMember("temperature_K")) {
                            double p = env["pressure_Pa"].asDouble();
                            double T = env["temperature_K"].asDouble();
                            if (T > 0.0) particle_density_m3 = p / (BOLTZMANN_CONSTANT * T);
                        }
                    }
                }
                if (particle_density_m3 <= 0.0) particle_density_m3 = LOSCHMIDT_CONSTANT; // fallback
                double actual_mobility = K * LOSCHMIDT_CONSTANT / particle_density_m3;
                double v_target_z = actual_mobility * E_z;
                // enforce position increment over last step to be consistent
                // with the mobility-driven velocity
                ion.pos.z = previous_ions_[i].pos.z + v_target_z * dt;
                // PARITY FIX: Comment out velocity override to match GPU approach
                // ion.vel.z = v_target_z;

                // If the instrument now considers the ion outside, attempt clip
                if (instrument_->isOutsideDomain(ion)) {
                    bool clipped = instrument_->clipToDomain(previous_ions_[i], ion);
                    if (clipped) {
                        ion.active = false;
                        ion.vel = Vec3{0.0, 0.0, 0.0};
                    } else {
                        ion.active = false;
                    }
                }
            }
        }
    } catch (...) {
        // ignore failures — fallback to integrator state
    }
}

namespace {

inline bool is_inside_domain_cpu_globalpos(const Vec3& global_pos,
                                           const DomainGPU& dom) {
    Vec3 diff;
    diff.x = global_pos.x - dom.geom.origin_m.x;
    diff.y = global_pos.y - dom.geom.origin_m.y;
    diff.z = global_pos.z - dom.geom.origin_m.z;

    Vec3 local;
    local.x = dom.geom.rot_row0.x * diff.x + dom.geom.rot_row0.y * diff.y + dom.geom.rot_row0.z * diff.z;
    local.y = dom.geom.rot_row1.x * diff.x + dom.geom.rot_row1.y * diff.y + dom.geom.rot_row1.z * diff.z;
    local.z = dom.geom.rot_row2.x * diff.x + dom.geom.rot_row2.y * diff.y + dom.geom.rot_row2.z * diff.z;

    double r = std::sqrt(local.x * local.x + local.y * local.y);
    const double tol = 0.0;

    if (dom.instrument != InstrumentGPU::Orbitrap) {
        // Use canonical boundary fields for CPU-side boundary checks.
        return (local.z >= tol && local.z < dom.boundary_length_m) && (r < dom.boundary_radius_m);
    }

    const double Rin = dom.geom.radius_in_m;
    const double Rout = dom.geom.radius_out_m;
    return (r >= Rin) && (r <= Rout);
}

}  // namespace

void SimulationEngine::deactivateOutsideIons() {
    // Multi-domain configuration: check all domains
    if (is_multi_domain_) {
        for (auto& ion : ions_) {
            if (!ion.active) continue;
            
            bool inside_any = false;
            for (const auto& inst : instruments_) {
                if (!inst->isOutsideDomain(ion)) {
                    inside_any = true;
                    break;
                }
            }
            
            if (!inside_any) {
                ion.active = false;
            }
        }
        return;
    }
    
    // Single-domain configuration: use instrument-level domain checks
    // Prefer instrument-level domain checks for CPU path (covers instruments
    // that don't provide GPU domain descriptors). Fall back to GPU-domain
    // geometry checks only when no instrument object exists.
    if (instrument_) {
        for (size_t idx = 0; idx < ions_.size(); ++idx) {
            auto& ion = ions_[idx];
            if (!ion.active) continue;
            if (instrument_->isOutsideDomain(ion)) {
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
                if (instrument_->type() == instrument::InstrumentType::Orbitrap) {
                    const double r = std::sqrt(ion.pos.x * ion.pos.x + ion.pos.y * ion.pos.y);
                    const double z = ion.pos.z;
                    double radius_in = 0.0;
                    double radius_out = 0.0;
                    double length = 0.0;
                    if (!instrument_domains_.empty()) {
                        const auto& dom = instrument_domains_.front();
                        radius_in = dom.geom.radius_in_m;
                        radius_out = dom.geom.radius_out_m;
                        length = dom.geom.length_m;
                    }
                    const bool axial_violation = (length > 0.0) ? (std::fabs(z) > 0.5 * length) : false;
                    const bool inner_violation = (radius_in > 0.0) ? (r < radius_in) : false;
                    const bool outer_violation = (radius_out > 0.0) ? (r > radius_out) : false;
                    std::cout << "[CPU Orbitrap Outside] ion=" << idx
                              << " r=" << r
                              << " z=" << z
                              << " axial_violation=" << axial_violation
                              << " inner_violation=" << inner_violation
                              << " outer_violation=" << outer_violation
                              << std::endl;
                }
#endif
                ion.active = false;
            }
        }
        return;
    }

    // No instrument wrapper available: use GPU domain geometry if present.
    if (instrument_domains_.empty()) {
        return;
    }

    for (auto& ion : ions_) {
        if (!ion.active) continue;
        bool inside_domain = false;
        for (const auto& dom : instrument_domains_) {
            if (is_inside_domain_cpu_globalpos(ion.pos, dom)) {
                inside_domain = true;
                break;
            }
        }
        if (!inside_domain) ion.active = false;
    }
}

void SimulationEngine::backupIonStates() {
    previous_ions_ = ions_;
}

void SimulationEngine::applyBoundaryIntersections(const std::vector<core::IonState>& prev,
                                                 std::vector<core::IonState>& next) {
    if (!instrument_) return;
    const size_t n = std::min(prev.size(), next.size());
    for (size_t i = 0; i < n; ++i) {
        if (!next[i].active) continue;
        // If instrument-specific outside test says ion is outside, attempt clipping
        if (instrument_->isOutsideDomain(next[i])) {
#ifdef ICARION_DEBUG_ORBITRAP_FIELDS
            if (instrument_->type() == instrument::InstrumentType::Orbitrap) {
                const double r = std::sqrt(next[i].pos.x * next[i].pos.x + next[i].pos.y * next[i].pos.y);
                const double z = next[i].pos.z;
                double radius_in = 0.0;
                double radius_out = 0.0;
                double length = 0.0;
                if (!instrument_domains_.empty()) {
                    const auto& dom = instrument_domains_.front();
                    radius_in = dom.geom.radius_in_m;
                    radius_out = dom.geom.radius_out_m;
                    length = dom.geom.length_m;
                }
                const bool axial_violation = (length > 0.0) ? (std::fabs(z) > 0.5 * length) : false;
                const bool inner_violation = (radius_in > 0.0) ? (r < radius_in) : false;
                const bool outer_violation = (radius_out > 0.0) ? (r > radius_out) : false;
                std::cout << "[CPU Orbitrap Clip] ion=" << i
                          << " prev_r=" << std::sqrt(prev[i].pos.x * prev[i].pos.x + prev[i].pos.y * prev[i].pos.y)
                          << " prev_z=" << prev[i].pos.z
                          << " new_r=" << r
                          << " new_z=" << z
                          << " axial_violation=" << axial_violation
                          << " inner_violation=" << inner_violation
                          << " outer_violation=" << outer_violation
                          << std::endl;
            }
#endif
            bool clipped = instrument_->clipToDomain(prev[i], next[i]);
            if (clipped) {
                // Instrument requested deactivation upon clipping
                next[i].active = false;
                next[i].vel = Vec3{0.0, 0.0, 0.0};
            } else {
                // Fallback: deactivate if still outside after clip attempt
                next[i].active = false;
            }
        }
    }
}

void SimulationEngine::updatePerformanceMetrics() {
    auto step_end_time = std::chrono::high_resolution_clock::now();
    auto step_duration = std::chrono::duration_cast<std::chrono::nanoseconds>(step_end_time - step_start_time_);
    
    accumulated_step_time_ += step_duration.count();
    performance_sample_count_++;
}

void SimulationEngine::updateEnergyConservation() {
    previous_energy_ = calculateTotalEnergy();
}

void SimulationEngine::writeOutput() {
    std::cout << "Writing output at time " << current_time_ << std::endl;
    if (!output_initialized_ || !hdf5_writer_) {
        return;
    }

#ifdef USE_GPU_ACCEL
    if (gpu_mode_hint_ && gpu_integrator_) {
        // Ensure host ions reflect latest GPU state before writing HDF5.
        try {
            gpu_integrator_->forceDownload();
        } catch (...) {
            // Don't let download failures break output generation; log and continue
            ICARION_DEBUG_LOG("gpuio", "[SimulationEngine] forceDownload() failed during writeOutput");
        }
    }
#endif

    hdf5_writer_->writeStep(current_time_, ions_);
}

void SimulationEngine::finalizeOutput() {
    if (!output_initialized_ || output_finalized_ || !hdf5_writer_) {
        return;
    }
    hdf5_writer_->finalize();
    output_finalized_ = true;
}

void SimulationEngine::updateOutputMetadata() {
    if (!output_initialized_ || !hdf5_writer_) {
        return;
    }

    std::ostringstream version_info;
    version_info << "ICARION Core";
#ifdef ICARION_VERSION
    version_info << " v" << ICARION_VERSION;
#endif
#ifdef GIT_HASH
    version_info << " (git " << GIT_HASH << ")";
#endif
#ifdef ICARION_BUILD_COMPILER
    version_info << " | compiler=" << ICARION_BUILD_COMPILER;
#endif
#ifdef ICARION_BUILD_FLAGS
    version_info << " | flags=" << ICARION_BUILD_FLAGS;
#endif
#ifdef ICARION_CUDA_VERSION
    version_info << " | cuda_runtime=" << ICARION_CUDA_VERSION;
#endif

    std::ostringstream perf;
    perf << "title=" << config_.title
         << ", timestep=" << config_.timestep << " s"
         << ", total_time=" << config_.total_time << " s"
         << ", steps=" << config_.total_steps
         << ", output_frequency=" << config_.output_frequency
         << ", rng_seed=" << config_.rng_seed;

    const std::string config_json = serialized_config_.empty() ? std::string("{}") : serialized_config_;
    hdf5_writer_->setMetadata(config_json, config_.rng_seed, version_info.str(), perf.str());
    if (instrument_) {
        hdf5_writer_->setInstrumentMetadata(instrument_->metadata());
    }
}

double SimulationEngine::calculateTotalEnergy() const {
    double kinetic = 0.0;
    
    for (const auto& ion : ions_) {
        double v_squared = ion.vel.x * ion.vel.x +
                           ion.vel.y * ion.vel.y +
                           ion.vel.z * ion.vel.z;
        kinetic += 0.5 * ion.mass_kg * v_squared;
    }
    
    return kinetic; // Simplified - only kinetic energy
}

std::pair<double, double> SimulationEngine::computeRmsErrors(
    const std::vector<core::IonState>& a,
    const std::vector<core::IonState>& b) const {
    if (a.size() != b.size() || a.empty()) {
        return {0.0, 0.0};
    }
    double pos_accum = 0.0;
    double vel_accum = 0.0;
    for (size_t i = 0; i < a.size(); ++i) {
        const auto& lhs = a[i];
        const auto& rhs = b[i];
        const double dx = lhs.pos.x - rhs.pos.x;
        const double dy = lhs.pos.y - rhs.pos.y;
        const double dz = lhs.pos.z - rhs.pos.z;
        pos_accum += dx * dx + dy * dy + dz * dz;

        const double dvx = lhs.vel.x - rhs.vel.x;
        const double dvy = lhs.vel.y - rhs.vel.y;
        const double dvz = lhs.vel.z - rhs.vel.z;
        vel_accum += dvx * dvx + dvy * dvy + dvz * dvz;
    }
    const double norm = static_cast<double>(a.size()) * 3.0;
    return {std::sqrt(pos_accum / norm), std::sqrt(vel_accum / norm)};
}

// Create a force function that captures the current force calculator
std::function<Vec3(const core::IonState&, size_t)> SimulationEngine::createForceFunction() {
    // Return zero force for inactive ions so integrators do not advance
    // deactivated particles (they may still be present in the ions_ vector).
    return [this](const core::IonState& ion, size_t ion_index) {
        if (!ion.active) return Vec3{0.0, 0.0, 0.0};
        return force_calculator_->calculateForce(ion, ion_index, ions_);
    };
}


std::function<Vec3(double, const core::IonState&, size_t)> SimulationEngine::createTimedForceFunction() {
    // Timed variant that respects ion.active to avoid advancing deactivated ions
    return [this](double time_s, const core::IonState& ion, size_t ion_index) {
        if (!ion.active) return Vec3{0.0, 0.0, 0.0};
        if (force_calculator_) {
            force_calculator_->setTime(time_s);
            return force_calculator_->calculateForce(ion, ion_index, ions_);
        }
        return Vec3{0.0, 0.0, 0.0};
    };
}

// Multi-domain support (v1.0+)
void SimulationEngine::checkDomainTransitions() {
    if (!is_multi_domain_ || instruments_.empty()) return;
    
    for (auto& ion : ions_) {
        if (!ion.active) continue;
        
        // Determine current domain from position
        size_t current_domain = getCurrentDomain(ion);
        
        // Check if ion has moved to next domain
        if (current_domain < instruments_.size() - 1) {
            // Check if ion crossed into next domain
            if (tryDomainTransition(ion, current_domain, current_domain + 1)) {
                // Transition successful - ion continues in next domain
                continue;
            }
        }
        
        // Check if ion is outside all domains
        bool inside_any = false;
        for (size_t d = 0; d < instruments_.size(); ++d) {
            if (!instruments_[d]->isOutsideDomain(ion)) {
                inside_any = true;
                break;
            }
        }
        
        if (!inside_any) {
            ion.active = false;
        }
    }
}

bool SimulationEngine::tryDomainTransition(core::IonState& ion, size_t from_domain, size_t to_domain) {
    if (to_domain >= instruments_.size()) return false;
    
    const auto& next_instrument = instruments_[to_domain];
    
    // Check if ion is inside the next domain
    if (next_instrument->isOutsideDomain(ion)) {
        return false; // Ion hasn't reached next domain yet
    }
    
    // Get aperture configuration from JSON if present
    const Json::Value& from_config = instruments_[from_domain]->instrumentConfig();
    double aperture_radius = -1.0; // -1 means no aperture (100% transmission)
    
    if (from_config.isMember("aperture") && from_config["aperture"].isObject()) {
        const Json::Value& aperture = from_config["aperture"];
        if (aperture.isMember("radius_m")) {
            aperture_radius = aperture["radius_m"].asDouble();
        }
    }
    
    // Check aperture transmission
    if (aperture_radius > 0.0) {
        double r = std::sqrt(ion.pos.x * ion.pos.x + ion.pos.y * ion.pos.y);
        
        if (r > aperture_radius) {
            // Ion blocked by aperture
            ion.active = false;
            return false;
        }
    }
    
    // Transition successful
    return true;
}

size_t SimulationEngine::getCurrentDomain(const core::IonState& ion) const {
    if (!is_multi_domain_ || instruments_.empty()) return 0;
    
    // Find which domain contains this ion based on z-position
    for (size_t d = 0; d < instruments_.size(); ++d) {
        if (!instruments_[d]->isOutsideDomain(ion)) {
            return d;
        }
    }
    
    // Ion is outside all domains - return last domain
    return instruments_.size() - 1;
}

} // namespace simulation
} // namespace ICARION

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        paramUtils.h
 *   @brief       Parameter loading utilities for ICARION.
 *
 * @details
 * Provides routines to read, validate, and initialize simulation parameters
 * from a JSON configuration file. Uses structured JSON format with three
 * main sections: simulation, physics, and output. Constructs derived 
 * quantities such as the simulation time array (`t_eval`) and maps string 
 * identifiers (collision models, integrators) to enumerated types.
 *
 *   @date        2025-11-09
 *   @version     1.0.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#pragma once

#include <string>
#include <vector>
#include <unordered_map>

#include "utils/constants.h"
#include "physics/collisions/collisionHelpers.h"
#include "core/types/IonState.h"
#include "core/utils/mathUtils.h"
#include "json/json.h"
#include "core/io/fieldArrayLoader.h"
#include "core/physics/spacecharge/spaceChargeSolver.h"
// RK45Settings is a small POD type used by GlobalParams; include its
// dedicated header to avoid circular includes with integrator.h.
#include "integrator/rk45_settings.h"
// Numerical safety configuration
#include "core/utils/safety/numericalSafetyGuards.h"
// Instrument type definitions
#include "instrument/InstrumentTypes.h"
// SSOT: CollisionModel defined in config::PhysicsEnums.h (Phase 2E)
#include "core/config/types/PhysicsEnums.h"

// Forward declarations
class FieldServer; // Phase 6.3 live field computation

namespace ICARION {
namespace core {

// -----------------------------
// Allowed instruments, collision models and solvers
// -----------------------------
// Use canonical InstrumentType from instrument namespace
using Instrument = ICARION::instrument::InstrumentType;

// Legacy core::CollisionModel deleted - use config namespace as single source of truth
using CollisionModel = ICARION::config::CollisionModel;

enum class SolverType { RK45, RK4, Boris};

/**
 * @brief Stores simulation-wide global parameters.
 * 
 * Contains parameters that are constant across all instrument domains.
 * Includes collision model, general simulation settings (time step, number of ions, 
 * parallelization, output files), and runtime flags.
 */
struct GlobalParams {
    //Simulation
    CollisionModel collisionModel;
    int num_ions;

    // Time parameters
    int write_interval;
    std::vector<double> t_eval;
    double dt_s;
    int sim_time_steps;

    // Runtime flags
    bool parallelization;
    bool enable_gpu = false;    // request to use GPU acceleration when available
    bool enable_reactions;
    bool enable_space_charge;
    bool print_results;

    // Reproducibility
    unsigned int rng_seed = 42;  // Random number generator seed for reproducibility

    // Output / files
    std::string output_file;
    std::string species_database_file;  // NEW: species database (new input system)
    std::string reaction_file;          // OLD: kept for backward compatibility
    std::string ion_cloud_file;
    std::string geometry_file;
    std::string input_file;

    // Thermalization OU kick for deterministic models (e.g., Friction, Langevin, HardSphere)
    bool enable_ou_thermalization = false;
    // Force OU thermalization for stochastic collision models when true
    // (some callers may set this to ensure consistent behaviour across models)
    bool force_ou_for_stochastic_models = false;

    // Continue mode
    std::string continue_from;           // Path to HDF5 file to continue from
    double continue_time_s = 0.0;        // Additional simulation time [s]
    bool auto_continue_if_active = false; // Auto-detect incomplete runs

    // Space charge solver
    SpaceChargeSolver* spaceChargeSolver = nullptr; // runtime pointer

    // Integrator RK45 settings (can be populated from adapter/JSON)
    RK45Settings rk45_settings;
    
    // Numerical safety settings (for NaN/Inf detection and bounds checking)
    ICARION::safety::NumericalSafetyConfig numerical_safety;

};

/**
 * @brief Contains DC voltage parameters.
 * 
 */
struct DCVoltages {
    double axial_V = 0.0;
    double EN_Td = 0.0;
    double EN_Vm2 = 0.0;
    double radial_V = 0.0;
    double quad_V = 0.0;
    bool enable_radial_voltage_sweep = false;
    double radial_slope_V_s = 0.0;
    double radial_start_time_s = 0.0;
    double radial_rise_time_s = 0.0;
};

/**
 * @brief Contains RF voltage parameters.
 * 
 */
struct RFVoltages {
    double voltage_V = 0.0;
    double frequency_Hz = 0.0;
    double angular_frequency_rad_s = 0.0;
    double phase_rad = 0.0; 
};

/**
 * @brief Contains AC voltage parameters (for excitation in LQIT).
 * 
 */
struct ACVoltages {
    double voltage_V = 0.0;
    double frequency_Hz = 0.0;
    double angular_frequency_rad_s = 0.0;
    bool enable_voltage_sweep = false;
    double amplitude_slope_V_s = 0.0;
    double start_time_s = 0.0;
    double rise_time_s = 0.0;
    // Frequency sweep (linear in Hz) support
    bool enable_frequency_sweep = false;
    double ac_start_freq_Hz = 0.0;
    double ac_sweep_slope_Hz_per_s = 0.0;
    // ICARION / LQIT locking parameters
    bool lqit_lock_enable = false;
    double lqit_lock_phase_rad = 0.0;
    double lqit_lock_bandwidth_Hz = 0.0;
    // Optional user-defined time->voltage table for arbitrary AC voltage waveforms.
    // Each entry is a pair {time_s, voltage_V}. The loader will populate this
    // vector if the JSON contains an array of [time, voltage] samples under
    // AC.voltage_time_table. The runtime will linearly interpolate between
    // samples.
    std::vector<std::pair<double,double>> voltage_time_table;
};

/**
 * @brief Contains magnetic field parameters.
 * 
 */
struct MagneticField {
    bool enabled = false;
    Vec3 field_strength_T = Vec3{0.0, 0.0, 0.0}; //constant B
    Vec3 field_gradient_T_m = Vec3{0.0, 0.0, 0.0}; //optional linear gradient
};

/**
 * @brief Contains geometry-related parameters.
 * 
 */
struct Geometry {
    double length_m;
    double radius_m;
    double radius_in_m;
    double radius_out_m;
    double radius_char_m;
    double acc_length_m;
    double end_aperture_m;

    // Domain bounds & transforms
    Vec3 origin_m;
    Vec3 min_bound;
    Vec3 max_bound;
};

/**
 * @brief Contains environment parameters.
 * 
 */
struct Environment {
    double pressure_Pa;
    double temperature_K;
    double particle_density_m_3;
    double mean_thermal_velocity_m_s;
    std::string neutral_species_id;
    double neutral_mass_kg;
    double neutral_polarizability_m3;
    double neutral_radius_m = 0.0; // Hard-sphere radius (m) for EHSS collisions
    Vec3 gas_velocity_m_s;
};

/**
 * @brief Stores instrument domain-specific parameters.
 * 
 * Contains parameters that are different across different instrument domains.
 * Includes geometry, environment and field settings (DC, RF and (for LQIT) AC excitation), 
 * and runtime flags.
 */
struct InstrumentDomain {
    int index = -1;
    Instrument instrument;

    // Solver selection
    SolverType solver_type;

    // Geometry
    Geometry geom;

    // Environment
    Environment env;

    // Potential arrays (precomputed fields)
    // Legacy single-file support
    std::string FA_file;
    FieldArray fieldArray;
    bool fieldArrayLoaded = false;

    // New: Superposition of multiple precomputed unit-voltage fields
    enum class FAScaleKind { Constant, DC_Axial, DC_Quad, DC_Radial, RF };
    struct FieldArrayTerm {
        std::string file;            // path to HDF5 file with unit-voltage field
        FieldArray field;            // loaded field data
        FAScaleKind kind = FAScaleKind::Constant; // scaling type
        double constant = 1.0;       // constant scaling factor (in Volts)
        double phase_rad = 0.0;      // additional phase offset for RF terms
        double frequency_Hz = 0.0;   // optional override for RF frequency (0 -> use dom.RF.frequency_Hz)
        bool loaded = false;         // load status
    };
    std::vector<FieldArrayTerm> FA_terms; // if non-empty, overrides FA_file

    // Voltages
    DCVoltages DC;
    RFVoltages RF;
    ACVoltages AC;

    // Magnetic field
    MagneticField B;

    // Live field computation 
    bool use_grid_field = false;  ///< If true, use FieldServer instead of analytic fields
    class FieldServer* fieldServer = nullptr; ///< Pointer to live field server (runtime)

    // Computed rotation matrices
    Mat3 rotation_global_to_local = Mat3::identity();
    Mat3 rotation_local_to_global = Mat3::identity();
};

// -----------------------------
// Required parameters per instrument
// -----------------------------
extern std::unordered_map<Instrument, std::vector<std::string>> required_params;

// -----------------------------
// Load and update functions
// -----------------------------

/**
 * @brief Load global simulation parameters from JSON configuration file
 * 
 * @param filename Path to JSON configuration file
 * 
 * @return GlobalParams structure with all simulation settings
 * 
 * Expects structured JSON format with three main sections:
 * - simulation: Time parameters, integrator selection, GPU/parallelization flags
 * - physics: Collision model, thermalization settings
 * - output: Output folder, trajectory file, progress printing
 * 
 * Additionally parses:
 * - domains: Array of instrument domain configurations
 * - ions: Ion cloud file and species information
 * - species_database: Path to species definitions
 * - reaction_database: Path to reaction definitions
 * 
 * Performs validation checks and throws std::runtime_error on:
 * - Missing required sections or parameters
 * - Invalid parameter types or values
 * - Unknown collision models or integrators
 */
GlobalParams load_global_params(const std::string& filename);

/**
 * @brief Load single instrument domain from JSON object
 * 
 * @param j JSON object containing domain configuration
 * 
 * @return InstrumentDomain structure with geometry, fields, and environment
 * 
 * Expects structured domain configuration with:
 * - instrument: Type identifier (LQIT, SIFDT-MS, IMS, Quadrupole, TOF, Orbitrap, FT-ICR)
 * - geometry: Dimensions (length_m, radius_m, optional origin_m)
 * - environment: Gas conditions (pressure_Pa, temperature_K, gas_species)
 * - fields: Electric/magnetic field configuration
 *   - DC: Axial voltage or field strength (field_strength_Td, axial_voltage_V)
 *   - RF: Voltage and frequency (voltage_V, frequency_Hz)
 *   - AC: Excitation parameters (voltage_V, frequency_Hz, sweep settings)
 *   - B: Magnetic field (enabled, field_strength_T, field_gradient_T_m)
 * - integrator: Solver selection (RK4, RK45, Boris) [optional]
 * - field_array: Precomputed field arrays from BEM/FEM solvers [optional]
 * 
 * Throws std::runtime_error on:
 * - Missing required fields (instrument, geometry, environment)
 * - Unknown instrument types or gas species
 * - Invalid parameter types or non-physical values
 * 
 * @note After loading, the domain is validated with `sanity_check_domain`.
 */
InstrumentDomain load_single_domain(const Json::Value& j);

// -----------------------------
// Guard functions
// -----------------------------

/**
 * @brief Validate global parameters and check for common configuration errors
 * 
 * @param gParams Global simulation parameters
 * 
 * Checks:
 * - Timespan and timestep are positive
 * - Output file paths are valid
 * - Species/reaction files exist (if specified)
 * - RNG seed is reasonable
 * 
 * Throws std::runtime_error on validation failure.
 */
void run_guard_check_global(const GlobalParams& gParams);

/**
 * @brief Validate instrument domain configuration
 * 
 * @param dom Instrument domain to validate
 * @param cm Collision model being used
 * 
 * Checks:
 * - Geometry dimensions are positive
 * - Voltages and frequencies are reasonable
 * - Environment parameters (pressure, temperature) are physical
 * - Collision model is compatible with domain settings
 * - Field arrays are loaded if required
 * 
 * Throws std::runtime_error on validation failure.
 */
void run_guard_check_domain(const InstrumentDomain& dom, CollisionModel cm);

/**
 * @brief Additional domain sanity checks (non-fatal warnings)
 * 
 * @param dom Instrument domain to check
 * 
 * Prints warnings for:
 * - Unusual parameter combinations
 * - Potentially inefficient settings
 * - Missing optional features
 */
void sanity_check_domain(const InstrumentDomain& dom);

// -----------------------------
// Utilities to find current instrument domain
// -----------------------------

/**
 * @brief Check if position is inside domain boundaries
 * 
 * @param dom Instrument domain
 * @param pos Position vector [m]
 * 
 * @return true if position is within domain bounds
 */
bool isInsideDomain(const InstrumentDomain& dom, const Vec3& pos);

/**
 * @brief Find which domain contains given position
 * 
 * @param pos Position vector [m]
 * @param domains Vector of all instrument domains
 * 
 * @return Domain index, or -1 if position is outside all domains
 */
int find_domain_index(const Vec3& pos, const std::vector<InstrumentDomain>& domains);

/**
 * @brief Find instrument domain containing given position
 * 
 * @param domains Vector of all instrument domains
 * @param pos Position vector [m]
 * 
 * @return Pointer to domain, or nullptr if position is outside all domains
 */
const InstrumentDomain* findInstrumentDomain(const std::vector<InstrumentDomain>& domains, const Vec3& pos);

/**
 * @brief Calculate Orbitrap radial position for given axial position
 * 
 * @param z Axial position [m]
 * @param R Outer electrode radius [m]
 * @param Rm Central electrode radius [m]
 * 
 * @return Radial position [m] on Orbitrap electrode surface
 * 
 * Computes r(z) for hyperlogarithmic Orbitrap geometry.
 */
double orbitrap_r_for_z(double z, double R, double Rm);

}  // namespace core
}  // namespace ICARION

// Bring core types into global namespace for backward compatibility
using ICARION::core::Instrument;
using ICARION::core::CollisionModel;
using ICARION::core::SolverType;
using ICARION::core::GlobalParams;
using ICARION::core::DCVoltages;
using ICARION::core::RFVoltages;
using ICARION::core::ACVoltages;
using ICARION::core::MagneticField;
using ICARION::core::Environment;
using ICARION::core::Geometry;
using ICARION::core::InstrumentDomain;
using ICARION::core::required_params;
using ICARION::core::load_global_params;
using ICARION::core::load_single_domain;
using ICARION::core::isInsideDomain;
using ICARION::core::find_domain_index;
using ICARION::core::findInstrumentDomain;
using ICARION::core::orbitrap_r_for_z;
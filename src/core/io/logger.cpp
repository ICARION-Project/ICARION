/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        logger.cpp
 *   @brief       Logs simulation parameters.
 *
 *   @details
 *   Provides functions to log simulation parameters separated by 
 *   global and instrument domain settings. Supports runtime appends.
 *
 *
 *   @date        2025-10-16
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#include "core/io/logger.h"
#include "core/io/speciesLoader.h"
#include <chrono>
#include <iomanip>
#include <iostream>
#include <unistd.h>
#include <limits.h>

namespace ICARION {
namespace io {

#ifdef ICARION_VERSION
#  define PROJECT_VERSION_STR ICARION_VERSION
#elif defined(ICARION_VERSION)
#  define PROJECT_VERSION_STR ICARION_VERSION
#else
#  define PROJECT_VERSION_STR "unknown"
#endif

#ifdef GIT_COMMIT_HASH
#  define PROJECT_GIT_HASH_STR GIT_COMMIT_HASH
#else
#  define PROJECT_GIT_HASH_STR "unknown"
#endif

// --- Helpers ---
std::string RunLogger::timestamp() const {
    auto now = std::chrono::system_clock::now();
    std::time_t now_c = std::chrono::system_clock::to_time_t(now);
    std::tm tm = *std::localtime(&now_c);
    std::ostringstream ss;
    ss << std::put_time(&tm, "%Y-%m-%d %H:%M:%S");
    return ss.str();
}

// --- Constructor ---
RunLogger::RunLogger(const std::string& output_file_base) {
    filepath_ = output_file_base + ".log";
    file_.open(filepath_, std::ios::out | std::ios::trunc);
    if (!file_.is_open()) throw std::runtime_error("Failed to create log file: " + filepath_);
}

// --- String helpers for enums ---
static std::string solver_to_string(SolverType s) {
    switch (s) {
        case SolverType::RK4:    return "RK4";
        case SolverType::RK45:   return "RK45";
        case SolverType::Boris:  return "Boris";
        default: return "Unknown";
    }
}

static std::string collision_model_to_string(CollisionModel c) {
    switch (c) {
        case CollisionModel::NoCollisions: return "NoCollisions";
        case CollisionModel::Friction:     return "Friction";
        case CollisionModel::HardSphere:   return "HardSphere";
        case CollisionModel::Langevin:     return "Langevin";
        case CollisionModel::EHSS:         return "EHSS";
        case CollisionModel::HSMC:         return "HSMC";
        case CollisionModel::UnknownCollisionModel: return "Unknown collision model";
        default: return "Unknown";
    }
}

static std::string instrument_to_string(Instrument i) {
    switch (i) {
        case Instrument::FTICR:   return "FT-ICR";
        case Instrument::Orbitrap: return "Orbitrap";
        case Instrument::LQIT:     return "LQIT";
        case Instrument::IMS:    return "IMS";
        case Instrument::QuadrupoleRF: return "Quadrupole";
        case Instrument::TOF:     return "TOF";
        case Instrument::NoFixedInstrument: return "No fixed instrument";
        case Instrument::UnknownInstrument: return "Unknown instrument";
        default: return "Unknown";
    }
}

// --- Write banner ---
void RunLogger::writeHeader() {
    file_ << "============================================================\n";
    file_ << "             ICARION Ion Collision And Reaction IntegratiON\n";
    file_ << "============================================================\n";
    file_ << "Simulation started: " << timestamp() << "\n\n";
}

// --- Write global parameters ---
void RunLogger::writeGlobalParams(const GlobalParams& g, const std::vector<IonState>& ions,
                                  const ICARION::io::SpeciesDatabase& speciesDB,
                                  const std::vector<ReactionEntry>& reaction_list) {
    char buf[HOST_NAME_MAX];

    file_ << "User               : " << getenv("USER") << "\n";
    if (gethostname(buf, sizeof(buf)) == 0)
        file_ << "Hostname           : " << buf << "\n";
    else
        file_ << "Hostname           : (unknown)\n";
    file_ << "ICARION Version    : " << PROJECT_VERSION_STR << "\n";
    file_ << "Git Commit Hash    : " << PROJECT_GIT_HASH_STR << "\n";
    file_ << "JSON Config File   : " << g.input_file << "\n";
    file_ << "------------------------------------------------------------\n";

    file_ << "[Global simulation Parameters]\n";
    file_ << "------------------------------------------------------------\n";
    file_ << "Time step [s]           : " << g.dt_s << "\n";
    file_ << "Integration end [s]     : " << g.dt_s*g.sim_time_steps << "\n";
    file_ << "Collision model         : " << collision_model_to_string(g.collisionModel) << "\n";
    file_ << "Write interval          : " << g.write_interval << "\n";
    file_ << "Output base file        : " << g.output_file << "\n";
    file_ << "Print results           : " << (g.print_results ? "true" : "false") << "\n";
    file_ << "RNG seed (reproducibility): " << g.rng_seed << "\n";
    file_ << "Number of ions        : " << ions.size() << "\n";
    file_ << "Species loaded        : " << speciesDB.size() << "\n";
    file_ << "Reactions loaded      : " << reaction_list.size() << "\n";
    if (g.enable_reactions) {
        file_ << "Reactions enabled      : true\n";
    } else {
        file_ << "Reactions enabled      : false\n";
    }
    file_ << "------------------------------------------------------------\n\n";
    file_.flush();
}

// --- Write instrument domains ---
void RunLogger::writeInstrumentDomains(const std::vector<InstrumentDomain>& domains) {
    file_ << "[Instrument Domains]\n";
    file_ << "------------------------------------------------------------\n";
    for (size_t i = 0; i < domains.size(); ++i) {
        const auto& d = domains[i];
        file_ << ">>> Domain #" << i << " (" << instrument_to_string(d.instrument) << ")\n";
        file_ << "Solver                  : " << solver_to_string(d.solver_type) << "\n";
        file_ << "  Geometry:\n";
        file_ << "    Origin [m]      : (" << d.geom.origin_m.x << ", "
                                             << d.geom.origin_m.y << ", "
                                             << d.geom.origin_m.z << ")\n";
        file_ << "    Length [m]      : " << d.geom.length_m << "\n";
        file_ << "    Radius [m]      : " << d.geom.radius_m << "\n";
        if (d.geom.end_aperture_m > 0.0) {
            file_ << "  End aperture [m] : " << d.geom.end_aperture_m << "\n";
        }
        if (i < domains.size()-1) {
            file_ << "  -> Connected to next domain through aperture #" << i << "\n";
        }
        if (d.fieldArrayLoaded) {
            file_ << "  Field Array loaded from file: " << d.FA_file << "\n";
        }
        if (d.instrument == Instrument::Orbitrap) {
            file_ << "    Inner radius [m]: " << d.geom.radius_in_m << "\n";
            file_ << "    Outer radius [m]: " << d.geom.radius_out_m << "\n";
            file_ << "    Char radius [m] : " << d.geom.radius_char_m << "\n";
        }
        file_ << "    Voltages and Fields:\n";
        file_ << "    DC axial [V]    : " << d.DC.axial_V << "\n";
        file_ << "    DC radial [V]   : " << d.DC.radial_V << "\n";
        if (d.DC.enable_radial_voltage_sweep) {
            file_ << "    DC radial sweep enabled: slope [V/s]: " << d.DC.radial_slope_V_s
                  << ", start time [s]: " << d.DC.radial_start_time_s
                  << ", rise time [s]: " << d.DC.radial_rise_time_s << "\n";
        }
        file_ << "    DC EN_Td [Td]   : " << d.DC.EN_Td << "\n";
        file_ << "    RF voltage [V]  : " << d.RF.voltage_V << "\n";
        file_ << "    RF frequency[Hz]: " << d.RF.frequency_Hz << "\n";
        file_ << "    AC excitation[V]: " << d.AC.voltage_V << "\n";
        file_ << "    AC frequency[Hz]: " << d.AC.frequency_Hz << "\n";
        if (d.AC.enable_voltage_sweep) {
            file_ << "    AC voltage sweep enabled: slope [V/s]: " << d.AC.amplitude_slope_V_s
                  << ", start time [s]: " << d.AC.start_time_s
                  << ", rise time [s]: " << d.AC.rise_time_s << "\n";
        }
        if (d.B.enabled) {
            file_ << "    Magnetic [T]    : (" << d.B.field_strength_T.x << ", "
                                                 << d.B.field_strength_T.y << ", "
                                                 << d.B.field_strength_T.z << ")\n";
        }
        file_ << "  Environment:\n";
        file_ << "    Pressure [Pa]   : " << d.env.pressure_Pa << "\n";
        file_ << "    Temperature [K] : " << d.env.temperature_K << "\n";
        file_ << "    Particle density [m^-3]: " << d.env.particle_density_m_3 << "\n";
        file_ << "    Neutral species : " << d.env.neutral_species_id << "\n";
        file_ << "    Gas velocity [m/s]: (" << d.env.gas_velocity_m_s.x << ", "
                                                 << d.env.gas_velocity_m_s.y << ", "
                                                 << d.env.gas_velocity_m_s.z << ")\n";
        file_ << "------------------------------------------------------------\n";
    }
    file_ << "\nLive Log:\n";
    file_ << "------------------------------------------------------------\n";
    file_.flush();
}

// --- Runtime message ---
void RunLogger::log(const std::string& msg) {
    file_ << "[" << timestamp() << "] " << msg << "\n";
    file_.flush();
}

// --- Final block ---
void RunLogger::finalize(const SimulationResult& result, const std::string& output_file) {
    // --- Compute ion stats ---
    size_t total   = result.ions.size();
    size_t active  = std::count_if(result.ions.begin(), result.ions.end(),
                                   [](const IonState& i){ return i.active; });
    size_t lost    = total - active;
    double frac_lost = total > 0 ? 100.0 * static_cast<double>(lost) / total : 0.0;

    // --- Compute output file size ---
    double file_size_MB = 0.0;
    std::ifstream f(output_file, std::ifstream::ate | std::ifstream::binary);
    if (f.is_open()) {
        file_size_MB = static_cast<double>(f.tellg()) / (1024.0 * 1024.0);
    }

    // --- Write block ---
    file_ << "------------------------------------------------------------\n";
    file_ << "Summary:\n";
    file_ << "  Active ions remaining    : " << active << "\n";
    file_ << "  Lost ions (boundary)     : " << lost
          << " (" << std::fixed << std::setprecision(1) << frac_lost << " %)\n";
    file_ << "  Output file size         : "
          << std::setprecision(1) << file_size_MB << " MB\n";
    file_ << "------------------------------------------------------------\n\n";

    file_ << "------------------------------------------------------------\n";
    file_ << "Simulation finished: " << timestamp() << "\n";
    file_ << "------------------------------------------------------------\n";
    file_ << "Project repository:\n";
    file_ << "  https://github.com/chsch95/Ion_Motion_Modelling\n";
    file_ << "Please cite ICARION once the reference paper is available.\n";
    file_ << "------------------------------------------------------------\n";
    file_.close();
    std::cout << "Log written to " << filepath_ << std::endl;
}

// Simple runtime debug logger (defines the free function declared in logger.h)
void debug_log(const std::string& msg) {
    std::cerr << msg << std::endl;
}

}  // namespace io
}  // namespace ICARION

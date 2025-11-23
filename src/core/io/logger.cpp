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
#include "core/log/Logger.h"  // For structured logging in debug_log()
#include <chrono>
#include <iomanip>
#include <iostream>
#include <fstream>
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

// The new config system uses ICARION::config enums with EnumMapper
// Note: solver_to_string, collision_model_to_string, instrument_to_string removed
// (legacy - were only used by writeInstrumentDomains which is now deleted)

// --- Write banner ---
void RunLogger::writeHeader() {
    file_ << "============================================================\n";
    file_ << "             ICARION Ion Collision And Reaction IntegratiON\n";
    file_ << "============================================================\n";
    file_ << "Simulation started: " << timestamp() << "\n\n";
}

// --- Write global parameters (DEPRECATED - disabled) ---
#if 0
void RunLogger::writeGlobalParams(const GlobalParams& g, const std::vector<IonState>& ions,
                                  const ICARION::io::SpeciesDatabase& speciesDB,
                                  const std::vector<int>& reaction_list) {  // Changed to int to avoid ReactionEntry
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
#endif

// Note: writeInstrumentDomains() removed (legacy GlobalParams/InstrumentDomain)

// --- Runtime message ---
void RunLogger::log(const std::string& msg) {
    file_ << "[" << timestamp() << "] " << msg << "\n";
    file_.flush();
}

// --- Final block ---
void RunLogger::finalize(const std::vector<IonState>& ions, const std::string& output_file) {
    // --- Compute ion stats ---
    size_t total   = ions.size();
    size_t active  = std::count_if(ions.begin(), ions.end(),
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
    file_ << "  https://github.com/ICARION-Project/ICARION\n";
    file_ << "Please cite ICARION once the reference paper is available.\n";
    file_ << "------------------------------------------------------------\n";
    file_.close();
    std::cout << "Log written to " << filepath_ << std::endl;
}

// Runtime debug logger - respects log level configuration
// Only outputs when log level is DEBUG (e.g., --verbose or --log-level DEBUG)
void debug_log(const std::string& msg) {
    // Use structured logger if initialized, otherwise fall back to stderr
    auto logger = ICARION::log::Logger::main();
    if (logger) {
        logger->debug("{}", msg);
    } else {
        // Fallback for early initialization (before Logger::init())
        std::cerr << "[DEBUG] " << msg << std::endl;
    }
}

}  // namespace io
}  // namespace ICARION

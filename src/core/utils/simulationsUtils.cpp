/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        simulationsUtils.cpp
 *   @brief       Utility functions for TRACE simulations.
 *
 * Provides routines to initialize ion clouds, optionally print simulation results,
 * and manage solver selection based on instrument type.
 *
 * @details
 * - Initializes ions from a JSON ion cloud file with species lookup.
 * - Populates ion physical properties (mass, charge, CCS, mobility) and
 *   environmental parameters based on the domain they are in.
 * - Optional printing of results (positions, velocities, arrival times) for
 *   a subset of ions.
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 *
 * =====================================================================
 */

#include "core/utils/simulationsUtils.h"
#include "core/config/conversion/EnumMapper.h"
#include "core/log/Logger.h"
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include "core/param/paramUtils.h"
#include "core/io/fieldArrayLoader.h"
#include "json/json.h"

namespace ICARION {
namespace utils {


// [REMOVED] init_ions() function - replaced by FullConfig::generate_ions() in config system
// See: src/core/config/types/FullConfig.h and src/core/config/loader/IonLoader.h

// [REMOVED] color_text() and instrument_name() - replaced by EnumMapper in config system

/**
 * @brief Prints a formatted summary of all instrument domains in TRACE.
 *
 * For each domain, the table shows:
 * - Index and instrument name
 * - Whether a precomputed field array (FA_file) is used
 * - Core field parameters (RF/DC)
 * - Load status of precomputed arrays
 *
 * Automatically loads PA_field files if defined.
 * 
 * @param domains Vector of instrument domains (modified in place if PA_field loaded)
 */
void print_domain_summary(std::vector<InstrumentDomain>& domains) {
    using ICARION::log::Logger;

    Logger::main()->info("=== Instrument Domains Loaded ===");
    Logger::main()->info("{:<6} {:<15} {:<8} {:<12} {:<12} {:<20} {:<20}", 
                         "Idx", "Instrument", "PA", "RF_V [V]", "DC_V [V]", "Status", "Grid Info");
    Logger::main()->info("{}", std::string(85, '-'));

    for (auto& dom : domains) {
        std::string FA_file = dom.FA_file.empty() ? "—" : "yes";
        std::string status = "OK";
        std::string grid_info = "—";

        // --- Try loading field array if specified ---
        if (!dom.FA_file.empty()) {
            try {
                dom.fieldArray = load_field_array(dom.FA_file);
                if (!dom.fieldArray.is_valid()) {
                    dom.fieldArrayLoaded = false;
                    status = "⚠ invalid";
                } else {
                    dom.fieldArrayLoaded = true;
                    status = "✓ loaded";

                    grid_info = std::to_string(dom.fieldArray.nx) + "×" +
                                std::to_string(dom.fieldArray.ny) + "×" +
                                std::to_string(dom.fieldArray.nz);
                }
            } catch (const std::exception& e) {
                status = "⚠ error";
                dom.fieldArrayLoaded = false;
                Logger::config()->error("Error loading PA field ({}): {}", dom.FA_file, e.what());
            }
        } else {
            dom.fieldArrayLoaded = false;
        }

        Logger::main()->info("{:<6} {:<15} {:<8} {:<12} {:<12} {:<20} {:<20}",
                             dom.index,
                             ICARION::config::EnumMapper::instrument_to_string(dom.instrument),
                             FA_file,
                             dom.RF.voltage_V,
                             dom.DC.axial_V,
                             status,
                             grid_info);
    }

    Logger::main()->info("{}\n", std::string(85, '='));
}

/**
 * @brief Print final ion positions, velocities, and arrival times.
 *
 * @param[in] result Final simulation results.
 * @param[in] max_nr_ions Maximum number of ions to print (default 100).
 *
 * @note Intended for human-readable console output; not for data export.
 *       Positions are printed in mm, velocities in m/s, and arrival times in µs.
 */
void print_results(const std::vector<IonState>& ions, size_t max_nr_ions = 100) {
    std::ostringstream oss;
    oss << "Simulation complete.\n";
    oss << "Number of ions: " << ions.size() << "\n";

    size_t count = 0;
    for (const auto& ion : ions) {
        if (count++ >= max_nr_ions)
            break;  

        oss << "Ion " << count << " final state:\n";
        oss << "  Position: X = " << ion.pos.x * 1e3 << " mm, Y = " << ion.pos.y * 1e3
            << " mm, Z = " << ion.pos.z * 1e3 << " mm\n";
        oss << "  Velocity: X = " << ion.vel.x << " m/s, Y = " << ion.vel.y
            << " m/s, Z = " << ion.vel.z << " m/s\n";
        oss << "Arrival time: " << std::fixed << std::setprecision(6) << std::max(0.0, ion.t) * 1e6 << " µs\n";
    }

    if (ions.size() > max_nr_ions) {
        oss << "... output truncated. Showing first " << max_nr_ions << " ions.\n";
    }

    std::cout << oss.str();
}

}  // namespace utils
}  // namespace ICARION

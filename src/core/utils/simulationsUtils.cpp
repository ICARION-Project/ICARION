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
#include <iomanip>
#include <iostream>
#include <vector>
#include <string>
#include "core/param/paramUtils.h"
#include "core/io/fieldArrayLoader.h"
#include "json/json.h"
#include "core/physics/reactions/reactionUtils.h"

namespace ICARION {
namespace utils {

/**
 * @brief Initialize ion states from a JSON ion cloud file.
 *
 * @param[in,out] gParams Simulation parameters (updates num_ions).
 * @param[in] speciesDB Database of species properties for reference.
 * @param[in] domains Vector of instrument domains for domain assignment.
 * @return std::vector<IonState> Fully initialized ion states.
 *
 * @throws std::runtime_error if the ion cloud JSON file cannot be opened.
 *
 * @details
 * Each ion is assigned:
 * - Position and velocity from JSON.
 * - Species properties (mass, charge, mobility, CCS).
 * - Domain-specific environmental parameters (temperature, neutral mass,
 *   particle density, polarizability, gas velocity) based on position.
 * - Birth time (if provided) and "born" flag.
 */

std::vector<IonState> init_ions(GlobalParams& gParams,
                                const ICARION::io::SpeciesDatabase& speciesDB,
                                const std::vector<InstrumentDomain>& domains) {
    std::string   filename = gParams.ion_cloud_file;
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open ion cloud JSON file: " + filename);
    }

    Json::Value root;
    file >> root;
    file.close();

    std::vector<IonState> ions;
    const Json::Value     ionsArray = root["ions"];
    ions.reserve(ionsArray.size());

    for (const auto& ionNode : ionsArray) {
        std::string species_name = ionNode["species"].asString();
        
        // NEW: Use SpeciesDatabase::has() and get()
        if (!speciesDB.has(species_name)) {
            std::cerr << "Warning: Species '" << species_name << "' not found in database, skipping ion\n";
            continue;
        }

        const ICARION::io::Species& sp = speciesDB.get(species_name);

        IonState ion;
        ion.species_id              = sp.id;
        ion.mass_kg                 = sp.mass_kg;
        ion.reduced_mobility_cm2_Vs = sp.mobility_m2Vs * 1e4;  // Convert m²/Vs to cm²/Vs
        ion.ion_charge_C            = sp.charge_C;
        ion.CCS_m2                  = sp.CCS_m2;
        ion.t = 0;
        ion.pos = Vec3(ionNode["pos"][0].asDouble(), ionNode["pos"][1].asDouble(),
                       ionNode["pos"][2].asDouble());

        ion.vel = Vec3(ionNode["vel"][0].asDouble(), ionNode["vel"][1].asDouble(),
                       ionNode["vel"][2].asDouble());

        // optional birth time
        if (ionNode.isMember("birth_time")) {
            ion.birth_time_s = ionNode["birth_time"].asDouble();
        } else {
            ion.birth_time_s = 0.0;
        }
        ion.born = (ion.birth_time_s <= 0.0);

        // Assign domain properties
    ion.current_domain_index = -1; // default: outside all domains
    // Debug prints removed to avoid excessive log spam during automated runs.
    // If you need them, re-enable manually for local debugging.
        for (size_t i = 0; i < domains.size(); ++i) {
            bool inside = isInsideDomain(domains[i], ion.pos);
            if (inside) {
                ion.current_domain_index = static_cast<int>(i);
                const InstrumentDomain& dom = domains[i];  

                ion.domain_neutral_mass_kg           = dom.env.neutral_mass_kg;
                ion.domain_particle_density_m3       = dom.env.particle_density_m_3;
                ion.domain_neutral_polarizability_m3 = dom.env.neutral_polarizability_m3;
                ion.domain_temperature_K             = dom.env.temperature_K;
                ion.domain_gas_velocity_m_s          = dom.env.gas_velocity_m_s;  
                break;
            }
        }
        ion.validate();
        ions.push_back(std::move(ion));
    }
    gParams.num_ions = ions.size();

    return ions;
}

// optional: simple colored output
static std::string color_text(const std::string& text, const std::string& color) {
    static const std::unordered_map<std::string, std::string> c = {
        {"red", "\033[31m"}, {"green", "\033[32m"}, {"yellow", "\033[33m"}, {"reset", "\033[0m"}
    };
    auto it = c.find(color);
    if (it == c.end()) return text;
    return it->second + text + c.at("reset");
}

// map enum → string
static std::string instrument_name(Instrument instr) {
    switch (instr) {
        case Instrument::LQIT:         return "LQIT";
        case Instrument::IMS:          return "IMS";
        case Instrument::Orbitrap:     return "Orbitrap";
        case Instrument::QuadrupoleRF: return "Quadrupole";
        case Instrument::TOF:          return "TOF";
        default:                       return "Unknown";
    }
}

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

    std::cout << "\n=== Instrument Domains Loaded ===\n";
    std::cout << std::left << std::setw(6) << "Idx"
              << std::setw(15) << "Instrument"
              << std::setw(8) << "PA"
              << std::setw(12) << "RF_V [V]"
              << std::setw(12) << "DC_V [V]"
              << std::setw(20) << "Status" 
              << std::setw(20) << "Grid Info" << "\n";
    std::cout << std::string(85, '-') << "\n";

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
                    status = color_text("⚠ invalid", "yellow");
                } else {
                    dom.fieldArrayLoaded = true;
                    status = color_text("✓ loaded", "green");

                    grid_info = std::to_string(dom.fieldArray.nx) + "×" +
                                std::to_string(dom.fieldArray.ny) + "×" +
                                std::to_string(dom.fieldArray.nz);
                }
            } catch (const std::exception& e) {
                status = color_text("⚠ error", "yellow");
                dom.fieldArrayLoaded = false;
                std::cerr << "Error loading PA field (" << dom.FA_file << "): " << e.what() << "\n";
            }
        } else {
            dom.fieldArrayLoaded = false;
        }

        std::cout << std::left << std::setw(6)  << dom.index
                  << std::setw(15) << instrument_name(dom.instrument)
                  << std::setw(8)  << FA_file
                  << std::setw(12) << dom.RF.voltage_V
                  << std::setw(12) << dom.DC.axial_V
                  << std::setw(20) << status
                  << std::setw(20) << grid_info
                  << "\n";
    }

    std::cout << std::string(85, '=') << "\n\n";
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
void print_results(const SimulationResult& result, size_t max_nr_ions = 100) {
    std::ostringstream oss;
    oss << "Simulation complete.\n";
    oss << "Number of ions: " << result.ions.size() << "\n";

    size_t count = 0;
    for (const auto& ion : result.ions) {
        if (count++ >= max_nr_ions)
            break;  

        oss << "Ion " << count << " final state:\n";
        oss << "  Position: X = " << ion.pos.x * 1e3 << " mm, Y = " << ion.pos.y * 1e3
            << " mm, Z = " << ion.pos.z * 1e3 << " mm\n";
        oss << "  Velocity: X = " << ion.vel.x << " m/s, Y = " << ion.vel.y
            << " m/s, Z = " << ion.vel.z << " m/s\n";
        oss << "Arrival time: " << std::fixed << std::setprecision(6) << std::max(0.0, ion.t) * 1e6 << " µs\n";
    }

    if (result.ions.size() > max_nr_ions) {
        oss << "... output truncated. Showing first " << max_nr_ions << " ions.\n";
    }

    std::cout << oss.str();
}

}  // namespace utils
}  // namespace ICARION

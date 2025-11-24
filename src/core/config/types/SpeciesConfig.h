// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_SPECIES_CONFIG_H
#define ICARION_CONFIG_SPECIES_CONFIG_H

#include "../validation/ValidationResult.h"
#include "../../../utils/constants.h"
#include <string>
#include <optional>
#include <unordered_map>

namespace ICARION::config {

/**
 * @brief Physical properties of a single species
 * 
 * Supports ions, neutrals, and clusters with simple constant properties.
 * No temperature-dependent behavior (reserved for v2.0).
 */
struct SpeciesProperties {
    // === Required fields ===
    std::string id;                         ///< Unique identifier (e.g., "H3O+", "He")
    double mass_amu;                        ///< Molecular mass [amu]
    int charge;                             ///< Charge state (-1, 0, +1, +2, ...)
    
    // === Optional fields (ions only) ===
    std::optional<double> mobility_cm2Vs;   ///< Reduced mobility [cm²/(V·s)] at STP
    std::optional<double> CCS_A2;           ///< Collision cross-section [Ų]
    
    // === Optional fields (neutrals only) ===
    std::optional<double> polarizability_A3; ///< Electric polarizability [ų]
    
    // === Optional metadata ===
    std::optional<std::string> name;        ///< Human-readable name
    std::optional<std::string> geometry_file; ///< Path to molecular geometry JSON
    
    // === Reference conditions (optional) ===
    std::optional<double> reference_temperature_K; ///< Temperature for mobility/CCS [K]
    std::optional<double> reference_pressure_Pa;   ///< Pressure for mobility [Pa]
    std::optional<std::string> ccs_method;         ///< CCS determination method
    
    // === Derived quantities (computed after load) ===
    double mass_kg = 0.0;                   ///< Mass [kg] (computed)
    double charge_C = 0.0;                  ///< Charge [C] (computed)
    double mobility_m2Vs = 0.0;             ///< Mobility [m²/(V·s)] (computed)
    double CCS_m2 = 0.0;                    ///< CCS [m²] (computed)
    double polarizability_m3 = 0.0;         ///< Polarizability [m³] (computed)
    
    /**
     * @brief Convert user-friendly units to SI
     * 
     * Converts:
     * - mass_amu → mass_kg
     * - charge (elementary) → charge_C
     * - mobility_cm2Vs → mobility_m2Vs
     * - CCS_A2 → CCS_m2
     * - polarizability_A3 → polarizability_m3
     * 
     * Must be called after loading from JSON.
     */
    void convert_to_SI() {
        // Use global constants from constants.h
        // Mass
        mass_kg = mass_amu * AMU_TO_KG;
        
        // Charge
        charge_C = charge * ELEM_CHARGE_C;
        
        // Mobility (if present)
        if (mobility_cm2Vs) {
            mobility_m2Vs = *mobility_cm2Vs * CM2_TO_M2;
        }
        
        // CCS (if present)
        if (CCS_A2) {
            CCS_m2 = *CCS_A2 * ANGSTROM2_TO_M2;
        }

        // Polarizability (if present)
        if (polarizability_A3) {
            polarizability_m3 = *polarizability_A3 * ANGSTROM3_TO_M3;
        }
    }
    
    /**
     * @brief Validate species parameters
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        if (id.empty()) {
            result.add_error("Species ID cannot be empty");
        }
        
        if (mass_amu <= 0.0) {
            result.add_error("Species '" + id + "': mass_amu must be positive");
        }
        
        // Validate ion properties
        if (charge != 0) {
            if (!mobility_cm2Vs && !CCS_A2) {
                result.add_warning("Ion '" + id + "': Neither mobility nor CCS specified");
            }
            if (polarizability_A3) {
                result.add_warning("Ion '" + id + "': polarizability specified but ions typically use CCS");
            }
        }
        
        // Validate neutral properties
        if (charge == 0) {
            if (mobility_cm2Vs) {
                result.add_warning("Neutral '" + id + "': mobility specified for uncharged species");
            }
        }
        
        // Validate physical ranges
        if (mobility_cm2Vs && *mobility_cm2Vs <= 0.0) {
            result.add_error("Species '" + id + "': mobility must be positive");
        }
        if (CCS_A2 && *CCS_A2 <= 0.0) {
            result.add_error("Species '" + id + "': CCS must be positive");
        }
        if (polarizability_A3 && *polarizability_A3 < 0.0) {
            result.add_error("Species '" + id + "': polarizability cannot be negative");
        }

        return result;
    }
};

/**
 * @brief Complete species database
 * 
 * Dictionary-style storage for O(1) lookup by species ID.
 */
struct SpeciesDatabase {
    std::unordered_map<std::string, SpeciesProperties> species;
    
    /**
     * @brief Get species by ID
     * 
     * @throws std::runtime_error if species not found
     */
    const SpeciesProperties& get(const std::string& id) const {
        auto it = species.find(id);
        if (it == species.end()) {
            throw std::runtime_error("Species '" + id + "' not found in database");
        }
        return it->second;
    }
    
    /**
     * @brief Check if species exists
     */
    bool has(const std::string& id) const {
        return species.find(id) != species.end();
    }
    
    /**
     * @brief Get number of species
     */
    size_t size() const { 
        return species.size(); 
    }
    
    /**
     * @brief Validate database
     * 
     * Checks:
     * - No duplicate IDs (implicit via map structure)
     * - All species have valid properties
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        if (species.empty()) {
            result.add_warning("Species database is empty");
        }
        
        for (const auto& [id, props] : species) {
            result.merge(props.validate());
        }
        
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_SPECIES_CONFIG_H

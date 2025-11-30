// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <string>
#include <unordered_map>
#include <vector>
#include <optional>

#include "utils/constants.h"

namespace ICARION {
namespace io {

/**
 * @brief Species database entry
 * 
 * Contains all physical properties needed for trajectory simulation
 * and optional reference to detailed molecular geometry.
 */
struct Species {
    std::string id;                      ///< Unique identifier (e.g., "H3O+", "Pentanal")
    std::string name;                    ///< Human-readable name
    double mass_u;                       ///< Mass in atomic mass units
    double mass_kg;                      ///< Mass in kilograms (derived)
    double charge_e;                     ///< Charge in elementary charges
    double charge_C;                     ///< Charge in Coulombs (derived)
    double CCS_m2;                       ///< Collision cross-section [m²]
    double mobility_m2Vs;                ///< Ion mobility [m²/(V·s)] (optional)
    std::optional<std::string> geometry_file;  ///< Path to molecular geometry JSON
    
    // Derived quantities
    double reduced_mass_kg;              ///< Reduced mass with background gas [kg]
    
    /**
     * @brief Calculate mobility from CCS using Mason-Schamp equation
     * 
     * @param temperature_K Gas temperature [K]
     * @param neutral_mass_kg Mass of background gas [kg]
     * @return Ion mobility [m²/(V·s)]
     */
    double calculate_mobility(double temperature_K, double neutral_mass_kg) const;
    
    /**
     * @brief Calculate CCS from mobility using Mason-Schamp equation
     * 
     * @param temperature_K Gas temperature [K]
     * @param neutral_mass_kg Mass of background gas [kg]
     * @return Collision cross-section [m²]
     */
    double calculate_CCS(double temperature_K, double neutral_mass_kg) const;
    
    /**
     * @brief Validate species parameters
     * 
     * @throws std::runtime_error if parameters are non-physical
     */
    void validate() const;
};

/**
 * @brief Species database with fast lookup
 * 
 * Provides O(1) lookup by species ID and maintains ordering.
 */
class SpeciesDatabase {
public:
    SpeciesDatabase() = default;
    
    /**
     * @brief Add species to database
     * 
     * @param species Species to add
     * @throws std::runtime_error if species ID already exists
     */
    void add(const Species& species);
    
    /**
     * @brief Get species by ID
     * 
     * @param id Species identifier
     * @return Reference to species
     * @throws std::runtime_error if species not found
     */
    const Species& get(const std::string& id) const;
    
    /**
     * @brief Check if species exists
     * 
     * @param id Species identifier
     * @return true if species exists
     */
    bool has(const std::string& id) const;
    
    /**
     * @brief Get all species
     * 
     * @return Vector of all species in insertion order
     */
    const std::vector<Species>& all() const { return species_list_; }
    
    /**
     * @brief Get number of species
     * 
     * @return Number of species in database
     */
    size_t size() const { return species_list_.size(); }
    
    /**
     * @brief Calculate derived quantities for all species
     * 
     * @param temperature_K Gas temperature [K]
     * @param neutral_mass_kg Mass of background gas [kg]
     * 
     * Calculates missing mobility or CCS values and reduced masses
     * for all species in the database.
     */
    void calculate_derived_quantities(double temperature_K, double neutral_mass_kg);

private:
    std::vector<Species> species_list_;                    ///< Ordered list of species
    std::unordered_map<std::string, size_t> id_to_index_;  ///< Fast ID lookup
};

/**
 * @brief Load species database from JSON file
 * 
 * @param filepath Path to species JSON file
 * @return Species database
 * @throws std::runtime_error on file access errors or validation failures
 * 
 * Expected JSON format:
 * ```json
 * {
 *   "species": [
 *     {
 *       "id": "H3O+",
 *       "name": "Hydronium Ion",
 *       "mass_u": 19.02,
 *       "charge_e": +1,
 *       "CCS_m2": 7.8e-19,
 *       "geometry_file": "molecules/H3O+.json"
 *     }
 *   ]
 * }
 * ```
 * 
 * Required fields: id, mass_u, charge_e
 * Optional fields: name, CCS_m2, mobility_m2Vs, geometry_file
 * 
 * If both CCS and mobility are missing, throws error.
 * If only one is provided, the other is calculated.
 */
SpeciesDatabase load_species(const std::string& filepath);

/**
 * @brief Load species database with environment parameters
 * 
 * @param filepath Path to species JSON file
 * @param temperature_K Gas temperature [K] for mobility/CCS calculation
 * @param neutral_mass_kg Mass of background gas [kg]
 * @return Species database with calculated derived quantities
 * @throws std::runtime_error on file access errors or validation failures
 * 
 * Convenience function that loads species and automatically calculates
 * mobility from CCS (or vice versa) and reduced masses.
 */
SpeciesDatabase load_species(const std::string& filepath, 
                             double temperature_K, 
                             double neutral_mass_kg);

}  // namespace io
}  // namespace ICARION

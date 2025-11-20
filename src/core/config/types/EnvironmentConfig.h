// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_ENVIRONMENT_CONFIG_H
#define ICARION_CONFIG_ENVIRONMENT_CONFIG_H

#include "core/utils/mathUtils.h"
#include "utils/constants.h"
#include <string>
#include <cmath>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief Environment configuration (gas conditions)
 * 
 * Input parameters + derived thermodynamic quantities.
 * Clean separation: input params are loaded from JSON, derived quantities
 * are computed via compute_derived_properties().
 * 
 * Future: Support time-varying temperature/pressure for gradient simulations.
 */
struct EnvironmentConfig {
    // === Input parameters ===
    double pressure_Pa = 101325.0;          ///< Pressure [Pa] (default: 1 atm)
    double temperature_K = 300.0;           ///< Temperature [K] (default: room temp)
    std::string gas_species = "He";         ///< Buffer gas species ID
    Vec3 gas_velocity_m_s = {0.0, 0.0, 0.0}; ///< Bulk gas flow velocity [m/s]
    
    // === Derived quantities (computed from input) ===
    double particle_density_m_3 = 0.0;      ///< Number density n = P/(kB·T) [m⁻³]
    double mean_thermal_velocity_m_s = 0.0; ///< v_th = sqrt(8kBT/πm) [m/s]
    double mean_free_path_m = 0.0;          ///< λ = 1/(√2 π d² n) [m]
    
    // Gas properties (set by compute_derived_properties)
    double gas_mass_kg = 0.0;               ///< Molecular mass [kg]
    double gas_polarizability_m3 = 0.0;     ///< Polarizability [m³]
    double gas_radius_m = 0.0;              ///< Hard-sphere radius [m]
    
    /**
     * @brief Compute derived thermodynamic quantities
     * 
     * Must be called after loading from JSON.
     * Uses constants from utils/constants.h for gas properties.
     */
    void compute_derived_properties() {
        // Lookup gas properties from constants.h
        if (gas_species == "N2") {
            gas_mass_kg = MOLAR_MASS_N2_KG;
            gas_radius_m = RADIUS_N2_M;
            gas_polarizability_m3 = POLARIZABILITY_N2_SI;
        } else if (gas_species == "He") {
            gas_mass_kg = MOLAR_MASS_HE_KG;
            gas_radius_m = RADIUS_HE_M;
            gas_polarizability_m3 = POLARIZABILITY_HE_SI;
        } else if (gas_species == "Ar") {
            gas_mass_kg = MOLAR_MASS_AR_KG;
            gas_radius_m = RADIUS_AR_M;
            gas_polarizability_m3 = POLARIZABILITY_AR_SI;
        } else if (gas_species == "O2") {
            gas_mass_kg = MOLAR_MASS_O2_KG;
            gas_radius_m = RADIUS_O2_M;
            gas_polarizability_m3 = POLARIZABILITY_O2_SI;
        } else if (gas_species == "CO2") {
            gas_mass_kg = MOLAR_MASS_CO2_KG;
            gas_radius_m = RADIUS_CO2_M;
            gas_polarizability_m3 = POLARIZABILITY_CO2_SI;
        } else if (gas_species == "Ne") {
            gas_mass_kg = MOLAR_MASS_NE_KG;
            gas_radius_m = RADIUS_NE_M;
            gas_polarizability_m3 = POLARIZABILITY_NE_SI;
        } else if (gas_species == "Air") {
            // Approximate as mean of N2, O2, Ar
            gas_mass_kg = 0.78084 * MOLAR_MASS_N2_KG +
                           0.20946 * MOLAR_MASS_O2_KG +
                           0.00934 * MOLAR_MASS_AR_KG;
            gas_radius_m = 0.78084 * RADIUS_N2_M +
                           0.20946 * RADIUS_O2_M +
                           0.00934 * RADIUS_AR_M;
            gas_polarizability_m3 = 0.78084 * POLARIZABILITY_N2_SI +
                                    0.20946 * POLARIZABILITY_O2_SI +
                                    0.00934 * POLARIZABILITY_AR_SI;
        } else if (gas_species == "Default" || gas_species.empty()) 
        {
            // Default to He if unknown
            gas_mass_kg = MOLAR_MASS_HE_KG;
            gas_radius_m = RADIUS_HE_M;
            gas_polarizability_m3 = POLARIZABILITY_HE_SI;
        }
        
        // Compute derived quantities
        particle_density_m_3 = pressure_Pa / (BOLTZMANN_CONSTANT * temperature_K);
        
        mean_thermal_velocity_m_s = std::sqrt(
            8.0 * BOLTZMANN_CONSTANT * temperature_K / (M_PI * gas_mass_kg)
        );
        
        mean_free_path_m = 1.0 / (
            std::sqrt(2.0) * M_PI * gas_radius_m * gas_radius_m * particle_density_m_3
        );
    }
    
    /**
     * @brief Validate environment parameters
     * 
     * @throws std::runtime_error if invalid
     */
    void validate() const {
        if (pressure_Pa <= 0.0) {
            throw std::runtime_error("Pressure must be positive");
        }
        if (temperature_K <= 0.0) {
            throw std::runtime_error("Temperature must be positive");
        }
        if (temperature_K < 1.0) {
            throw std::runtime_error("Temperature unrealistically low (< 1 K)");
        }
        if (temperature_K > 10000.0) {
            throw std::runtime_error("Temperature unrealistically high (> 10000 K)");
        }
        if (gas_species.empty()) {
            throw std::runtime_error("Gas species cannot be empty");
        }
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_ENVIRONMENT_CONFIG_H

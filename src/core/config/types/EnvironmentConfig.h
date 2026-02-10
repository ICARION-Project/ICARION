// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_ENVIRONMENT_CONFIG_H
#define ICARION_CONFIG_ENVIRONMENT_CONFIG_H

#include "core/utils/mathUtils.h"
#include "utils/constants.h"
#include "WaveformConfig.h"
#include "../validation/ValidationResult.h"
#include <string>
#include <cmath>
#include <stdexcept>
#include <vector>
#include <unordered_map>

namespace ICARION::config {

/**
 * @brief Single gas component for mixtures
 */
struct GasMixtureComponent {
    std::string species;          ///< Species ID (e.g., "N2")
    double mole_fraction = 0.0;   ///< Mole fraction (0..1)
    double cross_section_m2 = -1.0;      ///< Optional override for σ [m²], <0 = use default/ion CCS
    double polarizability_m3 = -1.0;     ///< Optional override for α [m³], <0 = use default
    bool participates_in_collisions = true; ///< If false, ignored for collision rates
    bool participates_in_reactions = true;  ///< If false, ignored for reaction concentrations

    // Derived
    double mass_kg = 0.0;
    double radius_m = 0.0;
    double density_m3 = 0.0;
    double polarizability_m3_derived = 0.0;
};

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
    ValueOrWaveform pressure_Pa_waveform{101325.0}; ///< Optional pressure waveform (evaluated at runtime)
    double temperature_K = 300.0;           ///< Temperature [K] (default: room temp)
    std::string gas_species = "He";         ///< Buffer gas species ID
    Vec3 gas_velocity_m_s = {0.0, 0.0, 0.0}; ///< Bulk gas flow velocity [m/s]
    std::vector<GasMixtureComponent> gas_mixture; ///< Optional mixture (if non-empty, overrides gas_species)
    
    // === Derived quantities (computed from input) ===
    double particle_density_m_3 = 0.0;      ///< Number density n = P/(kB·T) [m⁻³]
    double mean_thermal_velocity_m_s = 0.0; ///< v_th = sqrt(8kBT/πm) [m/s]
    double mean_free_path_m = 0.0;          ///< λ = 1/(√2 π d² n) [m]
    
    // Gas properties (set by compute_derived_properties)
    double gas_mass_kg = 0.0;               ///< Molecular mass [kg]
    double gas_polarizability_m3 = 0.0;     ///< Polarizability [m³]
    double gas_radius_m = 0.0;              ///< Hard-sphere radius [m]

    bool has_dynamic_pressure() const {
        return pressure_Pa_waveform.is_time_varying();
    }

    void update_time_dependent(double t_s, const std::map<std::string, Waveform>& waveform_library = {}) {
        if (!has_dynamic_pressure()) {
            return;
        }

        pressure_Pa = pressure_Pa_waveform.evaluate(t_s, waveform_library);
        if (pressure_Pa <= 0.0) {
            throw std::runtime_error("EnvironmentConfig: pressure waveform evaluated to non-positive value");
        }

        compute_derived_properties();
    }
    
    /**
     * @brief Compute derived thermodynamic quantities
     * 
     * Must be called after loading from JSON.
     * Uses constants from utils/constants.h for gas properties.
     */
    void compute_derived_properties() {
        // Total number density from ideal gas law
        particle_density_m_3 = pressure_Pa / (BOLTZMANN_CONSTANT * temperature_K);

        // Mixture handling
        if (!gas_mixture.empty()) {
            double frac_sum = 0.0;
            for (const auto& comp : gas_mixture) {
                frac_sum += comp.mole_fraction;
            }
            // Normalize fractions if they don't sum to 1 (avoid divide by zero)
            double norm = (frac_sum > 0.0) ? frac_sum : 1.0;

            // Derived mixture properties
            gas_mass_kg = 0.0;
            gas_radius_m = 0.0;
            gas_polarizability_m3 = 0.0;

            for (auto& comp : gas_mixture) {
                double f = comp.mole_fraction / norm;
                comp.density_m3 = f * particle_density_m_3;

                // Lookup gas properties from constants.h
                if (comp.species == "N2") {
                    comp.mass_kg = MOLAR_MASS_N2_KG;
                    comp.radius_m = RADIUS_N2_M;
                    comp.polarizability_m3_derived = POLARIZABILITY_N2_SI;
                } else if (comp.species == "He") {
                    comp.mass_kg = MOLAR_MASS_HE_KG;
                    comp.radius_m = RADIUS_HE_M;
                    comp.polarizability_m3_derived = POLARIZABILITY_HE_SI;
                } else if (comp.species == "Ar") {
                    comp.mass_kg = MOLAR_MASS_AR_KG;
                    comp.radius_m = RADIUS_AR_M;
                    comp.polarizability_m3_derived = POLARIZABILITY_AR_SI;
                } else if (comp.species == "O2") {
                    comp.mass_kg = MOLAR_MASS_O2_KG;
                    comp.radius_m = RADIUS_O2_M;
                    comp.polarizability_m3_derived = POLARIZABILITY_O2_SI;
                } else if (comp.species == "CO2") {
                    comp.mass_kg = MOLAR_MASS_CO2_KG;
                    comp.radius_m = RADIUS_CO2_M;
                    comp.polarizability_m3_derived = POLARIZABILITY_CO2_SI;
                } else if (comp.species == "Ne") {
                    comp.mass_kg = MOLAR_MASS_NE_KG;
                    comp.radius_m = RADIUS_NE_M;
                    comp.polarizability_m3_derived = POLARIZABILITY_NE_SI;
                } else if (comp.species == "H2O") {
                    comp.mass_kg = MOLAR_MASS_H2O_KG;
                    comp.radius_m = RADIUS_H2O_M;
                    comp.polarizability_m3_derived = POLARIZABILITY_H2O_SI;
                } else {
                    // Fallback to He-like properties
                    comp.mass_kg = MOLAR_MASS_HE_KG;
                    comp.radius_m = RADIUS_HE_M;
                    comp.polarizability_m3_derived = POLARIZABILITY_HE_SI;
                }

                if (comp.polarizability_m3 > 0.0) {
                    comp.polarizability_m3_derived = comp.polarizability_m3;
                }

                gas_mass_kg += f * comp.mass_kg;
                gas_radius_m += f * comp.radius_m;
                gas_polarizability_m3 += f * comp.polarizability_m3_derived;
            }

            // Mean thermal velocity using mixture-averaged mass
            mean_thermal_velocity_m_s = std::sqrt(
                8.0 * BOLTZMANN_CONSTANT * temperature_K / (M_PI * gas_mass_kg)
            );

            // Mean free path approximate with weighted radius (rough)
            mean_free_path_m = 1.0 / (
                std::sqrt(2.0) * M_PI * gas_radius_m * gas_radius_m * particle_density_m_3
            );

            return;
        }

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
        mean_thermal_velocity_m_s = std::sqrt(
            8.0 * BOLTZMANN_CONSTANT * temperature_K / (M_PI * gas_mass_kg)
        );
        
        mean_free_path_m = 1.0 / (
            std::sqrt(2.0) * M_PI * gas_radius_m * gas_radius_m * particle_density_m_3
        );
    }
    
    /**
     * @brief Validate environment parameters
     */
    ValidationResult validate() const {
        ValidationResult result;

        if (!pressure_Pa_waveform.is_valid()) {
            result.add_error("Pressure waveform configuration is invalid");
        }
        
        if (pressure_Pa <= 0.0) {
            result.add_error("Pressure must be positive");
        }
        if (temperature_K <= 0.0) {
            result.add_error("Temperature must be positive");
        }
        if (temperature_K < 1.0) {
            result.add_warning("Temperature unrealistically low (< 1 K)");
        }
        if (temperature_K > 10000.0) {
            result.add_warning("Temperature unrealistically high (> 10000 K)");
        }
        if (gas_species.empty()) {
            result.add_error("Gas species cannot be empty");
        }
        double frac_sum = 0.0;
        for (const auto& comp : gas_mixture) {
            if (comp.species.empty()) {
                result.add_error("Gas mixture component species cannot be empty");
            }
            if (comp.mole_fraction < 0.0) {
                result.add_error("Gas mixture component mole_fraction must be >= 0");
            }
            if (comp.cross_section_m2 >= 0.0 && comp.cross_section_m2 == 0.0) {
                result.add_warning("Gas mixture component cross_section_m2 is zero");
            }
            frac_sum += comp.mole_fraction;
        }
        if (!gas_mixture.empty()) {
            if (frac_sum <= 0.0) {
                result.add_error("Gas mixture mole fractions must sum to > 0");
            } else if (std::abs(frac_sum - 1.0) > 1e-3) {
                result.add_warning("Gas mixture mole fractions do not sum to 1 (sum=" + std::to_string(frac_sum) + ")");
            }
        }
        
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_ENVIRONMENT_CONFIG_H

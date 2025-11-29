// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file BoundaryConfig.h
 * @brief Configuration for domain boundary actions
 * 
 * Defines what happens when ion hits a domain boundary:
 * - Absorption: ion removed (detector, lossy walls)
 * - Specular Reflection: elastic bounce
 * - Diffuse Reflection: cosine-weighted random
 * - Thermal Reflection: Maxwell-Boltzmann re-emission
 */

#pragma once

#include "../validation/ValidationResult.h"
#include <string>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief Boundary action type
 */
enum class BoundaryActionType {
    Absorption,         ///< Ion absorbed at boundary (removed from simulation)
    SpecularReflection, ///< Mirror reflection: v' = v - 2(v·n)n
    DiffuseReflection,  ///< Cosine-weighted random direction + partial thermalization
    ThermalReflection   ///< Full Maxwell-Boltzmann re-emission at wall temperature
};

/**
 * @brief Configuration for boundary actions
 */
struct BoundaryConfig {
    BoundaryActionType type = BoundaryActionType::Absorption;  ///< Action type
    
    /**
     * Energy accommodation coefficient [0,1]
     * - 0 = fully specular (elastic)
     * - 1 = fully thermal (complete energy exchange with wall)
     * 
     * Used only for DiffuseReflection:
     * v_out = (1-α)*v_in + α*v_thermal
     */
    double accommodation_coeff = 1.0;
    
    /**
     * Wall temperature [K]
     * Used for DiffuseReflection and ThermalReflection to compute thermal speed.
     * If not specified, uses domain environment temperature.
     */
    double temperature_K = 300.0;
    
    /**
     * @brief Convert string to BoundaryActionType
     */
    static BoundaryActionType parse_type(const std::string& type_str) {
        if (type_str == "absorption" || type_str == "Absorption") {
            return BoundaryActionType::Absorption;
        } else if (type_str == "specular" || type_str == "SpecularReflection") {
            return BoundaryActionType::SpecularReflection;
        } else if (type_str == "diffuse" || type_str == "DiffuseReflection") {
            return BoundaryActionType::DiffuseReflection;
        } else if (type_str == "thermal" || type_str == "ThermalReflection") {
            return BoundaryActionType::ThermalReflection;
        } else {
            throw std::invalid_argument("Unknown boundary action type: " + type_str);
        }
    }
    
    /**
     * @brief Convert BoundaryActionType to string
     */
    static std::string type_to_string(BoundaryActionType type) {
        switch (type) {
            case BoundaryActionType::Absorption:         return "Absorption";
            case BoundaryActionType::SpecularReflection: return "SpecularReflection";
            case BoundaryActionType::DiffuseReflection:  return "DiffuseReflection";
            case BoundaryActionType::ThermalReflection:  return "ThermalReflection";
            default: return "Unknown";
        }
    }
    
    /**
     * @brief Validate boundary configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        if (accommodation_coeff < 0.0 || accommodation_coeff > 1.0) {
            result.add_error("Accommodation coefficient must be in [0,1], got " + 
                           std::to_string(accommodation_coeff));
        }
        
        if (temperature_K <= 0.0) {
            result.add_error("Wall temperature must be positive, got " + 
                           std::to_string(temperature_K) + " K");
        }
        
        if (temperature_K < 1.0) {
            result.add_warning("Wall temperature unusually low: " + 
                             std::to_string(temperature_K) + " K. Is this intentional?");
        }
        
        if (temperature_K > 10000.0) {
            result.add_warning("Wall temperature unusually high: " + 
                             std::to_string(temperature_K) + " K. Is this intentional?");
        }
        
        return result;
    }
};

}  // namespace ICARION::config

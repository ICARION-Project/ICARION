// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   ------------------------------------------------
 *   Modular framework for simulating ion trajectories in custom
 *   electric fields and background gas environments.
 *
 *   @file       configurationLoader.cpp
 *   @brief      Configuration loading implementation
 *
 *   @details
 *   Loads and validates simulation configuration from JSON files.  
 *
 *   @date       2025-11-09
 *   @version    1.0.0
 *   @authors    ICARION Development Team
 *               
 *
 * =====================================================================
 */

#include "core/io/configurationLoader.h"
#include "core/log/Logger.h"
#include <fstream>
#include <stdexcept>
#include <iostream>

namespace ICARION {
namespace core {

Json::Value ConfigurationLoader::loadFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open configuration file: " + filename);
    }
    
    Json::Value config;
    Json::Reader reader;
    
    if (!reader.parse(file, config)) {
        throw std::runtime_error("Failed to parse JSON: " + reader.getFormattedErrorMessages());
    }
    
    return config;
}

bool ConfigurationLoader::validate(const Json::Value& config) {
    return validateConfig(config);
}

bool ConfigurationLoader::validateConfig(const Json::Value& config) const {
    using ICARION::log::Logger;
    
    if (!config.isMember("simulation") || !config["simulation"].isObject()) {
        Logger::config()->error("Validation error: missing 'simulation' object (v1.0 schema)");
        return false;
    }

    const auto& sim = config["simulation"];
    if (!sim.isMember("timestep_ns") || !sim["timestep_ns"].isNumeric()) {
        Logger::config()->error("Validation error: simulation.timestep_ns must be a numeric value");
        return false;
    }
    if (!sim.isMember("max_time_ns") || !sim["max_time_ns"].isNumeric()) {
        Logger::config()->error("Validation error: simulation.max_time_ns must be a numeric value");
        return false;
    }
    if (!config.isMember("instrument")) {
        Logger::config()->error("Validation error: missing 'instrument' section");
        return false;
    }
    if (!config.isMember("ions") || !config["ions"].isObject()) {
        Logger::config()->error("Validation error: missing 'ions' object");
        return false;
    }
    const auto& ions = config["ions"];
    if (!ions.isMember("species") || !ions["species"].isArray() || ions["species"].empty()) {
        Logger::config()->error("Validation error: ions.species array is required in v1.0 schema");
        return false;
    }

    return true;
}

} // namespace core
} // namespace ICARION

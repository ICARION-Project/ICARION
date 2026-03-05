// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_CONFIG_SIMULATION_LOADER_H
#define ICARION_CONFIG_SIMULATION_LOADER_H

#include "../types/SimulationConfig.h"
#include <json/json.h>

namespace ICARION::config {

/**
 * @brief Loader for simulation configuration section
 */
class SimulationConfigLoader {
public:
    /**
     * @brief Load simulation config from JSON
     * 
     * @param json JSON object for simulation section
     * @return SimulationConfig Parsed configuration
     */
    static SimulationConfig load(const Json::Value& json);
    
private:
    // Helper: Get double with default
    static double get_double(const Json::Value& json, const char* key, double default_val);
    
    // Helper: Get int with default
    static int get_int(const Json::Value& json, const char* key, int default_val);
    
    // Helper: Get bool with default
    static bool get_bool(const Json::Value& json, const char* key, bool default_val);
    
    // Helper: Get string with default
    static std::string get_string(const Json::Value& json, const char* key, const std::string& default_val);
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_SIMULATION_LOADER_H

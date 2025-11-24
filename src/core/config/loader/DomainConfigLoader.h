// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_DOMAIN_LOADER_H
#define ICARION_CONFIG_DOMAIN_LOADER_H

#include "../types/DomainConfig.h"
#include "../conversion/EnumMapper.h"
#include <json/json.h>

namespace ICARION::config {

/**
 * @brief Loader for instrument domain configuration
 * 
 * Handles complex domain parsing with geometry, environment, and fields.
 */
class DomainConfigLoader {
public:
    /**
     * @brief Load single domain from JSON object
     * 
     * @param json JSON object for one domain
     * @param default_integrator Fallback integrator from simulation.integrator (default: "RK4")
     * @return DomainConfig Parsed domain configuration
     */
    static DomainConfig load(const Json::Value& json, const std::string& default_integrator = "RK4");
    
private:
    // Sub-loaders
    static GeometryConfig load_geometry(const Json::Value& json);
    static EnvironmentConfig load_environment(const Json::Value& json);
    static FieldsConfig load_fields(const Json::Value& json);
    
    // Field sub-components (v1.1: now accept waveform library for reference resolution)
    static DCFieldConfig load_dc_fields(const Json::Value& json, const std::map<std::string, Waveform>& waveform_library);
    static RFFieldConfig load_rf_fields(const Json::Value& json, const std::map<std::string, Waveform>& waveform_library);
    static ACFieldConfig load_ac_fields(const Json::Value& json, const std::map<std::string, Waveform>& waveform_library);
    static MagneticFieldConfig load_magnetic_fields(const Json::Value& json);
    
    // Helpers
    static Vec3 load_vec3(const Json::Value& json, const Vec3& default_val = {0,0,0});
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_DOMAIN_LOADER_H

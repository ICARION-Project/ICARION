// SPDX-License-Identifier: Apache-2.0
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
     * @return DomainConfig Parsed domain configuration
     */
    static DomainConfig load(const Json::Value& json);
    
private:
    // Sub-loaders
    static GeometryConfig load_geometry(const Json::Value& json);
    static EnvironmentConfig load_environment(const Json::Value& json);
    static FieldsConfig load_fields(const Json::Value& json);
    
    // Field sub-components
    static DCFieldConfig load_dc_fields(const Json::Value& json);
    static RFFieldConfig load_rf_fields(const Json::Value& json);
    static ACFieldConfig load_ac_fields(const Json::Value& json);
    static MagneticFieldConfig load_magnetic_fields(const Json::Value& json);
    
    // Helpers
    static Vec3 load_vec3(const Json::Value& json, const Vec3& default_val = {0,0,0});
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_DOMAIN_LOADER_H

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
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
     * @param global_waveforms Global waveform library (SSOT, optional)
     * @return DomainConfig Parsed domain configuration
     */
    static DomainConfig load(
        const Json::Value& json, 
        const std::string& default_integrator = "RK4",
        const std::map<std::string, Waveform>& global_waveforms = {}
    );
    
private:
    // Sub-loaders
    static GeometryConfig load_geometry(const Json::Value& json);
    static EnvironmentConfig load_environment(const Json::Value& json, const std::map<std::string, Waveform>& global_waveforms);
    static FieldsConfig load_fields(const Json::Value& json, const std::map<std::string, Waveform>& global_waveforms, const GeometryConfig& geometry, const EnvironmentConfig& environment);
    static BoundaryConfig load_boundary(const Json::Value& json, double domain_temperature_K);
    
    // Field sub-components (v1.1: now accept waveform library for reference resolution + geometry/env for EN_Td conversion)
    static DCFieldConfig load_dc_fields(const Json::Value& json, const std::map<std::string, Waveform>& local_library, const std::map<std::string, Waveform>& global_library, const GeometryConfig& geometry, const EnvironmentConfig& environment);
    static RFFieldConfig load_rf_fields(const Json::Value& json, const std::map<std::string, Waveform>& local_library, const std::map<std::string, Waveform>& global_library);
    static ACFieldConfig load_ac_fields(const Json::Value& json, const std::map<std::string, Waveform>& local_library, const std::map<std::string, Waveform>& global_library);
    static MagneticFieldConfig load_magnetic_fields(const Json::Value& json);
    
    // Helpers
    static Vec3 load_vec3(const Json::Value& json, const Vec3& default_val = {0,0,0});
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_DOMAIN_LOADER_H

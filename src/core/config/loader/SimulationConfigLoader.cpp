// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "SimulationConfigLoader.h"
#include <stdexcept>

namespace ICARION::config {

SimulationConfig SimulationConfigLoader::load(const Json::Value& json) {
    SimulationConfig config;
    
    // Time parameters - REQUIRED fields
    if (!json.isMember("total_time_s") || !json["total_time_s"].isNumeric()) {
        throw std::runtime_error("Missing required field 'total_time_s' in simulation config");
    }
    config.total_time_s = json["total_time_s"].asDouble();
    
    if (!json.isMember("dt_s") || !json["dt_s"].isNumeric()) {
        throw std::runtime_error("Missing required field 'dt_s' in simulation config");
    }
    config.dt_s = json["dt_s"].asDouble();
    
    config.write_interval = get_int(json, "write_interval", 100);
    
    // Integrator
    config.integrator = get_string(json, "integrator", "RK4");
    if (json.isMember("rk45_settings") && json["rk45_settings"].isObject()) {
        const auto& rk = json["rk45_settings"];
        if (rk.isMember("dt_min") && rk["dt_min"].isNumeric()) {
            config.rk45_min_step_s = rk["dt_min"].asDouble();
        }
    }
    // Backward-compatible direct key
    config.rk45_min_step_s = get_double(json, "rk45_min_step_s", config.rk45_min_step_s);
    
    // Execution mode
    config.enable_gpu = get_bool(json, "enable_gpu", false);
    config.enable_openmp = get_bool(json, "enable_openmp", false);
    config.rng_seed = get_int(json, "rng_seed", 42);

    // Numerical safety (optional)
    if (json.isMember("numerical_safety") && json["numerical_safety"].isObject()) {
        const auto& safety = json["numerical_safety"];
        config.enable_safety_logging = get_bool(safety, "enable_logging", config.enable_safety_logging);
        config.verbose_safety = get_bool(safety, "verbose_mode", config.verbose_safety);
        config.safety_checks.enable_nan_checks =
            get_bool(safety, "enable_nan_checks", config.safety_checks.enable_nan_checks);
        config.safety_checks.enable_bounds_checks =
            get_bool(safety, "enable_bounds_checks", config.safety_checks.enable_bounds_checks);
        config.safety_checks.throw_on_violation =
            get_bool(safety, "throw_on_violation", config.safety_checks.throw_on_violation);
        config.safety_checks.attempt_recovery =
            get_bool(safety, "attempt_recovery", config.safety_checks.attempt_recovery);
        config.safety_checks.max_position_m =
            get_double(safety, "max_position_m", config.safety_checks.max_position_m);
        config.safety_checks.max_velocity_ms =
            get_double(safety, "max_velocity_ms", config.safety_checks.max_velocity_ms);
        config.safety_checks.max_acceleration_ms2 =
            get_double(safety, "max_acceleration_ms2", config.safety_checks.max_acceleration_ms2);
    }
    
    // Checkpointing
    config.continue_from = get_string(json, "continue_from", "");
    config.continue_time_s = get_double(json, "continue_time_s", 0.0);
    config.auto_continue_if_active = get_bool(json, "auto_continue_if_active", false);
    
    // Compute derived quantities
    config.compute_derived();
    
    // Validate
    config.validate();
    
    return config;
}

double SimulationConfigLoader::get_double(const Json::Value& json, const char* key, double default_val) {
    if (json.isMember(key) && json[key].isNumeric()) {
        return json[key].asDouble();
    }
    return default_val;
}

int SimulationConfigLoader::get_int(const Json::Value& json, const char* key, int default_val) {
    if (json.isMember(key) && json[key].isInt()) {
        return json[key].asInt();
    }
    return default_val;
}

bool SimulationConfigLoader::get_bool(const Json::Value& json, const char* key, bool default_val) {
    if (json.isMember(key) && json[key].isBool()) {
        return json[key].asBool();
    }
    return default_val;
}

std::string SimulationConfigLoader::get_string(const Json::Value& json, const char* key, 
                                               const std::string& default_val) {
    if (json.isMember(key) && json[key].isString()) {
        return json[key].asString();
    }
    return default_val;
}

} // namespace ICARION::config

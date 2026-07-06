// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "PhysicsConfigLoader.h"

namespace ICARION::config {

PhysicsConfig PhysicsConfigLoader::load(const Json::Value& json) {
    PhysicsConfig config;
    
    // Collision model (required)
    if (json.isMember("collision_model")) {
        std::string model_str = json["collision_model"].asString();
        config.collision_model = EnumMapper::parse_collision_model(model_str);
    } else {
        throw std::runtime_error("Missing required 'collision_model' in physics section");
    }
    
    // Feature flags
    if (json.isMember("enable_reactions") && json["enable_reactions"].isBool()) {
        config.enable_reactions = json["enable_reactions"].asBool();
    }
    
    if (json.isMember("enable_space_charge") && json["enable_space_charge"].isBool()) {
        config.enable_space_charge = json["enable_space_charge"].asBool();
    }
    
    if (json.isMember("enable_space_charge_gpu") && json["enable_space_charge_gpu"].isBool()) {
        config.enable_space_charge_gpu = json["enable_space_charge_gpu"].asBool();
    }
    // Stochastic collision event handling
    if (json.isMember("collision_multi_event_mode") && json["collision_multi_event_mode"].isBool()) {
        config.collision_multi_event_mode = json["collision_multi_event_mode"].asBool();
    }
    if (json.isMember("collision_max_events_per_step") && json["collision_max_events_per_step"].isInt()) {
        config.collision_max_events_per_step = json["collision_max_events_per_step"].asInt();
    }
    if (json.isMember("collision_time_centered") && json["collision_time_centered"].isBool()) {
        config.collision_time_centered = json["collision_time_centered"].asBool();
    }
    if (json.isMember("collision_time_randomized") && json["collision_time_randomized"].isBool()) {
        config.collision_time_randomized = json["collision_time_randomized"].asBool();
    }
    if (json.isMember("collision_subcycles_per_step") && json["collision_subcycles_per_step"].isInt()) {
        config.collision_subcycles_per_step = json["collision_subcycles_per_step"].asInt();
    }
    if (json.isMember("ipm_orientation_mode") && json["ipm_orientation_mode"].isString()) {
        config.ipm_orientation_mode_type = EnumMapper::parse_ipm_orientation_mode(json["ipm_orientation_mode"].asString());
    }
    if (json.isMember("ipm_fixed_orientation_index") && json["ipm_fixed_orientation_index"].isInt()) {
        config.ipm_fixed_orientation_index = json["ipm_fixed_orientation_index"].asInt();
    }
    if (json.isMember("ipm_energy_accommodation")) {
        throw std::runtime_error(
            "physics.ipm_energy_accommodation has been removed "
            "because inelastic accommodation is not implemented in the runtime collision dynamics");
    }
    if (json.isMember("ipm_vrel_log_prefix") && json["ipm_vrel_log_prefix"].isString()) {
        config.ipm_vrel_log_prefix = json["ipm_vrel_log_prefix"].asString();
    }
    if (json.isMember("ipm_momentum_log_prefix") && json["ipm_momentum_log_prefix"].isString()) {
        config.ipm_momentum_log_prefix = json["ipm_momentum_log_prefix"].asString();
    }
    
    // Thermalization
    if (json.isMember("enable_ou_thermalization") && json["enable_ou_thermalization"].isBool()) {
        config.enable_ou_thermalization = json["enable_ou_thermalization"].asBool();
    }
    
    if (json.isMember("force_ou_for_stochastic") && json["force_ou_for_stochastic"].isBool()) {
        config.force_ou_for_stochastic = json["force_ou_for_stochastic"].asBool();
    }
    
    // Validate
    config.validate();

    config.ipm_orientation_mode = EnumMapper::ipm_orientation_mode_to_string(config.ipm_orientation_mode_type);

    return config;
}

} // namespace ICARION::config

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

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
    
    // Thermalization
    if (json.isMember("enable_ou_thermalization") && json["enable_ou_thermalization"].isBool()) {
        config.enable_ou_thermalization = json["enable_ou_thermalization"].asBool();
    }
    
    if (json.isMember("force_ou_for_stochastic") && json["force_ou_for_stochastic"].isBool()) {
        config.force_ou_for_stochastic = json["force_ou_for_stochastic"].asBool();
    }
    
    // Validate
    config.validate();
    
    return config;
}

} // namespace ICARION::config

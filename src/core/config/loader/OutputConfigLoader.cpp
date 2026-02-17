// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include "OutputConfigLoader.h"

namespace ICARION::config {

OutputConfig OutputConfigLoader::load(const Json::Value& json) {
    OutputConfig config;
    
    // Output folder
    if (json.isMember("folder") && json["folder"].isString()) {
        config.folder = json["folder"].asString();
    }
    
    // Trajectory file
    if (json.isMember("trajectory_file") && json["trajectory_file"].isString()) {
        config.trajectory_file = json["trajectory_file"].asString();
    }
    
    // Print progress
    if (json.isMember("print_progress") && json["print_progress"].isBool()) {
        config.print_progress = json["print_progress"].asBool();
    }

    // Optional RAM cap for trajectory buffer
    if (json.isMember("buffer_byte_cap") && json["buffer_byte_cap"].isNumeric()) {
        auto cap = json["buffer_byte_cap"].asInt64();
        if (cap < 0) {
            throw std::runtime_error("output.buffer_byte_cap must be >= 0");
        }
        config.buffer_byte_cap = static_cast<size_t>(cap);
    }
    
    // Validate
    config.validate();
    
    return config;
}

} // namespace ICARION::config

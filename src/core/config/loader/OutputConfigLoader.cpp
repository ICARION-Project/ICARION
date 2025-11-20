// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

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
    
    // Validate
    config.validate();
    
    return config;
}

} // namespace ICARION::config

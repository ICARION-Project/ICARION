// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#ifndef ICARION_CONFIG_OUTPUT_LOADER_H
#define ICARION_CONFIG_OUTPUT_LOADER_H

#include "../types/OutputConfig.h"
#include <json/json.h>

namespace ICARION::config {

/**
 * @brief Loader for output configuration section
 */
class OutputConfigLoader {
public:
    /**
     * @brief Load output config from JSON
     * 
     * @param json JSON object for output section
     * @return OutputConfig Parsed configuration
     */
    static OutputConfig load(const Json::Value& json);
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_OUTPUT_LOADER_H

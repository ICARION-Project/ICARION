// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_CONFIG_PHYSICS_LOADER_H
#define ICARION_CONFIG_PHYSICS_LOADER_H

#include "../types/PhysicsConfig.h"
#include "../conversion/EnumMapper.h"
#include <json/json.h>

namespace ICARION::config {

/**
 * @brief Loader for physics configuration section
 */
class PhysicsConfigLoader {
public:
    /**
     * @brief Load physics config from JSON
     * 
     * @param json JSON object for physics section
     * @return PhysicsConfig Parsed configuration
     */
    static PhysicsConfig load(const Json::Value& json);
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_PHYSICS_LOADER_H

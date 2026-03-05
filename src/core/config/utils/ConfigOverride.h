// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_CONFIG_OVERRIDE_H
#define ICARION_CONFIG_OVERRIDE_H

#include "../types/FullConfig.h"
#include <string>
#include <map>
#include <stdexcept>
#include <sstream>
#include <cstdint>

namespace ICARION::config {

/**
 * @brief Apply command-line config overrides to FullConfig
 * 
 * Supports dot-notation paths like:
 * - simulation.dt_s=1e-10
 * - physics.collision_model=EHSS
 * - simulation.rng_seed=42
 * - output.print_results=true
 */
class ConfigOverride {
public:
    /**
     * @brief Apply all overrides from CLI
     * 
     * @param config FullConfig to modify
     * @param overrides Map of key=value pairs from --set
     */
    static void apply(FullConfig& config, const std::map<std::string, std::string>& overrides);
    
private:
    // Helper: Parse double value
    static double parse_double(const std::string& value, const std::string& key);
    
    // Helper: Parse int value
    static int parse_int(const std::string& value, const std::string& key);
    
    // Helper: Parse bool value
    static bool parse_bool(const std::string& value, const std::string& key);
    
    // Helper: Parse unsigned int value
    static unsigned int parse_uint(const std::string& value, const std::string& key);

    // Helper: Parse uint64 value
    static uint64_t parse_uint64(const std::string& value, const std::string& key);
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_OVERRIDE_H

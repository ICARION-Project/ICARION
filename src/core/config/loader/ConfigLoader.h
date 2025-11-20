// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_LOADER_H
#define ICARION_CONFIG_LOADER_H

#include "../types/FullConfig.h"
#include <filesystem>
#include <json/json.h>
#include <fstream>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief Top-level configuration loader
 * 
 * Orchestrates loading of all config sections from JSON.
 * Entry point for the new config system.
 */
class ConfigLoader {
public:
    /**
     * @brief Load complete configuration from JSON file
     * 
     * @param config_path Path to JSON config file
     * @return FullConfig Complete validated configuration
     * @throws std::runtime_error on parse/validation errors
     */
    static FullConfig load(const std::filesystem::path& config_path);
    
    /**
     * @brief Load configuration from parsed JSON
     * 
     * @param root Pre-parsed JSON root object
     * @param base_path Base path for relative file references (defaults to config file dir)
     * @return FullConfig Complete validated configuration
     */
    static FullConfig load_from_json(const Json::Value& root, 
                                     const std::filesystem::path& base_path = "");
    
private:
    // Helper: Parse JSON file
    static Json::Value parse_json_file(const std::filesystem::path& path);
    
    // Helper: Resolve file paths relative to config file
    static std::string resolve_path(const std::string& path, 
                                    const std::filesystem::path& base_path);
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_LOADER_H

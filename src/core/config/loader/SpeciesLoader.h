// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_CONFIG_SPECIES_LOADER_H
#define ICARION_CONFIG_SPECIES_LOADER_H

#include "../types/SpeciesConfig.h"
#include <filesystem>
#include <json/json.h>

namespace ICARION::config {

/**
 * @brief Species database loader
 * 
 * Loads species from JSON with user-friendly units and converts to SI.
 */
class SpeciesLoader {
public:
    /**
     * @brief Load species database from JSON file
     * 
     * @param filepath Path to species JSON file
     * @return Validated species database with SI-converted units
     * @throws std::runtime_error on parse/validation errors
     */
    static SpeciesDatabase load(const std::filesystem::path& filepath);
    
    /**
     * @brief Load species database from parsed JSON
     * 
     * @param json Parsed JSON root object
     * @param base_path Base path for resolving relative file references (optional)
     * @return Validated species database
     */
    static SpeciesDatabase load_from_json(
        const Json::Value& json,
        const std::filesystem::path& base_path = {}
    );
    
private:
    /**
     * @brief Parse single species entry
     * 
     * @param id Species identifier (dictionary key)
     * @param json Species properties JSON object
     * @return Parsed species (units NOT yet converted to SI)
     */
    static SpeciesProperties parse_species(const std::string& id, const Json::Value& json);
    
    /**
     * @brief Helper: Get optional double field
     */
    static std::optional<double> get_optional_double(const Json::Value& json, const char* key);
    
    /**
     * @brief Helper: Get optional string field
     */
    static std::optional<std::string> get_optional_string(const Json::Value& json, const char* key);
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_SPECIES_LOADER_H

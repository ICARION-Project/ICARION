// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_REACTION_LOADER_H
#define ICARION_CONFIG_REACTION_LOADER_H

#include "../types/ReactionConfig.h"
#include "../types/SpeciesConfig.h"
#include <filesystem>
#include <json/json.h>

namespace ICARION::config {

/**
 * @brief Reaction database loader
 */
class ReactionLoader {
public:
    /**
     * @brief Load reaction database from JSON file
     * 
     * @param filepath Path to reaction JSON file
     * @param species_db Species database for validation (optional)
     * @return Validated reaction database
     * @throws std::runtime_error on parse/validation errors
     */
    static ReactionDatabase load(const std::filesystem::path& filepath,
                                 const SpeciesDatabase* species_db = nullptr);
    
    /**
     * @brief Load reaction database from parsed JSON
     * 
     * @param json Parsed JSON root object
     * @param species_db Species database for validation (optional)
     * @return Validated reaction database
     */
    static ReactionDatabase load_from_json(const Json::Value& json,
                                           const SpeciesDatabase* species_db = nullptr);
    
private:
    /**
     * @brief Parse single reaction entry
     */
    static Reaction parse_reaction(const Json::Value& json);
    
    /**
     * @brief Parse order term
     */
    static ReactionOrderTerm parse_order_term(const Json::Value& json);
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_REACTION_LOADER_H

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "FullConfig.h"
#include "../loader/SpeciesLoader.h"
#include "../loader/ReactionLoader.h"
#include <iostream>

namespace ICARION::config {

// Helper function to find global database by searching up the directory tree
static std::filesystem::path find_global_database(
    const std::filesystem::path& base_path,
    const std::string& relative_db_path) {
    
    std::filesystem::path search_path = base_path;
    
    // Search up to 5 levels up the directory tree
    for (int i = 0; i < 5; ++i) {
        std::filesystem::path candidate = search_path / relative_db_path;
        if (std::filesystem::exists(candidate)) {
            return candidate;
        }
        
        // Go up one level
        std::filesystem::path parent = search_path.parent_path();
        if (parent == search_path) {
            break; // Reached root
        }
        search_path = parent;
    }
    
    return {}; // Not found
}

void FullConfig::load_databases(const std::filesystem::path& base_path) {
    // Global fallback database paths (relative to workspace root)
    const std::string global_species_db = "data/species_database_v1.json";
    const std::string global_reactions_db = "data/reactions_database_v1.json";
    
    // Load species database
    if (!species_database_path.empty()) {
        std::filesystem::path species_path = species_database_path;
        
        // Resolve relative paths
        if (species_path.is_relative()) {
            species_path = base_path / species_path;
        }
        
        std::cout << "[DatabaseLoader] Loading species from: " << species_path << "\n";
        species_db = SpeciesLoader::load(species_path);
    } else {
        // Try global fallback by searching up the directory tree
        std::filesystem::path fallback_species = find_global_database(base_path, global_species_db);
        if (!fallback_species.empty()) {
            std::cout << "[DatabaseLoader] ℹ No species database specified, using global fallback: " 
                      << fallback_species.filename() << "\n";
            species_db = SpeciesLoader::load(fallback_species);
        } else {
            std::cout << "[DatabaseLoader] No species database specified (optional)\n";
        }
    }
    
    // Load reactions (with species validation)
    if (!reaction_database_path.empty()) {
        std::filesystem::path reaction_path = reaction_database_path;
        
        // Resolve relative paths
        if (reaction_path.is_relative()) {
            reaction_path = base_path / reaction_path;
        }
        
        std::cout << "[DatabaseLoader] Loading reactions from: " << reaction_path << "\n";
        
        // Pass species DB for validation if available
        const SpeciesDatabase* species_ptr = (species_db.size() > 0) ? &species_db : nullptr;
        reaction_db = ReactionLoader::load(reaction_path, species_ptr);
    } else {
        // Try global fallback (only if species DB exists)
        if (species_db.size() > 0) {
            std::filesystem::path fallback_reactions = find_global_database(base_path, global_reactions_db);
            if (!fallback_reactions.empty()) {
                std::cout << "[DatabaseLoader] ℹ No reaction database specified, using global fallback: " 
                          << fallback_reactions.filename() << "\n";
                reaction_db = ReactionLoader::load(fallback_reactions, &species_db);
            } else {
                std::cout << "[DatabaseLoader] No reaction database specified (optional)\n";
            }
        } else {
            std::cout << "[DatabaseLoader] No reaction database specified (optional)\n";
        }
    }
}

} // namespace ICARION::config

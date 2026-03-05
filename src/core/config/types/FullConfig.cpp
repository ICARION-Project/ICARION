// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "FullConfig.h"
#include "../loader/SpeciesLoader.h"
#include "../loader/ReactionLoader.h"
#include "../loader/IonLoader.h"
#include "core/types/IonState.h"
#include <iostream>

namespace ICARION::config {

// Helper function to find global database by searching up the directory tree
static std::filesystem::path find_global_database(
    const std::filesystem::path& base_path,
    const std::string& relative_db_path) {
    
    // First try: data/ subdirectory in current working directory (for build/ runs)
    std::filesystem::path cwd = std::filesystem::current_path();
    std::filesystem::path cwd_data = cwd / "data" / relative_db_path;
    if (std::filesystem::exists(cwd_data)) {
        return cwd_data;
    }
    
    // Second try: ../data/ from current working directory (for build/ runs)
    std::filesystem::path parent_data = cwd.parent_path() / "data" / relative_db_path;
    if (std::filesystem::exists(parent_data)) {
        return parent_data;
    }
    
    // Third try: search up from base_path (for config file locations)
    std::filesystem::path search_path = base_path;
    for (int i = 0; i < 5; ++i) {
        std::filesystem::path candidate = search_path / "data" / relative_db_path;
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
    // Global fallback database paths (relative to workspace root or build/)
    const std::string global_species_db = "species_database_v1.json";
    const std::string global_reactions_db = "reactions_database_v1.json";
    
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

std::vector<IonState> FullConfig::generate_ions(std::mt19937& rng) const {
    // Use new IonConfig if specified
    if (ions.is_valid()) {
        IonGenerationResult result = IonLoader::generate_ions(ions, species_db, rng);
        
        if (!result.validation.valid) {
            std::cerr << "Ion generation validation failed:\n";
            result.validation.print();
            throw std::runtime_error("Ion generation failed - check species database and configuration");
        }
        
        if (!result.validation.warnings.empty()) {
            std::cout << "Warnings during ion generation:\n";
            for (const auto& warning : result.validation.warnings) {
                std::cout << "  ⚠  " << warning << "\n";
            }
        }
        
        return result.ions;
    }
    
    // No ion configuration
    throw std::runtime_error("No ion configuration specified ('ions' field required)");
}

} // namespace ICARION::config

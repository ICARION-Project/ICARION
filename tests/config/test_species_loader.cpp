// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
// Test: Load species and reaction databases

#include "core/config/loader/SpeciesLoader.h"
#include "core/config/loader/ReactionLoader.h"
#include <iostream>

int main() {
    using namespace ICARION::config;
    
    try {
        // Load species
        std::cout << "=== Loading Species Database ===\n";
        auto species_db = SpeciesLoader::load("../data/species_database_v1.json");
        
        std::cout << "\nLoaded " << species_db.size() << " species:\n";
        for (const auto& [id, props] : species_db.species) {
            std::cout << "  - " << id;
            if (props.name) {
                std::cout << " (" << *props.name << ")";
            }
            std::cout << ": " << props.mass_amu << " amu, charge=" << props.charge;
            if (props.mobility_cm2Vs) {
                std::cout << ", K₀=" << *props.mobility_cm2Vs << " cm²/Vs";
            }
            if (props.CCS_A2) {
                std::cout << ", CCS=" << *props.CCS_A2 << " Ų";
            }
            std::cout << "\n";
        }
        
        // Load reactions
        std::cout << "\n=== Loading Reactions Database ===\n";
        auto reaction_db = ReactionLoader::load("../data/reactions_database_v1.json", &species_db);
        
        std::cout << "\nLoaded " << reaction_db.size() << " reactions:\n";
        for (const auto& rxn : reaction_db.reactions) {
            std::cout << "  - " << rxn.id << ": " << rxn.reactant << " → " << rxn.product;
            std::cout << " (k=" << rxn.rate_constant << " m³/s)\n";
            
            if (!rxn.order_terms.empty()) {
                for (const auto& term : rxn.order_terms) {
                    std::cout << "      [" << term.species << "]^" << term.exponent;
                    if (term.concentration_m3 > 0) {
                        std::cout << " = " << term.concentration_m3 << " m⁻³";
                    }
                    std::cout << "\n";
                }
            }
        }
        
        std::cout << "\n✓ All databases loaded successfully!\n";
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
        return 1;
    }
}

// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <iostream>
#include <filesystem>
#include "core/config/loader/ConfigLoader.h"

int main() {
    using namespace ICARION::config;
    
    try {
        // Test 0: Load config WITHOUT database paths (should use fallback)
        std::cout << "=== Test 0: IMS with Fallback Databases ===\n";
        std::filesystem::path config_path0 = "/home/chsch95/ICARION/examples/ims_with_fallback_db.json";
        auto config0 = ConfigLoader::load(config_path0);
        
        std::cout << "✓ Config loaded successfully\n";
        std::cout << "  - Domains: " << config0.domains.size() << "\n";
        std::cout << "  - Species: " << config0.species_db.size() << " (loaded from global fallback)\n";
        std::cout << "  - Reactions: " << config0.reaction_db.size() << " (loaded from global fallback)\n";
        
        // Test 1: Load config with species database
        std::cout << "\n=== Test 1: IMS with Species Database ===\n";
        std::filesystem::path config_path1 = "/home/chsch95/ICARION/examples/ims_with_species_db.json";
        auto config1 = ConfigLoader::load(config_path1);
        
        std::cout << "✓ Config loaded successfully\n";
        std::cout << "  - Domains: " << config1.domains.size() << "\n";
        std::cout << "  - Species: " << config1.species_db.size() << "\n";
        std::cout << "  - Reactions: " << config1.reaction_db.size() << "\n";
        
        if (config1.species_db.size() > 0) {
            std::cout << "  - Example species: ";
            size_t count = 0;
            for (const auto& [id, props] : config1.species_db.species) {
                std::cout << id << " ";
                if (++count >= 3) break; // Show first few
            }
            std::cout << "...\n";
        }
        
        // Test 2: Load config with reactions
        std::cout << "\n=== Test 2: Reaction Demo ===\n";
        std::filesystem::path config_path2 = "/home/chsch95/ICARION/examples/reaction_demo.json";
        auto config2 = ConfigLoader::load(config_path2);
        
        std::cout << "✓ Config loaded successfully\n";
        std::cout << "  - Domains: " << config2.domains.size() << "\n";
        std::cout << "  - Species: " << config2.species_db.size() << "\n";
        std::cout << "  - Reactions: " << config2.reaction_db.size() << "\n";
        
        if (config2.reaction_db.size() > 0) {
            std::cout << "  - Example reactions:\n";
            for (const auto& rxn : config2.reaction_db.reactions) {
                std::cout << "      " << rxn.id << ": " << rxn.reactant 
                         << " → " << rxn.product << "\n";
            }
        }
        
        std::cout << "\n✓ All integration tests passed!\n";
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "❌ Error: " << e.what() << "\n";
        return 1;
    }
}

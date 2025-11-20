// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors
//
// Smoke test for new config loader system

#include "core/config/loader/ConfigLoader.h"
#include "core/config/conversion/EnumMapper.h"
#include <iostream>
#include <iomanip>

using namespace ICARION::config;

int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <config.json>" << std::endl;
        return 1;
    }
    
    std::string config_path = argv[1];
    
    std::cout << "=== ICARION Config Loader Smoke Test ===" << std::endl;
    std::cout << "Loading: " << config_path << std::endl << std::endl;
    
    try {
        // Load config
        auto config = ConfigLoader::load(config_path);
        
        // Print simulation config
        std::cout << "--- Simulation ---" << std::endl;
        std::cout << "  Total time: " << config.simulation.total_time_s << " s" << std::endl;
        std::cout << "  Time step: " << config.simulation.dt_s << " s" << std::endl;
        std::cout << "  Total steps: " << config.simulation.total_steps << std::endl;
        std::cout << "  Write interval: " << config.simulation.write_interval << std::endl;
        std::cout << std::endl;
        
        // Print physics config
        std::cout << "--- Physics ---" << std::endl;
        std::cout << "  Collision model: " << EnumMapper::collision_model_to_string(config.physics.collision_model) << std::endl;
        std::cout << "  Reactions: " << (config.physics.enable_reactions ? "enabled" : "disabled") << std::endl;
        std::cout << "  Space charge: " << (config.physics.enable_space_charge ? "enabled" : "disabled") << std::endl;
        std::cout << std::endl;
        
        // Print output config
        std::cout << "--- Output ---" << std::endl;
        std::cout << "  Folder: " << config.output.folder << std::endl;
        std::cout << "  Trajectory file: " << config.output.trajectory_file << std::endl;
        std::cout << "  Print progress: " << (config.output.print_progress ? "yes" : "no") << std::endl;
        std::cout << std::endl;
        
        // Print domains
        std::cout << "--- Domains (" << config.domains.size() << ") ---" << std::endl;
        for (size_t i = 0; i < config.domains.size(); ++i) {
            const auto& domain = config.domains[i];
            std::cout << "  [" << i << "] " << domain.name << std::endl;
            std::cout << "      Instrument: " << EnumMapper::instrument_to_string(domain.instrument) << std::endl;
            std::cout << "      Solver: " << EnumMapper::solver_to_string(domain.solver) << std::endl;
            
            // Geometry
            std::cout << "      Geometry: length=" << domain.geometry.length_m << " m, "
                      << "radius=" << domain.geometry.radius_m << " m" << std::endl;
            
            // Environment
            std::cout << "      Environment: P=" << domain.environment.pressure_Pa << " Pa, "
                      << "T=" << domain.environment.temperature_K << " K, "
                      << "gas=" << domain.environment.gas_species << std::endl;
            
            // Fields
            if (domain.fields.dc.axial_V != 0.0 || domain.fields.dc.EN_Td > 0.0) {
                std::cout << "      DC fields: axial=" << domain.fields.dc.axial_V << " V";
                if (domain.fields.dc.EN_Td > 0.0) {
                    std::cout << ", E/N=" << domain.fields.dc.EN_Td << " Td";
                }
                std::cout << std::endl;
            }
            
            if (domain.fields.rf.voltage_V != 0.0) {
                std::cout << "      RF fields: V=" << domain.fields.rf.voltage_V << " V, "
                          << "f=" << domain.fields.rf.frequency_Hz << " Hz" << std::endl;
            }
            
            if (domain.fields.magnetic.enabled) {
                std::cout << "      Magnetic: B=(" 
                          << domain.fields.magnetic.field_strength_T.x << ", "
                          << domain.fields.magnetic.field_strength_T.y << ", "
                          << domain.fields.magnetic.field_strength_T.z << ") T" << std::endl;
            }
            
            std::cout << std::endl;
        }
        
        // Validate
        std::cout << "--- Validation ---" << std::endl;
        try {
            auto validation = config.validate();
            if (validation.valid) {
                std::cout << "  ✓ Configuration is valid" << std::endl;
            }
            if (!validation.warnings.empty()) {
                std::cout << "  Warnings:" << std::endl;
                for (const auto& warn : validation.warnings) {
                    std::cout << "    ⚠ " << warn << std::endl;
                }
            }
        } catch (const std::exception& e) {
            std::cerr << "  ✗ Validation failed: " << e.what() << std::endl;
            return 1;
        }
        
        std::cout << std::endl << "=== Test completed successfully ===" << std::endl;
        return 0;
        
    } catch (const std::exception& e) {
        std::cerr << "ERROR: " << e.what() << std::endl;
        return 1;
    }
}

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "geometryUtils.h"
#include "core/io/logger.h"
#include <filesystem>
#include <stdexcept>

namespace ICARION::physics {

GeometryMap load_geometry_map(
    const std::unordered_set<std::string>& species_ids,
    const std::string& geometry_file,
    bool strict_mode
) {
    GeometryMap geometry_map;
    
    namespace fs = std::filesystem;
    
    // Check if path is file or directory
    bool is_directory = fs::exists(geometry_file) && fs::is_directory(geometry_file);
    
    for (const auto& species_id : species_ids) {
        try {
            std::string file_path;
            
            if (is_directory) {
                // Try: <geometry_file>/<species_id>.json
                file_path = (fs::path(geometry_file) / (species_id + ".json")).string();
            } else {
                // Single file containing all species
                file_path = geometry_file;
            }
            
            // SSOT: Use central MoleculeLoader
            io::Molecule molecule = io::load_molecule(file_path);
            
            // SSOT: Use central conversion function
            GeometryData geom = convert_molecule_to_geometry(molecule);
            
            if (!geom.first.empty()) {
                io::debug_log(
                    "[GeometryUtils] Loaded geometry for species '" + species_id + 
                    "': " + std::to_string(geom.first.size()) + " atoms"
                );
                geometry_map[species_id] = std::move(geom);
            } else {
                std::string msg = "Empty geometry for species '" + species_id + "'";
                if (strict_mode) {
                    throw std::runtime_error("[GeometryUtils] " + msg);
                }
                io::debug_log("[GeometryUtils] " + msg + " - will fallback to CCS");
                // Empty geometry signals fallback
                geometry_map[species_id] = GeometryData{{}, {}};
            }
            
        } catch (const std::exception& e) {
            std::string msg = "Failed to load geometry for species '" + species_id + 
                            "': " + std::string(e.what());
            if (strict_mode) {
                throw std::runtime_error("[GeometryUtils] " + msg);
            }
            io::debug_log("[GeometryUtils] " + msg + " - will fallback to CCS");
            // Insert empty geometry (signals fallback)
            geometry_map[species_id] = GeometryData{{}, {}};
        }
    }
    
    return geometry_map;
}

} // namespace ICARION::physics

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file geometryUtils.h
 * @brief Utility functions for molecular geometry conversion
 * 
 * SSOT: Single place to convert io::Molecule → GeometryData for collision handlers
 */

#pragma once

#include "EHSSCollisionHandler.h"  // GeometryData, GeometryMap
#include "core/io/moleculeLoader.h"
#include <string>
#include <unordered_set>

namespace ICARION::physics {

/**
 * @brief Convert Molecule to GeometryData for collision handler
 * 
 * @param molecule Molecule loaded from JSON via io::load_molecule()
 * @return GeometryData pair of (atom centers [m], atom radii [m])
 * 
 * **SSOT:** Conversion logic defined ONCE
 * 
 * **Usage:**
 * ```cpp
 * auto molecule = io::load_molecule("H3O+.json");
 * GeometryData geom = convert_molecule_to_geometry(molecule);
 * ```
 */
inline GeometryData convert_molecule_to_geometry(const io::Molecule& molecule,
                                                 bool center_on_com = true) {
    // Conversion factor from Lennard-Jones sigma to hard-sphere radius
    // Physical basis: LJ potential V(r) = 4ε[(σ/r)¹² - (σ/r)⁶]
    // At r=σ, V(σ)=0 (zero-crossing), defining the "contact distance"
    // Hard-sphere model approximates this as sphere with radius r_HS = σ/2
    // Reference: Hirschfelder et al., "Molecular Theory of Gases and Liquids" (1954)
    constexpr double LJ_SIGMA_TO_HS_RADIUS = 0.5;
    
    std::vector<Vec3> centers;
    std::vector<double> radii;
    
    centers.reserve(molecule.atoms.size());
    radii.reserve(molecule.atoms.size());
    
    Vec3 center{0.0, 0.0, 0.0};
    if (center_on_com && !molecule.atoms.empty()) {
        double mass_sum = 0.0;
        for (const auto& atom : molecule.atoms) {
            center = center + atom.pos_m * atom.mass_u;
            mass_sum += atom.mass_u;
        }
        if (mass_sum > 0.0) {
            center = center / mass_sum;
        } else {
            center = Vec3{0.0, 0.0, 0.0};
        }
    }

    for (const auto& atom : molecule.atoms) {
        // NOTE: MoleculeLoader already converts JSON positions (Å) to SI units (m)
        // in atom.pos_m, so no further conversion needed here.
        centers.push_back(atom.pos_m - center);
        
        radii.push_back(LJ_SIGMA_TO_HS_RADIUS * atom.LJ_sigma_m);
    }
    
    return GeometryData{std::move(centers), std::move(radii)};
}

/**
 * @brief Load geometry for multiple species from molecule file(s)
 * 
 * @param species_ids Set of species to load (e.g., {"H3O+", "NH4+"})
 * @param geometry_file Path to single JSON file or directory containing <species>.json
 * @param strict_mode If true, throw on load failure; if false (default), return empty geometry
 * @return GeometryMap with loaded geometries
 * 
 * **Behavior:**
 * - If geometry_file is a single file: Loads all species from that file
 * - If geometry_file is directory: Loads each species from <dir>/<species_id>.json
 * - Missing species (strict_mode=false): Returns empty GeometryData (handler will fallback to CCS)
 * - Missing species (strict_mode=true): Throws std::runtime_error
 * 
 * **SSOT:** Uses io::load_molecule() + convert_molecule_to_geometry()
 * 
 * **Usage:**
 * ```cpp
 * // Silent fallback (default)
 * auto geometry_map = load_geometry_map({"H3O+", "N2"}, "molecules/");
 * 
 * // Strict mode (throw on error)
 * auto geometry_map = load_geometry_map({"H3O+", "N2"}, "molecules/", true);
 * ```
 */
GeometryMap load_geometry_map(
    const std::unordered_set<std::string>& species_ids,
    const std::string& geometry_file,
    bool strict_mode = false
);

} // namespace ICARION::physics

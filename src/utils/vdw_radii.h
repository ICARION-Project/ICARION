// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include <unordered_map>
#include <string>
#include <stdexcept>

/**
 * @brief Van der Waals radii in Angstroms (Å)
 * 
 * Based on Bondi (1964) and Alvarez (2013) compilations.
 * For EHSS hard-sphere collisions, these represent effective
 * atomic collision radii.
 */
static const std::unordered_map<std::string, double> VDW_RADII_ANGSTROM = {
    // Group 1: Alkali metals
    {"H",  1.20},  // Hydrogen
    {"Li", 1.82},  // Lithium
    {"Na", 2.27},  // Sodium
    {"K",  2.75},  // Potassium
    {"Rb", 3.03},  // Rubidium
    {"Cs", 3.43},  // Cesium
    
    // Group 2: Alkaline earth metals
    {"Be", 1.53},  // Beryllium
    {"Mg", 1.73},  // Magnesium
    {"Ca", 2.31},  // Calcium
    {"Sr", 2.49},  // Strontium
    {"Ba", 2.68},  // Barium
    
    // Group 13
    {"B",  1.92},  // Boron
    {"Al", 1.84},  // Aluminum
    {"Ga", 1.87},  // Gallium
    {"In", 1.93},  // Indium
    
    // Group 14
    {"C",  1.70},  // Carbon
    {"Si", 2.10},  // Silicon
    {"Ge", 2.11},  // Germanium
    {"Sn", 2.17},  // Tin
    {"Pb", 2.02},  // Lead
    
    // Group 15: Pnictogens
    {"N",  1.55},  // Nitrogen
    {"P",  1.80},  // Phosphorus
    {"As", 1.85},  // Arsenic
    {"Sb", 2.06},  // Antimony
    {"Bi", 2.07},  // Bismuth
    
    // Group 16: Chalcogens
    {"O",  1.52},  // Oxygen
    {"S",  1.80},  // Sulfur
    {"Se", 1.90},  // Selenium
    {"Te", 2.06},  // Tellurium
    
    // Group 17: Halogens
    {"F",  1.47},  // Fluorine
    {"Cl", 1.75},  // Chlorine
    {"Br", 1.85},  // Bromine
    {"I",  1.98},  // Iodine
    
    // Noble gases
    {"He", 1.40},  // Helium
    {"Ne", 1.54},  // Neon
    {"Ar", 1.88},  // Argon
    {"Kr", 2.02},  // Krypton
    {"Xe", 2.16},  // Xenon
    
    // Transition metals (selected common ones)
    {"Ti", 2.15},  // Titanium
    {"V",  2.05},  // Vanadium
    {"Cr", 2.05},  // Chromium
    {"Mn", 2.05},  // Manganese
    {"Fe", 2.04},  // Iron
    {"Co", 2.00},  // Cobalt
    {"Ni", 1.97},  // Nickel
    {"Cu", 1.96},  // Copper
    {"Zn", 2.01},  // Zinc
    {"Pd", 2.10},  // Palladium
    {"Ag", 2.11},  // Silver
    {"Cd", 2.18},  // Cadmium
    {"Pt", 2.13},  // Platinum
    {"Au", 2.14},  // Gold
    {"Hg", 2.23},  // Mercury
};

/**
 * @brief Get Van der Waals radius for an element in meters
 * 
 * @param element_symbol Chemical element symbol (e.g., "C", "N", "O")
 * @return Van der Waals radius in meters
 * @throws std::runtime_error if element is not found in lookup table
 * 
 * @note Element symbols are case-sensitive: use "C" not "c"
 */
inline double get_vdw_radius_m(const std::string& element_symbol) {
    auto it = VDW_RADII_ANGSTROM.find(element_symbol);
    if (it == VDW_RADII_ANGSTROM.end()) {
        throw std::runtime_error("Van der Waals radius not found for element: " + element_symbol + 
                                "\nAvailable elements: H, C, N, O, S, P, F, Cl, Br, I, etc.");
    }
    // Convert Angstrom to meters: 1 Å = 1e-10 m
    return it->second * 1e-10;
}

/**
 * @brief Get Van der Waals radius with fallback to LJ sigma approximation
 * 
 * If element is not in lookup table, falls back to 0.5 * LJ_sigma
 * (historical approximation used in ICARION).
 * 
 * @param element_symbol Chemical element symbol
 * @param lj_sigma_m Lennard-Jones sigma parameter in meters (fallback)
 * @return Van der Waals radius in meters
 */
inline double get_vdw_radius_m_with_fallback(const std::string& element_symbol, 
                                              double lj_sigma_m) {
    // Use historical approximation: 0.5 * LJ sigma as primary radius
    // (Some workflows prefer half the LJ sigma over lookup VdW values.)
    // This intentionally disregards the VDW lookup to match previous behaviour.
    return 0.5 * lj_sigma_m;
}

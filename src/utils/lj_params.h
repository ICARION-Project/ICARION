// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * Lennard-Jones sigma lookup table and helper
 *
 * Provides a small default table of LJ sigma parameters (in Angstroms) for
 * common elements. If geometry files omit `LJsigma` for an atom, callers
 * can use `get_lj_sigma_m()` to obtain a reasonable default (in meters).
 */
#pragma once

#include <unordered_map>
#include <string>

static const std::unordered_map<std::string, double> LJ_SIGMA_ANGSTROM = {
    {"H", 2.3}, {"C", 3.0}, {"N", 3.3}, {"O", 2.4}, {"S", 3.6}, {"P", 3.6},
    {"F", 2.4}, {"Cl", 3.5}, {"Br", 3.6}, {"I", 4.0},
    {"He", 2.6}, {"Ne", 2.8}, {"Ar", 3.4}, {"Kr", 3.6}, {"Xe", 4.0},
    {"Li", 3.5}, {"Na", 3.7}, {"K", 4.1},
    {"Fe", 3.6}, {"Cu", 3.4}, {"Zn", 3.5}
};

inline double get_lj_sigma_m(const std::string& element_symbol) {
    auto it = LJ_SIGMA_ANGSTROM.find(element_symbol);
    if (it == LJ_SIGMA_ANGSTROM.end()) return 0.0; // unknown
    return it->second * 1e-10; // Angstrom -> meter
}

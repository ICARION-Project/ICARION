// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        reactionUtils.h
 *   @brief       Load chemical species and reaction definitions
 *
 *   @details
 * Provides routines to parse JSON files containing chemical species
 * and reaction information. Species entries include mass, charge, mobility,
 * and collision cross-section (CCS), while reactions define reactants,
 * products, rate constants, and reaction orders. 
 *
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "core/param/paramUtils.h"

// -----------------------------
// Structs
// -----------------------------

struct ReactionOrderTerm {
    std::string species;  // name of neutral partner
    int         exponent;
};

struct ReactionEntry {
    std::string                    reactant;       // Ion before reaction
    std::string                    product;        // Ion after reaction
    double                         rate_constant;  // SI: m^(3n-1)/s (depending on reaction order)
    double                         rate;
    std::vector<ReactionOrderTerm> order;
    double                         neutral_concentration;
};

struct Species {
    std::string name;
    double      mass_kg;
    double      mobility;
    double      charge;
    double      CCS;
};

// -----------------------------
// Loader functions
// -----------------------------

/**
 * @brief Load chemical species database from JSON file
 * @deprecated Use config::SpeciesLoader instead (Phase 3 refactor, 2025-11-22)
 * 
 * This legacy function will be removed in Phase 3D after database unification.
 * Use config::SpeciesLoader::load_species() for new code.
 * 
 * @param gParams Global simulation parameters (contains speciesDB_json path)
 * @param temperature_K Gas temperature [K] for mobility calculations
 * @param neutral_mass_kg Mass of background gas molecules [kg]
 * 
 * @return Map of species name to Species properties (mass, charge, mobility, CCS)
 * 
 * Parses JSON file specified in gParams.speciesDB_json containing species definitions.
 * Each species entry includes:
 * - name: Species identifier (e.g., "H3O+", "NO2+")
 * - mass_kg: Ion mass in kilograms
 * - charge: Ion charge in elementary charges
 * - mobility: Ion mobility [m²/(V·s)] (optional, can be calculated from CCS)
 * - CCS: Collision cross-section [Å²] (optional, can be calculated from mobility)
 * 
 * If mobility is not provided, it's calculated from CCS using Mason-Schamp equation.
 * If CCS is not provided, it's calculated from mobility.
 * 
 * @see src/core/config/loaders/SpeciesLoader.h
 * @see src/core/config/types/SpeciesConfig.h
 */
std::unordered_map<std::string, Species> load_speciesDB(
    const GlobalParams& gParams,
    double temperature_K = 300.0,
    double neutral_mass_kg = MOLAR_MASS_HE_KG
);

/**
 * @brief Load ion-molecule reactions from JSON file
 * @deprecated Use config::ReactionLoader instead (Phase 3 refactor, 2025-11-22)
 * 
 * This legacy function will be removed in Phase 3D after database unification.
 * Use config::ReactionLoader::load_reactions() for new code.
 * 
 * @param gParams Global simulation parameters (contains reactions_json path)
 * 
 * @return Vector of reaction entries with rate constants and products
 * 
 * @see src/core/config/loaders/ReactionLoader.h
 * @see src/core/config/types/ReactionConfig.h
 * 
 * Parses JSON file specified in gParams.reactions_json containing reaction definitions.
 * Each reaction entry includes:
 * - reactant: Reactant ion species name
 * - product: Product ion species name
 * - rate_constant: Reaction rate constant [m^(3n-1)/s] where n is total reaction order
 * - order: Vector of {species, exponent} for reaction order terms
 * 
 * Reaction rates are calculated at each timestep based on local neutral concentrations
 * and ion densities. Supports multi-step reactions with different neutrals.
 */
std::vector<ReactionEntry>               load_reactions(const GlobalParams& gParams);
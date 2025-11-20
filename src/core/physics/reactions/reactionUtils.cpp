/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        reactionUtils.cpp
 *   @brief       Load chemical species and reaction definitions
 *
 *   @details
 * Provides routines to parse JSON files containing chemical species
 * and reaction information. Species entries include mass, charge, mobility,
 * and collision cross-section (CCS), while reactions define reactants,
 * products, rate constants, and reaction orders. 
 *
 * These are used for stochastic reaction handling during ion trajectory simulations.
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include "reactionUtils.h"
#include <cmath>
#include <fstream>
#include <iostream>
#include <nlohmann/json.hpp>
#include "core/param/paramUtils.h"

using json = nlohmann::json;

/**
 * @brief Load species database from JSON and compute derived parameters.
 *
 * @param[in] gParams Global simulation parameters (contains file paths, defaults).
 * @param[in] temperature_K Temperature used for CCS calculation [K].
 * @param[in] neutral_mass_kg Mass of neutral particle used for CCS calculation [kg].
 * @return std::unordered_map<std::string, Species> Map of species name to Species object.
 *
 * @throws std::runtime_error If the JSON file cannot be opened or parsed.
 *
 * @details For each species:
 * - Converts mass from amu to kg.
 * - Computes collision cross-section (CCS) using the given temperature and neutral mass.
 * - Stores charge in Coulombs.
 * - Retains user-specified mobility.
 *
 * @note CCS formula is derived from standard ion mobility theory.
 */

std::unordered_map<std::string, Species> load_speciesDB(
    const GlobalParams& gParams,
    double temperature_K,
    double neutral_mass_kg
) {
    // If reaction_file is empty or missing, return empty DB
    if (gParams.reaction_file.empty()) {
        return std::unordered_map<std::string, Species>();
    }
    
    std::ifstream file(gParams.reaction_file);
    if (!file) {
        // Check if reactions are actually enabled
        if (gParams.enable_reactions) {
            throw std::runtime_error("Could not open species JSON file: " + gParams.reaction_file);
        }
        // If reactions disabled, just return empty
        return std::unordered_map<std::string, Species>();
    }
    json j;
    file >> j;

    std::unordered_map<std::string, Species> db;
    db.reserve(j["species"].size());
    
    for (auto& sp : j["species"]) {
        Species s;
        s.name          = sp["name"].get<std::string>();
        double mass_amu = sp["mass_amu"].get<double>();
        s.mass_kg       = mass_amu * AMU_TO_KG;
        s.mobility      = sp["mobility"].get<double>();
        s.charge        = sp["charge"].get<int>() * ELECTRON_CHARGE;
        s.CCS           = 3.0 / 16.0 / LOSCHMIDT_CONSTANT * s.charge *
                std::sqrt(2.0 * M_PI /
                          (BOLTZMANN_CONSTANT * temperature_K *
                           ((s.mass_kg * neutral_mass_kg) /
                            (s.mass_kg + neutral_mass_kg)))) *
                1.0 / (s.mobility * 1e-4);
        db.emplace(std::move(s.name), std::move(s));
    }
    return db;
}

/**
 * @brief Load reaction entries from JSON file.
 *
 * @param[in] gParams Global simulation parameters (provides reaction file path).
 * @return std::vector<ReactionEntry> Vector of reactions.
 *
 * @throws std::runtime_error If the reactions JSON file cannot be opened or parsed.
 *
 * @details Each reaction entry includes:
 * - Reactant and product species names.
 * - Rate constant (scaled as necessary in simulation).
 * - Neutral concentration (optional).
 * - Reaction order terms specifying species and exponents.
 */
std::vector<ReactionEntry> load_reactions(const GlobalParams& gParams) {
    // If no reaction_file was provided, return an empty list. Do not gate
    // parsing on `enable_reactions` — callers/tests expect the loader to
    // populate the reaction list when a file is present even if reactions
    // are later disabled at runtime.
    if (gParams.reaction_file.empty()) {
        return std::vector<ReactionEntry>();
    }
    
    std::ifstream file(gParams.reaction_file);
    if (!file) {
        if (gParams.enable_reactions) {
            throw std::runtime_error("Could not open reactions JSON file: " + gParams.reaction_file);
        }
        return std::vector<ReactionEntry>();
    }
    json j;
    file >> j;

    std::vector<ReactionEntry> reactions;
    reactions.reserve(j["reactions"].size());

    for (auto& r : j["reactions"]) {
        ReactionEntry rxn;
        rxn.reactant = r["reactant"].get<std::string>();
        rxn.product  = r["product"].get<std::string>();
    // Rate constant: keep value as provided in the JSON (tests expect SI value)
    rxn.rate_constant = r["rate_constant"].get<double>();
        rxn.neutral_concentration = r["neutral_concentration"].get<double>(); // given in volume fraction

        rxn.order.reserve(r["order"].size());
        for (auto& term : r["order"]) {
            ReactionOrderTerm t;
            t.species  = term["species"].get<std::string>();
            t.exponent = term["exponent"].get<int>();
            rxn.order.push_back(t);
        }

        reactions.push_back(std::move(rxn));
    }

    return reactions;
}
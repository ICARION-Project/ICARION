// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_REACTION_CONFIG_H
#define ICARION_CONFIG_REACTION_CONFIG_H

#include "SpeciesConfig.h"
#include "../validation/ValidationResult.h"
#include <string>
#include <vector>
#include <unordered_map>

namespace ICARION::config {

/**
 * @brief Concentration dependence term
 * 
 * Specifies how reaction rate depends on neutral concentration.
 * Rate = k * [species]^exponent
 */
struct ReactionOrderTerm {
    std::string species;            ///< Neutral species ID
    int exponent;                   ///< Concentration exponent (0, 1, or 2)
    double concentration_m3 = -1.0; ///< Number density [m⁻³] (-1.0 = use buffer gas)
    
    ValidationResult validate() const {
        ValidationResult result;
        
        if (species.empty()) {
            result.add_error("Order term: species cannot be empty");
        }
        
        if (exponent < 0 || exponent > 2) {
            result.add_error("Order term: exponent must be 0, 1, or 2 (got " + 
                           std::to_string(exponent) + ")");
        }
        
        // concentration_m3 validation:
        // -1.0 = "not specified" (use buffer gas density)
        // >= 0.0 = explicit concentration
        // Other negative values = error
        if (concentration_m3 < -1.0 || (concentration_m3 > -1.0 && concentration_m3 < 0.0)) {
            result.add_error("Order term: concentration must be >= 0 or -1.0 (buffer gas fallback), got " +
                           std::to_string(concentration_m3));
        }
        
        return result;
    }
};

/**
 * @brief Simple single-reactant reaction
 * 
 * Models: A⁺ + n·X → B⁺ + products
 * where A⁺ is reactant ion, X is neutral, B⁺ is product ion.
 * 
 * No temperature dependence, no multi-step reactions (reserved for v2.0).
 */
struct Reaction {
    // === Required fields ===
    std::string id;                         ///< Unique reaction identifier
    std::string reactant;                   ///< Reactant ion species ID
    std::string product;                    ///< Product ion species ID
    double rate_constant_m3s;               ///< Rate constant [m³/s] for 2-body
    
    // === Optional fields ===
    std::vector<ReactionOrderTerm> order_terms; ///< Concentration dependence
    
    /**
     * @brief Calculate effective rate constant
     * 
     * @param concentrations Map of species ID → concentration [m⁻³]
     * @return k_eff [s⁻¹] including concentration dependence
     * 
     * For 2-body: k_eff = k * [X]
     * For 3-body: k_eff = k * [X] * [M]
     */
    double effective_rate_s(const std::unordered_map<std::string, double>& concentrations) const {
        double k_eff = rate_constant_m3s;
        
        for (const auto& term : order_terms) {
            auto it = concentrations.find(term.species);
            double conc = (it != concentrations.end()) ? it->second : term.concentration_m3;
            
            // Apply concentration^exponent
            for (int i = 0; i < term.exponent; ++i) {
                k_eff *= conc;
            }
        }
        
        return k_eff;
    }
    
    /**
     * @brief Validate reaction
     * 
     * @param species_db Species database for ID validation
     */
    ValidationResult validate(const SpeciesDatabase& species_db) const {
        ValidationResult result;
        
        // Check required fields
        if (id.empty()) {
            result.add_error("Reaction: ID cannot be empty");
        }
        
        if (reactant.empty()) {
            result.add_error("Reaction '" + id + "': reactant cannot be empty");
        } else if (!species_db.has(reactant)) {
            result.add_error("Reaction '" + id + "': reactant '" + reactant + 
                           "' not found in species database");
        } else {
            // Check reactant is an ion
            const auto& reactant_species = species_db.get(reactant);
            if (reactant_species.charge == 0) {
                result.add_warning("Reaction '" + id + "': reactant '" + reactant + 
                                 "' is neutral (expected ion)");
            }
        }
        
        if (product.empty()) {
            result.add_error("Reaction '" + id + "': product cannot be empty");
        } else if (!species_db.has(product)) {
            result.add_error("Reaction '" + id + "': product '" + product + 
                           "' not found in species database");
        } else {
            // Check product is an ion
            const auto& product_species = species_db.get(product);
            if (product_species.charge == 0) {
                result.add_warning("Reaction '" + id + "': product '" + product + 
                                 "' is neutral (expected ion)");
            }
        }
        
        if (rate_constant_m3s <= 0.0) {
            result.add_error("Reaction '" + id + "': rate constant must be positive");
        }
        
        // Validate order terms
        for (const auto& term : order_terms) {
            result.merge(term.validate());
            
            // Check species exists
            if (!species_db.has(term.species)) {
                result.add_error("Reaction '" + id + "': order term species '" + 
                               term.species + "' not found in database");
            }
        }
        
        // Check charge conservation (optional, informational)
        if (species_db.has(reactant) && species_db.has(product)) {
            int reactant_charge = species_db.get(reactant).charge;
            int product_charge = species_db.get(product).charge;
            
            if (reactant_charge != product_charge) {
                result.add_warning("Reaction '" + id + "': Charge changes from " + 
                                 std::to_string(reactant_charge) + " to " + 
                                 std::to_string(product_charge) + 
                                 " (ensure neutral byproducts are accounted for)");
            }
        }
        
        return result;
    }
};

/**
 * @brief Complete reaction database
 */
struct ReactionDatabase {
    std::vector<Reaction> reactions;
    
    /**
     * @brief Get reactions for a specific reactant species
     * 
     * @param species_id Reactant species ID
     * @return Vector of pointers to matching reactions
     */
    std::vector<const Reaction*> get_reactions_for(const std::string& species_id) const {
        std::vector<const Reaction*> result;
        
        for (const auto& rxn : reactions) {
            if (rxn.reactant == species_id) {
                result.push_back(&rxn);
            }
        }
        
        return result;
    }
    
    /**
     * @brief Get number of reactions
     */
    size_t size() const {
        return reactions.size();
    }
    
    /**
     * @brief Validate database
     * 
     * Checks:
     * - All reactant/product IDs exist in species database
     * - No duplicate reaction IDs
     * - Rate constants are positive
     */
    ValidationResult validate(const SpeciesDatabase& species_db) const {
        ValidationResult result;
        
        if (reactions.empty()) {
            result.add_warning("Reaction database is empty");
            return result;
        }
        
        // Check for duplicate IDs
        std::unordered_map<std::string, int> id_counts;
        for (const auto& rxn : reactions) {
            id_counts[rxn.id]++;
        }
        
        for (const auto& [id, count] : id_counts) {
            if (count > 1) {
                result.add_error("Duplicate reaction ID: '" + id + "' (" + 
                               std::to_string(count) + " occurrences)");
            }
        }
        
        // Validate each reaction
        for (const auto& rxn : reactions) {
            result.merge(rxn.validate(species_db));
        }
        
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_REACTION_CONFIG_H

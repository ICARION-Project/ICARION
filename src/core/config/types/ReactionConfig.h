// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_CONFIG_REACTION_CONFIG_H
#define ICARION_CONFIG_REACTION_CONFIG_H

#include "SpeciesConfig.h"
#include "../validation/ValidationResult.h"
#include <algorithm>
#include <string>
#include <vector>
#include <unordered_map>

namespace ICARION::config {

inline bool is_reaction_order_placeholder(const std::string& species) {
    return species == "neutral" || species == "M";
}

/**
 * @brief Temperature dependence model for rate constant
 * 
 * Determines how k varies with temperature T [K].
 */
enum class RateModel {
    Constant,          ///< k(T) = k₀ (no temperature dependence)
    Arrhenius,         ///< k(T) = A × exp(-Eₐ / (kB·T))
    ModifiedArrhenius  ///< k(T) = A × (T/T₀)ⁿ × exp(-Eₐ / (kB·T))
};

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
 * @brief Single-reactant reaction with optional temperature dependence
 * 
 * Models: A⁺ + n·X → B⁺ + products
 * where A⁺ is reactant ion, X is a neutral concentration term, and B⁺ is
 * product ion. Neutral products are not explicit state variables.
 * 
 * Supports temperature-dependent rate constants (Arrhenius, modified Arrhenius).
 * Multi-step reactions reserved for v2.0.
 * 
 * ⚠️ DIMENSIONAL CONSISTENCY:
 * rate_constant has order-dependent dimensions:
 * - 0 order terms (spontaneous):     k [s⁻¹]
 * - 1 order term with exponent=1:    k [m³/s]
 * - 1 order term with exponent=2:    k [m⁶/s]
 * - 2 order terms with exponent=1:   k [m⁶/s]
 * 
 * Example: For A⁺ + 2X → B⁺, use exponent=2 and k in [m⁶/s]!
 */
struct Reaction {
    // === Required fields ===
    std::string id;                         ///< Unique reaction identifier
    std::string reactant;                   ///< Reactant ion species ID
    std::string product;                    ///< Product ion species ID
    double rate_constant;                   ///< Base rate constant (k₀ or A), dimensions depend on reaction order
    
    // === Temperature dependence (optional) ===
    RateModel rate_model = RateModel::Constant;  ///< Temperature dependence model
    double activation_energy_eV = 0.0;           ///< Activation energy Eₐ [eV] (Arrhenius)
    double temperature_exponent = 0.0;           ///< Temperature exponent n (modified Arrhenius)
    double reference_temperature_K = 300.0;      ///< Reference temperature T₀ [K] (modified Arrhenius)
    
    // === Concentration dependence (optional) ===
    std::vector<ReactionOrderTerm> order_terms; ///< Concentration dependence

    // === Thermochemistry / equilibrium metadata (optional) ===
    // All thermo values are stored internally in SI units.
    // The loader accepts SI keys (delta_r_H_J_mol, delta_r_S_J_molK) as primary
    // and legacy keys (delta_r_H_kcalmol, delta_r_S_cal_molK) with immediate conversion.
    // v1.1 equilibrium is a population-balance model: order_terms are also the
    // thermodynamic neutral partners, except third-body placeholder M.
    bool equilibrium = false;                     ///< Marks forward reaction as equilibrium-controlled
    bool has_thermo = false;                      ///< True if delta_r_H/S were provided in input
    double delta_r_H_J_mol = 0.0;                ///< Reaction enthalpy ΔrH° [J/mol] (SI, SSOT)
    double delta_r_S_J_molK = 0.0;               ///< Reaction entropy ΔrS° [J/(mol·K)] (SI, SSOT)

    // === Reverse-link metadata (internal/runtime) ===
    bool reverse_dynamic_from_equilibrium = false; ///< Compute reverse rate from linked forward + K_eq(T)
    std::string linked_forward_id;                 ///< Forward reaction ID used by dynamic reverse evaluation
    
    /**
     * @brief Compute temperature-dependent rate constant
     * 
     * @param temperature_K Gas temperature [K]
     * @return k(T) with correct dimensions [m³/s for 2nd-order, m⁶/s for 3rd-order]
     * 
     * **Models:**
     * - Constant: k(T) = k₀
     * - Arrhenius: k(T) = A × exp(-Eₐ / (kB·T))
     * - ModifiedArrhenius: k(T) = A × (T/T₀)ⁿ × exp(-Eₐ / (kB·T))
     */
    double compute_rate_constant(double temperature_K) const;
    
    /**
     * @brief Calculate effective rate (includes temperature + concentration)
     * 
     * @param temperature_K Gas temperature [K]
     * @param concentrations Map of species ID → concentration [m⁻³]
     * @return k_eff [s⁻¹] including both T-dependence and concentration
     * 
     * **Formula:** k_eff = k(T) × ∏ᵢ [Xᵢ]^nᵢ
     * 
     * For 2-body: k_eff = k(T) * [X]
     * For 3-body: k_eff = k(T) * [X] * [M]
     */
    double effective_rate_s(
        double temperature_K,
        const std::unordered_map<std::string, double>& concentrations
    ) const {
        // Step 1: Compute k(T)
        double k_T = compute_rate_constant(temperature_K);
        
        // Step 2: Apply concentration terms
        double k_eff = k_T;
        for (const auto& term : order_terms) {
            auto it = concentrations.find(term.species);
            // Fallback: use explicit concentration_m3 only when >= 0.
            // concentration_m3 == -1.0 is a sentinel meaning "use buffer gas density";
            // if the caller did not insert the species into `concentrations`, fall back to 0.0
            // rather than using the raw sentinel and producing a negative rate.
            double conc;
            if (it != concentrations.end()) {
                conc = it->second;
            } else if (term.concentration_m3 >= 0.0) {
                conc = term.concentration_m3;
            } else {
                conc = 0.0;  // sentinel -1.0: species absent from map → treat as zero
            }
            conc = std::max(0.0, conc);  // guard: densities must be non-negative
            
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
        
        if (rate_constant <= 0.0) {
            result.add_error("Reaction '" + id + "': rate constant must be positive");
        }
        
        // Validate order terms
        for (const auto& term : order_terms) {
            result.merge(term.validate());
            
            // Check species exists; "neutral" and "M" are runtime placeholders.
            if (!is_reaction_order_placeholder(term.species) && !species_db.has(term.species)) {
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

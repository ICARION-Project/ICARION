// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <string>
#include <vector>
#include <unordered_map>

#include "core/io/speciesLoader.h"

namespace ICARION {
namespace io {

/**
 * @brief Reaction type classification
 */
enum class ReactionType {
    ChargeTransfer,    ///< Ion + Neutral → Neutral + Ion
    ProtonTransfer,    ///< Ion + Neutral → Product Ion + Product Neutral
    Association,       ///< Ion + Neutral → Cluster Ion
    Dissociation,      ///< Cluster Ion → Ion + Neutral
    Switching,         ///< Ion + Neutral → Different Ion + Different Neutral
    Unknown
};

/**
 * @brief Reaction order term
 * 
 * Specifies concentration dependence for a particular species.
 * Rate equation: rate = k * [A]^n * [B]^m * ...
 */
struct ReactionOrderTerm {
    std::string species_id;    ///< Species identifier
    int exponent;              ///< Concentration exponent (typically 0, 1, or 2)
};

/**
 * @brief Complete reaction definition
 * 
 * Defines a single chemical reaction with rate constants and products.
 */
struct Reaction {
    std::string id;                              ///< Unique reaction identifier
    ReactionType type;                           ///< Classification of reaction type
    std::vector<std::string> reactants;          ///< Reactant species IDs
    std::vector<std::string> products;           ///< Product species IDs
    double rate_constant_SI;                     ///< Rate constant [SI units depend on order]
    double activation_energy_eV;                 ///< Activation energy [eV] (default: 0)
    double branching_ratio;                      ///< Branching ratio (0-1, default: 1)
    std::vector<ReactionOrderTerm> order_terms;  ///< Concentration dependence
    
    /**
     * @brief Calculate temperature-dependent rate constant
     * 
     * @param temperature_K Gas temperature [K]
     * @return Effective rate constant at given temperature [SI]
     * 
     * Uses Arrhenius equation: k(T) = k0 * exp(-Ea / (kB * T))
     */
    double effective_rate_constant(double temperature_K) const;
    
    /**
     * @brief Validate reaction definition
     * 
     * @param species_db Species database for ID validation
     * @throws std::runtime_error if validation fails
     * 
     * Checks:
     * - All reactant and product IDs exist in species database
     * - Rate constant is positive
     * - Branching ratio is between 0 and 1
     * - Order terms reference valid species
     */
    void validate(const SpeciesDatabase& species_db) const;
    
    /**
     * @brief Get total reaction order
     * 
     * @return Sum of all order term exponents
     */
    int total_order() const;
};

/**
 * @brief Reaction database with fast lookup
 * 
 * Provides efficient access to reactions by ID or reactant species.
 */
class ReactionDatabase {
public:
    ReactionDatabase() = default;
    
    /**
     * @brief Add reaction to database
     * 
     * @param reaction Reaction to add
     * @throws std::runtime_error if reaction ID already exists
     */
    void add(const Reaction& reaction);
    
    /**
     * @brief Get reaction by ID
     * 
     * @param id Reaction identifier
     * @return Reference to reaction
     * @throws std::runtime_error if reaction not found
     */
    const Reaction& get(const std::string& id) const;
    
    /**
     * @brief Check if reaction exists
     * 
     * @param id Reaction identifier
     * @return true if reaction exists
     */
    bool has(const std::string& id) const;
    
    /**
     * @brief Get all reactions involving a species
     * 
     * @param species_id Species identifier
     * @return Vector of pointers to reactions where species is a reactant
     */
    std::vector<const Reaction*> get_reactions_for_species(const std::string& species_id) const;
    
    /**
     * @brief Get all reactions
     * 
     * @return Vector of all reactions in insertion order
     */
    const std::vector<Reaction>& all() const { return reactions_list_; }
    
    /**
     * @brief Get number of reactions
     * 
     * @return Number of reactions in database
     */
    size_t size() const { return reactions_list_.size(); }
    
    /**
     * @brief Validate all reactions against species database
     * 
     * @param species_db Species database for validation
     * @throws std::runtime_error if any reaction is invalid
     */
    void validate_all(const SpeciesDatabase& species_db) const;

private:
    std::vector<Reaction> reactions_list_;                        ///< Ordered list of reactions
    std::unordered_map<std::string, size_t> id_to_index_;        ///< Fast ID lookup
    std::unordered_map<std::string, std::vector<size_t>> species_to_reactions_;  ///< Species → reactions lookup
    
    /**
     * @brief Update species-to-reactions index
     * 
     * @param reaction_idx Index of reaction in reactions_list_
     */
    void update_species_index(size_t reaction_idx);
};

/**
 * @brief Load reactions from JSON file
 * 
 * @param filepath Path to reactions JSON file
 * @return Reaction database
 * @throws std::runtime_error on file access errors or parsing failures
 * 
 * Expected JSON format:
 * ```json
 * {
 *   "reactions": [
 *     {
 *       "id": "rxn_001",
 *       "type": "charge_transfer",
 *       "reactants": ["H3O+", "Pentanal"],
 *       "products": ["PentanalH+", "H2O"],
 *       "rate_constant_SI": 1.2e-15,
 *       "activation_energy_eV": 0.05,
 *       "branching_ratio": 1.0,
 *       "order": [
 *         {"species": "Pentanal", "exponent": 1}
 *       ]
 *     }
 *   ]
 * }
 * ```
 * 
 * Required fields: id, reactants, products, rate_constant_SI
 * Optional fields: type, activation_energy_eV, branching_ratio, order
 * 
 * Note: Does not validate against species database at load time.
 * Call validate_all() after loading species database.
 */
ReactionDatabase load_reactions(const std::string& filepath);

/**
 * @brief Load and validate reactions against species database
 * 
 * @param filepath Path to reactions JSON file
 * @param species_db Species database for validation
 * @return Validated reaction database
 * @throws std::runtime_error on file errors, parsing failures, or validation errors
 * 
 * Convenience function that loads reactions and immediately validates
 * all species IDs against the provided species database.
 */
ReactionDatabase load_reactions(const std::string& filepath, 
                                const SpeciesDatabase& species_db);

/**
 * @brief Parse reaction type from string
 * 
 * @param type_str Reaction type string
 * @return ReactionType enum value
 */
ReactionType parse_reaction_type(const std::string& type_str);

/**
 * @brief Convert reaction type to string
 * 
 * @param type Reaction type enum
 * @return Human-readable type string
 */
std::string reaction_type_to_string(ReactionType type);

}  // namespace io
}  // namespace ICARION

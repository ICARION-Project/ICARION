// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "IReactionHandler.h"
#include <vector>
#include <unordered_map>

namespace ICARION {
namespace physics {

/**
 * @brief Stochastic reaction handler (discrete reaction events)
 * 
 * Implements ion-molecule reactions using Monte Carlo approach:
 * - Computes reaction probability from rate constant
 * - Samples random number to determine if reaction occurs
 * - Updates ion species/mass/charge if reaction happens
 * 
 * **SSOT Design:** Reads reactions directly from ReactionDatabase.
 * No intermediate ReactionEntry structs or parameter copies.
 * 
 * **Reaction Order Handling:**
 * Supports multi-order reactions (1st/2nd/3rd-order) through order_terms:
 * - k_eff = k₀ × ∏ᵢ [Xᵢ]^nᵢ
 * - If concentration_m3 == -1.0: uses env.particle_density_m_3 (buffer gas fallback)
 * - If concentration_m3 >= 0: uses explicit concentration (including 0 = no neutral!)
 * 
 * **Example:**
 * @code
 * // 3rd-order reaction: H3O+ + H2O + He → H5O2+ + He
 * // k_eff = 1.2e-28 [m⁶/s] × [H2O] [m⁻³] × [He] [m⁻³] = [s⁻¹]
 * 
 * StochasticReactionHandler handler(false);
 * bool rxn = handler.handle_reaction(ion, dt, rng, reaction_db, species_db, env);
 * @endcode
 * 
 * @see IReactionHandler
 * @see config::Reaction
 * @see config::ReactionOrderTerm
 */
class StochasticReactionHandler : public IReactionHandler {
public:
    /**
     * @brief Construct handler
     * @param enable_logging Enable debug logging (CSV output)
     */
    explicit StochasticReactionHandler(bool enable_logging = false);
    
    /**
     * @brief Handle reaction for single timestep (SoA)
     */
    bool handle_reaction(
        core::IonReactionData& view,
        double dt,
        PhysicsRng& rng,
        const config::ReactionDatabase& reaction_db,
        const config::SpeciesDatabase& species_db,
        const config::EnvironmentConfig& env
    ) override;
    
    std::string name() const override { return "Stochastic"; }
    
    ReactionStats get_stats() const override { return stats_; }
    void reset_stats() override;

    /**
     * @brief Compute effective rate constant for a reaction (public for GPU parity use)
     */
    double compute_effective_rate(
        const config::Reaction& reaction,
        double temperature,
        double particle_density,
        const std::unordered_map<std::string, double>& concentrations,
        const config::ReactionDatabase& reaction_db
    ) const;
    
private:
    bool enable_logging_;
    mutable ReactionStats stats_;
    
    /**
     * @brief Find applicable reactions for given ion species
     * 
     * @param species_id Current ion species (e.g., "H3O+")
     * @param reaction_db Reaction database (SSOT)
     * @return Vector of reactions where reactant matches species_id
     * 
     * **Performance:** O(N) linear search (reactions are pre-filtered in database)
     * **Future:** Could optimize with hash map (species_id → reactions)
     */
    std::vector<const config::Reaction*> find_applicable_reactions(
        const std::string& species_id,
        const config::ReactionDatabase& reaction_db
    ) const;
    
    /**
     * @brief Compute effective rate constant for reaction (private helper)
     * 
     * @param reaction Reaction definition (from database)
     * @param temperature Temperature [K] (from env)
     * @param particle_density Neutral density [m⁻³] (from env)
     * @return Effective rate constant [s⁻¹]
     * 
     * **Formula:** k_eff = k₀ × ∏ᵢ [Xᵢ]^nᵢ
     * 
     * **Order Term Handling:**
     * - If term.concentration_m3 == -1.0: Use particle_density (buffer gas fallback)
     * - If term.concentration_m3 >= 0: Use explicit concentration (including 0 = no neutral!)
     * 
     * **Example (2nd-order):**
     * - k₀ = 2.0e-9 [m³/s]
     * - [NH3] = 1e20 [m⁻³]
     * - k_eff = 2.0e-9 × 1e20 = 2e11 [s⁻¹]
     * 
     * **Example (3rd-order):**
     * - k₀ = 1.2e-28 [m⁶/s]
     * - [H2O] = 2.5e25 [m⁻³]
     * - [He] = 2.5e25 [m⁻³] (buffer gas)
     * - k_eff = 1.2e-28 × 2.5e25 × 2.5e25 = 7.5e22 [s⁻¹]
     */
    /**
     * @brief Update ion after reaction
     * 
     * @param ion Ion state (modified in-place)
     * @param product_id Product species ID (e.g., "NH4+")
     * @param species_db Species database (lookup product properties)
     * 
     * **Updates (SSOT - read from species_db):**
     * - ion.species_id = product_id
     * - ion.mass_kg = product mass
     * - ion.ion_charge_C = product charge
     * - ion.CCS_m2 = product CCS
     * - ion.reduced_mobility_cm2_Vs = product mobility
     * 
     * **Error Handling:**
     * If product species not found in database:
     * - Log warning
     * - Keep reactant species (no change)
     * - Increment failed_lookups counter
     */
    void update_ion_species(
        IonState& ion,
        const std::string& product_id,
        const config::SpeciesDatabase& species_db
    ) const;
};

} // namespace physics
} // namespace ICARION

// StochasticReactionHandler.cpp
// Implementation of stochastic reaction handler
//
// SSOT Design: Reads all parameters directly from databases and EnvironmentConfig.
// No parameter duplication!
//
// Created: 2025-11-22 (Phase 3 Refactor)

#include "StochasticReactionHandler.h"
#include "core/log/Logger.h"
#include "utils/constants.h"
#include <algorithm>
#include <cmath>

namespace ICARION {
namespace physics {

StochasticReactionHandler::StochasticReactionHandler(bool enable_logging)
    : enable_logging_(enable_logging)
{}

bool StochasticReactionHandler::handle_reaction(
    IonState& ion,
    double dt,
    EhssRng& rng,
    const config::ReactionDatabase& reaction_db,
    const config::SpeciesDatabase& species_db,
    const config::EnvironmentConfig& env
) {
    // SSOT: Read temperature/density directly from EnvironmentConfig
    const double T_K = env.temperature_K;
    const double n_m3 = env.particle_density_m_3;
    
    // Find reactions for current ion species
    auto applicable_reactions = find_applicable_reactions(ion.species_id, reaction_db);
    
    if (applicable_reactions.empty()) {
        return false;  // No reactions for this species
    }
    
    // === COMPETING CHANNELS ALGORITHM ===
    // When multiple reactions are possible, we need to:
    // 1. Compute total rate: k_total = Σ k_i
    // 2. Total probability: P_total = 1 - exp(-k_total * dt)
    // 3. Select channel with probability: P(channel i) = k_i / k_total
    
    // Step 1: Compute effective rates for all channels
    std::vector<double> k_effs;
    k_effs.reserve(applicable_reactions.size());
    double k_total = 0.0;
    
    for (const auto* reaction : applicable_reactions) {
        double k_eff = compute_effective_rate(*reaction, T_K, n_m3);
        k_effs.push_back(k_eff);
        k_total += k_eff;
    }
    
    // Early exit for negligible reaction rates
    // If k_total < 1e-60 s⁻¹, reaction probability is effectively zero
    if (k_total < 1e-60) {
        return false;  // No reaction (rate too slow)
    }
    
    // Step 2: Total reaction probability (exponential decay)
    // Numerical safety for large k_total*dt
    // If k*dt > 50, exp(-k*dt) < 2e-22 ≈ 0 → P_total ≈ 1.0
    double P_total;
    if (k_total * dt > 50.0) {
        P_total = 1.0;  // Certain reaction (avoid exp underflow)
    } else {
        P_total = 1.0 - std::exp(-k_total * dt);
    }
    
    if (rng.uniform01() >= P_total) {
        return false;  // No reaction occurs
    }
    
    // Step 3: Select reaction channel (weighted by k_eff)
    // Generate random number in [0, k_total)
    double r = rng.uniform01() * k_total;
    double cumulative = 0.0;
    
    for (size_t i = 0; i < applicable_reactions.size(); ++i) {
        cumulative += k_effs[i];
        if (r < cumulative) {
            // This reaction channel selected!
            const auto* selected_reaction = applicable_reactions[i];
            
            if (enable_logging_) {
                log::debug_log(
                    "[StochasticReactionHandler] Reaction: " + selected_reaction->reactant + 
                    " -> " + selected_reaction->product + 
                    " (k_eff=" + std::to_string(k_effs[i]) + " s⁻¹, " +
                    "k_total=" + std::to_string(k_total) + " s⁻¹)"
                );
            }
            
            // SSOT: Update ion from SpeciesDatabase
            update_ion_species(ion, selected_reaction->product, species_db);
            
            stats_.total_reactions++;
            return true;
        }
    }
    
    // Numerical safety: Should never reach here (cumulative == k_total)
    // But if we do due to floating-point errors, select last channel
    const auto* last_reaction = applicable_reactions.back();
    update_ion_species(ion, last_reaction->product, species_db);
    stats_.total_reactions++;
    return true;
}

std::vector<const config::Reaction*> StochasticReactionHandler::find_applicable_reactions(
    const std::string& species_id,
    const config::ReactionDatabase& reaction_db
) const {
    std::vector<const config::Reaction*> result;
    
    // SSOT: Read directly from reaction_db.reactions
    for (const auto& rxn : reaction_db.reactions) {
        if (rxn.reactant == species_id) {
            result.push_back(&rxn);
        }
    }
    
    return result;
}

double StochasticReactionHandler::compute_effective_rate(
    const config::Reaction& reaction,
    double temperature,
    double particle_density
) const {
    // ✅ STEP 1: Compute temperature-dependent rate constant k(T)
    // Models: Constant, Arrhenius, Modified Arrhenius
    double k_T = reaction.compute_rate_constant(temperature);
    
    // ✅ STEP 2: Apply order terms (concentration dependencies)
    // Optimization 3: Dimensional consistency check
    // ⚠️ IMPORTANT: rate_constant must have correct dimensions!
    // - 1st order (exponent=1): k [m³/s]   → k_eff = k(T) * [X]    [s⁻¹]
    // - 2nd order (exponent=2): k [m⁶/s]   → k_eff = k(T) * [X]²   [s⁻¹]
    // User is responsible for providing k with correct dimensional units!
    
    double k_eff = k_T;
    
    for (const auto& term : reaction.order_terms) {
        // Concentration handling:
        // - If concentration_m3 == -1.0: Use buffer gas density (fallback)
        // - Otherwise: Use explicit concentration (including 0 = no neutral)
        double conc_m3 = (term.concentration_m3 < 0.0)
            ? particle_density             // Fallback: buffer gas density
            : term.concentration_m3;       // Explicit concentration
        
        // k_eff *= [X]^n  (mathematically correct, but dimensional correctness depends on user!)
        k_eff *= std::pow(conc_m3, term.exponent);
    }
    
    return k_eff;  // [s⁻¹]
}

void StochasticReactionHandler::update_ion_species(
    IonState& ion,
    const std::string& product_id,
    const config::SpeciesDatabase& species_db
) const {
    // SSOT: Look up product properties in SpeciesDatabase
    if (!species_db.has(product_id)) {
        log::debug_log(
            "[StochasticReactionHandler] WARNING: Product species '" + product_id + 
            "' not found in database - keeping reactant"
        );
        stats_.failed_lookups++;
        return;
    }
    
    const auto& product = species_db.get(product_id);
    
    // Update ion properties
    ion.species_id = product_id;
    ion.mass_kg = product.mass_kg;
    ion.ion_charge_C = product.charge_C;
    ion.CCS_m2 = product.CCS_m2;
    ion.reduced_mobility_cm2_Vs = product.mobility_m2Vs / CM2_TO_M2;  // Convert back to reduced units
    
    if (enable_logging_) {
        log::debug_log(
            "[StochasticReactionHandler] Ion updated: " + product_id + 
            " mass=" + std::to_string(ion.mass_kg) + " kg"
        );
    }
}

} // namespace physics
} // namespace ICARION

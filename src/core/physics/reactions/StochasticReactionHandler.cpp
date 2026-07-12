// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "StochasticReactionHandler.h"
#include "core/log/Logger.h"
#include "utils/constants.h"
#include <algorithm>
#include <atomic>
#include <cmath>
#include <limits>
#include <mutex>
#include <unordered_set>

namespace {

constexpr double P0_PA = 1.0e5;
constexpr double PROBABILITY_SATURATION_X = 50.0;
constexpr double KEQ_EXPONENT_MIN = -700.0;
constexpr double KEQ_EXPONENT_MAX = 700.0;
constexpr double MIN_TOTAL_RATE_S_INV = 1e-60;
constexpr double MIN_VALID_EQ_RATIO = 1e-300;

double get_concentration_m3(
    const ICARION::config::ReactionOrderTerm& term,
    double particle_density,
    const std::unordered_map<std::string, double>& concentrations
) {
    double conc_m3 = term.concentration_m3;
    if (term.concentration_m3 < 0.0) {
        auto it = concentrations.find(term.species);
        if (it != concentrations.end()) {
            conc_m3 = it->second;
        } else {
            conc_m3 = (term.species == "M" || term.species == "neutral") ? particle_density : 0.0;
        }
    }
    return conc_m3;
}

double reaction_probability(double k_total, double dt) {
    const double x = k_total * dt;
    if (!std::isfinite(x) || x <= 0.0) {
        return 0.0;
    }
    if (x > PROBABILITY_SATURATION_X) {
        return 1.0;
    }
    // Numerically stable for small x: 1 - exp(-x) = -expm1(-x)
    return -std::expm1(-x);
}

void warn_reverse_dynamic_issue_once(const std::string& reaction_id, const std::string& reason) {
    static std::mutex warned_mutex;
    static std::unordered_set<std::string> warned_keys;
    const std::string key = reaction_id + "|" + reason;
    bool should_warn = false;
    {
        std::lock_guard<std::mutex> lock(warned_mutex);
        should_warn = warned_keys.insert(key).second;
    }
    if (should_warn) {
        ICARION::log::Logger::main()->warn(
            "StochasticReactionHandler: reverse_dynamic_from_equilibrium for '{}' disabled: {}",
            reaction_id,
            reason);
    }
}

const ICARION::config::Reaction* find_reaction_by_id(const ICARION::config::ReactionDatabase& reaction_db, const std::string& id) {
    for (const auto& rxn : reaction_db.reactions) {
        if (rxn.id == id) {
            return &rxn;
        }
    }
    return nullptr;
}

double compute_keq_from_thermo(double temperature_K, const ICARION::config::Reaction& reaction) {
    if (!reaction.has_thermo || temperature_K <= 0.0) {
        return std::numeric_limits<double>::quiet_NaN();
    }

    const double dH_J_mol = reaction.delta_r_H_J_mol;
    const double dS_J_molK = reaction.delta_r_S_J_molK;
    const double dG_J_mol = dH_J_mol - temperature_K * dS_J_molK;
    const double exponent = -dG_J_mol / (GAS_CONSTANT_J_MOLK * temperature_K);
    if (exponent < KEQ_EXPONENT_MIN) {
        return 0.0;
    }
    if (exponent > KEQ_EXPONENT_MAX) {
        return std::numeric_limits<double>::infinity();
    }
    return std::exp(exponent);
}

double compute_partner_pressure_ratio(
    const ICARION::config::Reaction& forward,
    double temperature_K,
    double particle_density,
    const std::unordered_map<std::string, double>& concentrations
) {
    double ratio = 1.0;
    for (const auto& term : forward.order_terms) {
        // M is a third-body placeholder and cancels from the population
        // equilibrium expression. "neutral" is rejected for equilibrium=true by
        // ReactionLoader; keep the guard here for backward compatible explicit
        // reverse channels loaded without that metadata.
        if (term.species == "M" || term.species == "neutral") {
            continue;
        }
        const double conc_m3 = std::max(0.0, get_concentration_m3(term, particle_density, concentrations));
        const double p_i_pa = conc_m3 * BOLTZMANN_CONSTANT * temperature_K;
        const double p_ratio_i = p_i_pa / P0_PA;
        ratio *= std::pow(std::max(0.0, p_ratio_i), term.exponent);
    }
    return ratio;
}

} // namespace

namespace ICARION {
namespace physics {

StochasticReactionHandler::StochasticReactionHandler(bool enable_logging)
    : enable_logging_(enable_logging)
{}

void StochasticReactionHandler::reset_stats() {
    stats_ = {};
}

bool StochasticReactionHandler::handle_reaction(
    core::IonReactionData& view,
    double dt,
    PhysicsRng& rng,
    const config::ReactionDatabase& reaction_db,
    const config::SpeciesDatabase& species_db,
    const config::EnvironmentConfig& env
) {
    IonState ion;
    ion.pos = view.kin.pos();
    ion.vel = view.kin.vel();
    ion.mass_kg = view.kin.get_mass();
    ion.ion_charge_C = view.kin.get_charge();
    ion.species_id = view.species_id();
    ion.CCS_m2 = view.get_CCS();
    ion.reduced_mobility_cm2_Vs = view.get_mobility();
    
    // SSOT: Read temperature/density directly from EnvironmentConfig
    const double T_K = env.temperature_K;
    const double n_m3 = env.particle_density_m_3;

    // Concentration map for mixture-aware rates
    std::unordered_map<std::string, double> concentrations;
    bool has_effective_reaction_component = false;
    for (const auto& comp : env.gas_mixture) {
        if (!comp.participates_in_reactions) continue;
        if (comp.density_m3 > 0.0) {
            has_effective_reaction_component = true;
        }
        concentrations[comp.species] = comp.density_m3;
    }
    if (!env.gas_mixture.empty() && !has_effective_reaction_component) {
        static std::once_flag no_active_reaction_mixture_warn_once;
        std::call_once(no_active_reaction_mixture_warn_once, []() {
            ICARION::log::Logger::main()->warn(
                "[Reactions] gas_mixture is configured but has no active reaction components; "
                "reaction rates depending on mixture species are effectively disabled until "
                "mixture flags/densities are fixed.");
        });
    }
    // Third-body and legacy buffer-gas convenience keys.
    concentrations["M"] = n_m3;
    concentrations["neutral"] = n_m3;
    // If gas_species is not explicitly listed in mixture map, provide total buffer fallback.
    if (env.gas_mixture.empty() && !env.gas_species.empty() && concentrations.find(env.gas_species) == concentrations.end()) {
        concentrations[env.gas_species] = n_m3;
    }
    
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
        // Equilibrium-linked reverse rates use the bath temperature so the
        // configured chemical population ratio is not drift energy dependent.
        double k_eff = compute_effective_rate(*reaction, T_K, n_m3, concentrations, reaction_db);
        if (!std::isfinite(k_eff) || k_eff < 0.0) {
            if (enable_logging_) {
                static std::atomic<bool> warned_invalid_rate{false};
                bool expected = false;
                if (warned_invalid_rate.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
                    log::Logger::main()->warn(
                        "StochasticReactionHandler: encountered invalid k_eff (non-finite/negative); "
                        "clamping to zero and continuing.");
                }
            }
            k_eff = 0.0;
        }
        k_effs.push_back(k_eff);
        k_total += k_eff;
    }
    
    // Early exit for negligible reaction rates
    // If k_total < MIN_TOTAL_RATE_S_INV, reaction probability is effectively zero.
    if (k_total < MIN_TOTAL_RATE_S_INV) {
        return false;  // No reaction (rate too slow)
    }
    
    // Step 2: Total reaction probability (numerically stable)
    const double P_total = reaction_probability(k_total, dt);
    
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
            view.set_species_id(ion.species_id);
            if (view.species_id() != ion.species_id) {
                static std::atomic<bool> warned_missing_species_slot{false};
                bool expected = false;
                if (warned_missing_species_slot.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
                    log::Logger::main()->warn(
                        "StochasticReactionHandler: product species '{}' is not present in IonEnsemble species map; "
                        "reaction update aborted to avoid inconsistent ion state.",
                        ion.species_id);
                }
                stats_.failed_lookups++;
                return false;
            }
            view.kin.set_mass(ion.mass_kg);
            view.kin.set_charge(ion.ion_charge_C);
            view.kin.set_pos(ion.pos);
            view.kin.set_vel(ion.vel);
            view.set_CCS(ion.CCS_m2);
            view.set_mobility(ion.reduced_mobility_cm2_Vs);

            stats_.total_reactions++;
            return true;
        }
    }
    
    // Numerical safety: Should never reach here (cumulative == k_total)
    // But if we do due to floating-point errors, select last channel
    const auto* last_reaction = applicable_reactions.back();
    update_ion_species(ion, last_reaction->product, species_db);
    view.set_species_id(ion.species_id);
    if (view.species_id() != ion.species_id) {
        static std::atomic<bool> warned_missing_species_slot{false};
        bool expected = false;
        if (warned_missing_species_slot.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
            log::Logger::main()->warn(
                "StochasticReactionHandler: product species '{}' is not present in IonEnsemble species map; "
                "reaction update aborted to avoid inconsistent ion state.",
                ion.species_id);
        }
        stats_.failed_lookups++;
        return false;
    }
    view.kin.set_mass(ion.mass_kg);
    view.kin.set_charge(ion.ion_charge_C);
    view.kin.set_pos(ion.pos);
    view.kin.set_vel(ion.vel);
    view.set_CCS(ion.CCS_m2);
    view.set_mobility(ion.reduced_mobility_cm2_Vs);
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
    double particle_density,
    const std::unordered_map<std::string, double>& concentrations,
    const config::ReactionDatabase& reaction_db
) const {
    // Dynamic reverse channel: k_reverse_eff(T) = k_forward_eff(T) / (K_eq(T) * Π p_partner/p0)
    // Forward rates remain constant as defined in DB; reverse is derived at runtime from thermo.
    if (reaction.reverse_dynamic_from_equilibrium) {
        if (reaction.linked_forward_id.empty()) {
            warn_reverse_dynamic_issue_once(reaction.id, "linked_forward_id is empty");
            return 0.0;
        }

        const auto* forward = find_reaction_by_id(reaction_db, reaction.linked_forward_id);
        if (forward == nullptr) {
            warn_reverse_dynamic_issue_once(
                reaction.id,
                "linked forward reaction '" + reaction.linked_forward_id + "' not found");
            return 0.0;
        }
        if (!forward->has_thermo) {
            warn_reverse_dynamic_issue_once(
                reaction.id,
                "linked forward reaction '" + reaction.linked_forward_id + "' has no thermo metadata");
            return 0.0;
        }

        double k_forward_T = forward->compute_rate_constant(temperature);
        double k_forward_eff = k_forward_T;
        for (const auto& term : forward->order_terms) {
            const double conc_m3 = get_concentration_m3(term, particle_density, concentrations);
            k_forward_eff *= std::pow(std::max(0.0, conc_m3), term.exponent);
        }

        const double K_eq = compute_keq_from_thermo(temperature, *forward);
        const double partner_pressure_ratio = compute_partner_pressure_ratio(
            *forward, temperature, particle_density, concentrations);
        const double eq_ratio = K_eq * partner_pressure_ratio;

        if (std::isfinite(k_forward_eff) && std::isfinite(eq_ratio) && eq_ratio > MIN_VALID_EQ_RATIO) {
            return std::max(0.0, k_forward_eff / eq_ratio);
        }
        warn_reverse_dynamic_issue_once(
            reaction.id,
            "invalid equilibrium ratio computed from linked forward reaction");
        return 0.0;
    }

    // STEP 1: Compute temperature-dependent rate constant k(T)
    // Models: Constant, Arrhenius, Modified Arrhenius
    double k_T = reaction.compute_rate_constant(temperature);
    
    // STEP 2: Apply order terms (concentration dependencies)
    // Optimization 3: Dimensional consistency check
    // IMPORTANT: rate_constant must have correct dimensions!
    // - 1st order (exponent=1): k [m³/s]   → k_eff = k(T) * [X]    [s⁻¹]
    // - 2nd order (exponent=2): k [m⁶/s]   → k_eff = k(T) * [X]²   [s⁻¹]
    // User is responsible for providing k with correct dimensional units!
    
    double k_eff = k_T;
    
    for (const auto& term : reaction.order_terms) {
        // Concentration handling:
        // - If concentration_m3 == -1.0: Use mixture concentration when available.
        //   Special cases: placeholders 'M' and 'neutral' fall back to total buffer density.
        // - Otherwise: Use explicit concentration (including 0 = no neutral)
        double conc_m3 = term.concentration_m3;
        if (term.concentration_m3 < 0.0) {
            auto it = concentrations.find(term.species);
            if (it != concentrations.end()) {
                conc_m3 = it->second;
            } else {
                static std::atomic<bool> warned{false};
                if (!warned.load(std::memory_order_relaxed) && enable_logging_) {
                    bool expected = false;
                    if (warned.compare_exchange_strong(expected, true, std::memory_order_relaxed)) {
                        log::debug_log(
                            "[StochasticReactionHandler] No concentration for '" + term.species +
                            "' in mixture; using zero concentration fallback");
                    }
                }
                conc_m3 = (term.species == "M" || term.species == "neutral") ? particle_density : 0.0;
            }
        }

        // k_eff *= [X]^n  (mathematically correct, but dimensional correctness depends on user!)
        k_eff *= std::pow(std::max(0.0, conc_m3), term.exponent);
    }

    if (!std::isfinite(k_eff) || k_eff < 0.0) {
        return 0.0;
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

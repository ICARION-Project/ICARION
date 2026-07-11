// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "HSSCollisionHandler.h"
#include "core/physics/collisions/core/CollisionKernels.h"
#include "core/physics/collisions/core/VelocitySampling.h"
#include "utils/constants.h"
#include "core/log/Logger.h"
#include <atomic>
#include <cmath>
#include <memory>
#include <mutex>
#include <stdexcept>
#include <unordered_map>
#include <vector>

namespace ICARION::physics {

namespace {
    constexpr double MIN_RELATIVE_VELOCITY = 1e-10;  ///< Minimum relative velocity to consider collision [m/s]

    struct alignas(64) HSSStatsSlot {
        std::uint64_t generation = 0;
        CollisionStats stats{};
        std::unordered_map<std::string, size_t> collisions_by_species;
    };

    struct HSSSlotCacheEntry {
        std::weak_ptr<HSSCollisionStatsState> state;
        HSSStatsSlot* slot = nullptr;
    };

    thread_local std::unordered_map<const HSSCollisionStatsState*, HSSSlotCacheEntry> hss_slot_cache;
}

struct HSSCollisionStatsState {
    mutable std::mutex mutex;
    std::atomic<std::uint64_t> generation{1};
    std::vector<std::unique_ptr<HSSStatsSlot>> slots;
    mutable std::unordered_map<std::string, size_t> species_snapshot;
};

HSSCollisionHandler::HSSCollisionHandler(bool enable_logging, const config::SpeciesDatabase* species_db)
    : enable_logging_(enable_logging)
    , species_db_(species_db)
    , stats_state_(std::make_shared<HSSCollisionStatsState>())
{}

void HSSCollisionHandler::record_collision(const std::string* gas_species) const {
    auto state = stats_state_;
    HSSStatsSlot* slot = nullptr;
    auto cache_it = hss_slot_cache.find(state.get());
    const auto cached_state = (cache_it == hss_slot_cache.end())
        ? std::shared_ptr<HSSCollisionStatsState>{}
        : cache_it->second.state.lock();
    if (cache_it == hss_slot_cache.end() || cached_state.get() != state.get()) {
        std::lock_guard<std::mutex> lock(state->mutex);
        auto owned_slot = std::make_unique<HSSStatsSlot>();
        slot = owned_slot.get();
        state->slots.push_back(std::move(owned_slot));
        cache_it = hss_slot_cache.find(state.get());
        if (cache_it == hss_slot_cache.end()) {
            cache_it = hss_slot_cache.emplace(state.get(), HSSSlotCacheEntry{state, slot}).first;
        } else {
            cache_it->second = HSSSlotCacheEntry{state, slot};
        }
    } else {
        slot = cache_it->second.slot;
    }

    const std::uint64_t generation = state->generation.load(std::memory_order_acquire);
    if (slot->generation != generation) {
        slot->generation = generation;
        slot->stats = {};
        slot->collisions_by_species.clear();
    }

    slot->stats.total_collisions++;
    if (gas_species) {
        slot->collisions_by_species[*gas_species]++;
    }
}

bool HSSCollisionHandler::mark_warning_once(const std::string& key) const {
    std::lock_guard<std::mutex> lock(warning_mutex_);
    return warned_missing_sigma_.insert(key).second;
}

bool HSSCollisionHandler::handle_collision(
    core::IonCollisionData& view,
    double dt,
    PhysicsRng& rng,
    const config::EnvironmentConfig& env,
    CollisionEventDiagnostics* diagnostics
) {
    // Build lightweight IonState view for reuse of scalar helpers
    IonState ion;
    ion.vel = view.kin.vel();
    ion.mass_kg = view.kin.get_mass();
    ion.ion_charge_C = view.kin.get_charge();
    ion.CCS_m2 = view.get_CCS();
    ion.species_id = view.species_id();

    // Mixture-aware path
    if (!env.gas_mixture.empty()) {
        // FIXED ALGORITHM (Dec 5, 2025):
        // 1. Sample neutral velocities for each component FIRST
        // 2. Calculate collision rates using ACTUAL relative velocities (not thermal average)
        // 3. This accounts for ion drift velocity increasing collision rates
        //
        // The key insight: collision rate depends on |v_ion - v_neutral|, not just <v_rel>
        // When ion drifts at high speed, collision rate increases proportionally!
        
        // Step 1: Sample neutral velocities and compute actual collision rates
        struct ComponentCollisionData {
            size_t component_idx;
            double k_i;              // Collision rate for this component
            Vec3 v_neutral;          // Sampled neutral velocity
            double sigma_i;          // Cross section
            double v_rel_mag;        // Relative speed used in rate
        };
        
        std::vector<ComponentCollisionData> collision_data;
        collision_data.reserve(env.gas_mixture.size());
        double k_total = 0.0;
        
        for (size_t i = 0; i < env.gas_mixture.size(); ++i) {
            const auto& comp = env.gas_mixture[i];
            if (!comp.participates_in_collisions) {
                continue;
            }
            double n_i = comp.density_m3;
            if (n_i <= 0.0) continue;
            
            // Get CCS
            double sigma_i = (comp.cross_section_m2 > 0.0) ? comp.cross_section_m2 : ion.CCS_m2;
            if (species_db_) {
                auto it_spec = species_db_->species.find(ion.species_id);
                if (it_spec != species_db_->species.end()) {
                    const auto& map = it_spec->second.ccs_hss_m2;
                    auto it_g = map.find(comp.species);
                    if (it_g != map.end() && it_g->second > 0.0) {
                        sigma_i = it_g->second;
                    } else {
                        // Try automatic derivation if no precomputed CCS
                        if (it_spec->second.CCS_m2 > 0.0 && 
                            it_spec->second.ccs_reference_gas.has_value()) {
                            double derived = derive_ccs_for_target_gas(
                                it_spec->second.CCS_m2,
                                *it_spec->second.ccs_reference_gas,
                                comp.species
                            );
                            if (derived > 0.0) {
                                sigma_i = derived;
                                std::string key = ion.species_id + ":" + comp.species;
                                if (mark_warning_once(key)) {
                                    ICARION::log::Logger::main()->info(
                                        "[HSS] Derived CCS for {}:{} = {:.1f} Å² (from {} reference)",
                                        ion.species_id, comp.species, derived * 1e20, 
                                        *it_spec->second.ccs_reference_gas);
                                }
                            }
                        } else if (sigma_i > 0.0) {
                            std::string key = ion.species_id + ":" + comp.species;
                            if (mark_warning_once(key)) {
                                ICARION::log::debug_log(
                                    "[HSSCollisionHandler] Warning: No CCS_HSS for gas '" + comp.species +
                                    "' and species '" + ion.species_id + "'; using ion.CCS_m2 fallback");
                            }
                        }
                    }
                }
            }
            
            if (sigma_i <= 0.0) {
                throw std::runtime_error(
                    "[HSS] Missing CCS for ion '" + ion.species_id +
                    "' and gas '" + comp.species + "'; set ion.CCS_m2 or gas cross_section_m2");
            }
            
            // Sample neutral velocity for this component
            const Vec3 v_neutral = collision_core::VelocitySampling::sample_neutral_velocity(
                env.temperature_K, comp.mass_kg, env.gas_velocity_m_s, rng
            );
            
            // Compute ACTUAL relative velocity (accounts for ion drift!)
            const Vec3 v_rel = ion.vel - v_neutral;
            const double v_rel_mag = norm(v_rel);
            
            if (v_rel_mag < MIN_RELATIVE_VELOCITY) {
                continue;  // Skip this component
            }
            
            // Collision rate using actual relative velocity (like single-gas path)
            double k_i = n_i * sigma_i * v_rel_mag;
            
            collision_data.push_back({i, k_i, v_neutral, sigma_i, v_rel_mag});
            k_total += k_i;
        }

        if (k_total <= 0.0 || collision_data.empty()) {
            return false;  // No valid collisions possible
        }

        // Step 2: Decide if collision occurs (using total rate)
        double P_total = 1.0;
        if (k_total * dt <= 50.0) {
            P_total = 1.0 - std::exp(-k_total * dt);
        }
        if (rng.uniform01() >= P_total) {
            return false;  // No collision
        }

        // Step 3: Select which component to collide with (weighted by k_i)
        double r = rng.uniform01() * k_total;
        size_t selected_idx = 0;
        double cum = 0.0;
        for (size_t j = 0; j < collision_data.size(); ++j) {
            cum += collision_data[j].k_i;
            if (r < cum) {
                selected_idx = j;
                break;
            }
        }
        
        const auto& selected = collision_data[selected_idx];
        const auto& comp = env.gas_mixture[selected.component_idx];

        // Step 4: Apply collision using the SAME neutral we sampled earlier
        // (This ensures consistency with the collision rate calculation)
        const Vec3 v_post = collision_core::CollisionKernels::hss_collision(
            ion.vel,                // Ion velocity (lab frame)
            selected.v_neutral,     // Neutral velocity we sampled earlier
            ion.mass_kg,            // Ion mass
            comp.mass_kg,           // Neutral mass
            rng                     // RNG
        );
        
        // Write back velocity to SoA view
        view.kin.set_vel(v_post);
        if (diagnostics) {
            diagnostics->v_rel_before_m_s = selected.v_rel_mag;
            diagnostics->sigma_mt_m2 = selected.sigma_i;
        }
        
        record_collision(&comp.species);
        
        return true;
    }

    // ===================================================================
    // Single-gas path (legacy HSS)
    // ===================================================================
    // ===================================================================
    // READ PARAMETERS DIRECTLY FROM ENV (SSOT!)
    // ===================================================================
    const double n = env.particle_density_m_3;
    const double T_K = env.temperature_K;
    const double m_neutral = env.gas_mass_kg;
    const Vec3 v_gas = env.gas_velocity_m_s;
    
    // Get gas-specific CCS (CRITICAL FIX: don't use ion.CCS_m2 which is for reference gas!)
    double sigma_eff = ion.CCS_m2;  // Fallback
    
    if (species_db_) {
        auto it_spec = species_db_->species.find(ion.species_id);
        if (it_spec != species_db_->species.end()) {
            const auto& map = it_spec->second.ccs_hss_m2;
            auto it_g = map.find(env.gas_species);
            if (it_g != map.end() && it_g->second > 0.0) {
                sigma_eff = it_g->second;  // Use gas-specific CCS
                static std::once_flag logged;
                std::call_once(logged, [&]() {
                    ICARION::log::Logger::main()->info(
                        "[HSS] Single-gas: Using CCS_HSS[{}][{}] = {:.1f} Å²",
                        ion.species_id, env.gas_species, sigma_eff * 1e20);
                });
            } else if (it_spec->second.CCS_m2 > 0.0 && 
                       it_spec->second.ccs_reference_gas.has_value() &&
                       *it_spec->second.ccs_reference_gas != env.gas_species) {
                // Try auto-derivation if gas doesn't match reference
                double derived = derive_ccs_for_target_gas(
                    it_spec->second.CCS_m2,
                    *it_spec->second.ccs_reference_gas,
                    env.gas_species
                );
                if (derived > 0.0) {
                    sigma_eff = derived;
                    const std::string key = ion.species_id + ":" + env.gas_species;
                    if (mark_warning_once(key)) {
                        ICARION::log::Logger::main()->info(
                            "[HSS] Single-gas path: Derived CCS for {}:{} = {:.1f} Å² (from {} reference)",
                            ion.species_id, env.gas_species, derived * 1e20,
                            *it_spec->second.ccs_reference_gas);
                    }
                }
            }
        } else {
            static std::once_flag logged;
            std::call_once(logged, [&]() {
                ICARION::log::Logger::main()->warn(
                    "[HSS] Single-gas: Species '{}' not found in database!", ion.species_id);
            });
        }
    } else {
        static std::once_flag logged;
        std::call_once(logged, []() {
            ICARION::log::Logger::main()->warn("[HSS] Single-gas: species_db_ is NULL!");
        });
    }
    
    if (sigma_eff <= 0.0) {
        return false;  // Invalid CCS
    }
    
    // Setup EHSS parameters for collision helper
    EHSSParams p;
    p.n = n;
    p.dt = dt;
    p.mi = ion.mass_kg;
    p.mn = m_neutral;
    p.kB = BOLTZMANN_CONSTANT;
    p.Tn = T_K;
    p.ubx = v_gas.x;
    p.uby = v_gas.y;
    p.ubz = v_gas.z;
    p.sigma_eff = sigma_eff;
    p.Rn = std::sqrt(sigma_eff / M_PI);  // Effective neutral radius from CCS
    
    // Sample neutral velocity from Maxwell-Boltzmann distribution FIRST
    const Vec3 v_neutral = collision_core::VelocitySampling::sample_neutral_velocity(
        T_K, m_neutral, v_gas, rng
    );
    
    // Compute relative velocity with the ACTUAL neutral we sampled
    const Vec3 v_rel = ion.vel - v_neutral;
    const double v_rel_mag = norm(v_rel);
    
    if (v_rel_mag < MIN_RELATIVE_VELOCITY) {
        return false;  // Ion stationary relative to neutral
    }
    
    // Collision probability (exponential distribution of free path)
    // P = 1 - exp(-n·σ·v_rel·dt)
    const double P = 1.0 - std::exp(-n * sigma_eff * v_rel_mag * dt);
    
    // Check if collision occurs
    if (rng.uniform01() >= P) {
        return false;  // No collision
    }
    
    // ===================================================================
    // COLLISION OCCURRED - Apply isotropic hard-sphere scattering
    // Use the SAME neutral we already sampled for consistency!
    // ==================================================================="
    
    // Apply isotropic hard-sphere collision
    const Vec3 v_post = collision_core::CollisionKernels::hss_collision(
        ion.vel,      // Ion velocity (lab frame)
        v_neutral,    // Neutral velocity (lab frame)
        ion.mass_kg,  // Ion mass
        m_neutral,    // Neutral mass
        rng           // RNG
    );
    
    // Update ion velocity
    ion.vel = v_post;
    view.kin.set_vel(ion.vel);
    if (diagnostics) {
        diagnostics->v_rel_before_m_s = v_rel_mag;
        diagnostics->sigma_mt_m2 = sigma_eff;
    }
    
    record_collision();
    
    return true;  // Collision occurred
}

CollisionStats HSSCollisionHandler::get_stats() const {
    auto state = stats_state_;
    std::lock_guard<std::mutex> lock(state->mutex);
    const std::uint64_t generation = state->generation.load(std::memory_order_acquire);

    CollisionStats reduced{};
    for (const auto& slot : state->slots) {
        if (slot && slot->generation == generation) {
            reduced.total_collisions += slot->stats.total_collisions;
            reduced.rejected_collisions += slot->stats.rejected_collisions;
            reduced.average_collision_rate += slot->stats.average_collision_rate;
        }
    }
    return reduced;
}

void HSSCollisionHandler::reset_stats() {
    auto state = stats_state_;
    std::lock_guard<std::mutex> lock(state->mutex);
    state->generation.fetch_add(1, std::memory_order_acq_rel);
    state->species_snapshot.clear();
}

const std::unordered_map<std::string, size_t>& HSSCollisionHandler::collisions_by_species() const {
    auto state = stats_state_;
    std::lock_guard<std::mutex> lock(state->mutex);
    const std::uint64_t generation = state->generation.load(std::memory_order_acquire);

    state->species_snapshot.clear();
    for (const auto& slot : state->slots) {
        if (!slot || slot->generation != generation) {
            continue;
        }
        for (const auto& [species, count] : slot->collisions_by_species) {
            state->species_snapshot[species] += count;
        }
    }
    return state->species_snapshot;
}

double HSSCollisionHandler::derive_ccs_for_target_gas(
    double sigma_ref_m2,
    const std::string& gas_ref,
    const std::string& gas_target
) const {
    // Gas radii lookup (same as EHSS)
    static const std::unordered_map<std::string, double> GAS_RADII = {
        {"He", RADIUS_HE_M},
        {"N2", RADIUS_N2_M},
        {"O2", RADIUS_O2_M},
        {"Ar", RADIUS_AR_M},
        {"CO2", RADIUS_CO2_M},
        {"Ne", RADIUS_NE_M},
        {"H2O", RADIUS_H2O_M}
    };
    
    auto it_ref = GAS_RADII.find(gas_ref);
    auto it_tgt = GAS_RADII.find(gas_target);
    
    if (it_ref == GAS_RADII.end() || it_tgt == GAS_RADII.end()) {
        return 0.0;  // Unknown gas
    }
    
    // HSS formula: extract ion radius, then compute for target gas
    double r_ref = it_ref->second;
    double r_ion = std::max(0.0, std::sqrt(sigma_ref_m2 / M_PI) - r_ref);
    double r_target = it_tgt->second;
    
    return M_PI * (r_ion + r_target) * (r_ion + r_target);
}

} // namespace ICARION::physics

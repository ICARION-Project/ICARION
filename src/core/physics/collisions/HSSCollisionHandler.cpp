// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "HSSCollisionHandler.h"
#include "core/physics/collisions/core/CollisionKernels.h"
#include "core/physics/collisions/core/VelocitySampling.h"
#include "utils/constants.h"
#include "core/log/Logger.h"
#include <cmath>
#include <vector>

namespace {
    constexpr double MIN_RELATIVE_VELOCITY = 1e-10;  ///< Minimum relative velocity to consider collision [m/s]
}

namespace ICARION::physics {

HSSCollisionHandler::HSSCollisionHandler(bool enable_logging, const config::SpeciesDatabase* species_db)
    : enable_logging_(enable_logging)
    , species_db_(species_db)
{}

bool HSSCollisionHandler::handle_collision_soa(
    core::IonCollisionData& view,
    double dt,
    PhysicsRng& rng,
    const config::EnvironmentConfig& env
) {
    // Build lightweight IonState view for reuse of scalar helpers
    IonState ion;
    ion.vel = view.kin.vel();
    ion.mass_kg = view.kin.get_mass();
    ion.ion_charge_C = view.kin.get_charge();
    ion.CCS_m2 = view.get_CCS();
    ion.species_id = view.species_id();
    
    bool occurred = handle_collision(ion, dt, rng, env);
    
    if (occurred) {
        view.kin.set_vel(ion.vel);
    }
    return occurred;
}

bool HSSCollisionHandler::handle_collision(
    IonState& ion,
    double dt,
    PhysicsRng& rng,
    const config::EnvironmentConfig& env
) {
    // Mixture-aware path
    if (!env.gas_mixture.empty()) {
        Vec3 v_rel_bulk = ion.vel - env.gas_velocity_m_s;
        double v_rel_mag = norm(v_rel_bulk);
        if (v_rel_mag < MIN_RELATIVE_VELOCITY) {
            return false;
        }

        std::vector<double> k_values;
        k_values.reserve(env.gas_mixture.size());
        double k_total = 0.0;
        for (const auto& comp : env.gas_mixture) {
            double sigma_i = (comp.cross_section_m2 > 0.0) ? comp.cross_section_m2 : ion.CCS_m2;
            if (species_db_) {
                auto it_spec = species_db_->species.find(ion.species_id);
                if (it_spec != species_db_->species.end()) {
                    const auto& map = it_spec->second.ccs_hss_m2;
                    auto it_g = map.find(comp.species);
                    if (it_g != map.end() && it_g->second > 0.0) {
                        sigma_i = it_g->second;
                        std::string key = ion.species_id + ":" + comp.species;
                        if (!warned_missing_sigma_.count(key)) {
                            ICARION::log::Logger::main()->info(
                                "[HSS] Using precomputed CCS for {}:{} = {:.1f} Å²",
                                ion.species_id, comp.species, sigma_i * 1e20);
                            warned_missing_sigma_.insert(key);
                        }
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
                                if (!warned_missing_sigma_.count(key)) {
                                    ICARION::log::Logger::main()->info(
                                        "[HSS] Derived CCS for {}:{} = {:.1f} Å² (from {} reference)",
                                        ion.species_id, comp.species, derived * 1e20, 
                                        *it_spec->second.ccs_reference_gas);
                                    warned_missing_sigma_.insert(key);
                                }
                            }
                        } else if (sigma_i > 0.0) {
                            std::string key = ion.species_id + ":" + comp.species;
                            if (!warned_missing_sigma_.count(key)) {
                                ICARION::log::debug_log(
                                    "[HSSCollisionHandler] Warning: No CCS_HSS for gas '" + comp.species +
                                    "' and species '" + ion.species_id + "'; using ion.CCS_m2 fallback");
                                warned_missing_sigma_.insert(key);
                            }
                        }
                    }
                }
            }
            double n_i = comp.density_m3;
            if (sigma_i <= 0.0 || n_i <= 0.0) {
                k_values.push_back(0.0);
                continue;
            }
            double k_i = n_i * sigma_i * v_rel_mag;
            k_values.push_back(k_i);
            k_total += k_i;
        }

        if (k_total <= 0.0) {
            throw std::runtime_error("[HSSCollisionHandler] No valid sigma in gas mixture for species '" + ion.species_id + "'");
        }

        double P_total = 1.0;
        if (k_total * dt <= 50.0) {
            P_total = 1.0 - std::exp(-k_total * dt);
        }
        if (rng.uniform01() >= P_total) {
            return false;
        }

        double r = rng.uniform01() * k_total;
        size_t idx = 0;
        double cum = 0.0;
        for (; idx < k_values.size(); ++idx) {
            cum += k_values[idx];
            if (r < cum) break;
        }
        if (idx >= env.gas_mixture.size()) {
            idx = env.gas_mixture.size() - 1;
        }
        const auto& comp = env.gas_mixture[idx];

        EHSSParams p;
        p.n = comp.density_m3;
        p.dt = dt;
        p.mi = ion.mass_kg;
        p.mn = comp.mass_kg;
        p.kB = BOLTZMANN_CONSTANT;
        p.Tn = env.temperature_K;
        p.ubx = env.gas_velocity_m_s.x;
        p.uby = env.gas_velocity_m_s.y;
        p.ubz = env.gas_velocity_m_s.z;
        // Get CCS: try precomputed, then derive, then fallback
        p.sigma_eff = 0.0;
        
        // 1. Try precomputed CCS_HSS map (BEST)
        if (species_db_) {
            auto it_spec = species_db_->species.find(ion.species_id);
            if (it_spec != species_db_->species.end()) {
                const auto& map = it_spec->second.ccs_hss_m2;
                auto it_g = map.find(comp.species);
                if (it_g != map.end() && it_g->second > 0.0) {
                    p.sigma_eff = it_g->second;  // ✅ Precomputed
                }
                // 2. Try automatic derivation (GOOD)
                else if (it_spec->second.CCS_m2 > 0.0 && 
                         it_spec->second.ccs_reference_gas.has_value()) {
                    double derived = derive_ccs_for_target_gas(
                        it_spec->second.CCS_m2,
                        *it_spec->second.ccs_reference_gas,
                        comp.species
                    );
                    if (derived > 0.0) {
                        p.sigma_eff = derived;  // ⚠️ Derived
                        std::string key = ion.species_id + ":" + comp.species;
                        if (!warned_missing_sigma_.count(key)) {
                            ICARION::log::Logger::get("collision")->info(
                                "[HSS] Derived CCS for {}:{} = {:.2f} Å² (from {} reference)",
                                ion.species_id, comp.species, derived * 1e20,
                                *it_spec->second.ccs_reference_gas
                            );
                            warned_missing_sigma_.insert(key);
                        }
                    }
                }
            }
        }
        
        // 3. Fallback: component CCS or ion reference CCS
        if (p.sigma_eff <= 0.0) {
            p.sigma_eff = (comp.cross_section_m2 > 0.0) ? comp.cross_section_m2 : ion.CCS_m2;
            if (p.sigma_eff > 0.0) {
                std::string key = ion.species_id + ":" + comp.species;
                if (!warned_missing_sigma_.count(key)) {
                    ICARION::log::Logger::get("collision")->warn(
                        "[HSS] Using fallback CCS ({:.2f} Å²) for gas {} - may be inaccurate!",
                        p.sigma_eff * 1e20, comp.species
                    );
                    warned_missing_sigma_.insert(key);
                }
            }
        }
        if (p.sigma_eff <= 0.0) {
            throw std::runtime_error("[HSSCollisionHandler] Missing CCS for species '" + ion.species_id +
                                     "' in gas '" + comp.species + "'");
        }
        p.Rn = comp.radius_m > 0.0 ? comp.radius_m : std::sqrt(std::max(p.sigma_eff, 0.0) / M_PI);

        const Vec3 v_neutral = collision_core::VelocitySampling::sample_neutral_velocity(
            env.temperature_K, comp.mass_kg, env.gas_velocity_m_s, rng
        );
        const Vec3 v_rel = ion.vel - v_neutral;
        const double v_rel_mag_actual = norm(v_rel);
        if (v_rel_mag_actual < MIN_RELATIVE_VELOCITY) {
            return false;
        }

        const Vec3 v_post = collision_core::CollisionKernels::hss_collision(
            ion.vel, v_neutral, ion.mass_kg, comp.mass_kg, rng
        );
        ion.vel = v_post;
        // For single-threaded tests, update collision statistics
        stats_.total_collisions++;
        collisions_by_species_[comp.species]++;
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
                static bool logged = false;
                if (!logged) {
                    ICARION::log::Logger::main()->info(
                        "[HSS] Single-gas: Using CCS_HSS[{}][{}] = {:.1f} Å²",
                        ion.species_id, env.gas_species, sigma_eff * 1e20);
                    logged = true;
                }
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
                    ICARION::log::Logger::main()->info(
                        "[HSS] Single-gas path: Derived CCS for {}:{} = {:.1f} Å² (from {} reference)",
                        ion.species_id, env.gas_species, derived * 1e20, 
                        *it_spec->second.ccs_reference_gas);
                }
            }
        } else {
            static bool logged = false;
            if (!logged) {
                ICARION::log::Logger::main()->warn(
                    "[HSS] Single-gas: Species '{}' not found in database!", ion.species_id);
                logged = true;
            }
        }
    } else {
        static bool logged = false;
        if (!logged) {
            ICARION::log::Logger::main()->warn("[HSS] Single-gas: species_db_ is NULL!");
            logged = true;
        }
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
    
    // Update statistics (for single-threaded tests)
    stats_.total_collisions++;
    
    return true;  // Collision occurred
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

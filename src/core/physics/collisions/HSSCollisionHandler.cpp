// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "HSSCollisionHandler.h"
#include "collisionHelpers.h"
#include "utils/constants.h"
#include <cmath>
#include <vector>

namespace {
    constexpr double MIN_RELATIVE_VELOCITY = 1e-10;  ///< Minimum relative velocity to consider collision [m/s]
}

namespace ICARION::physics {

HSSCollisionHandler::HSSCollisionHandler(bool enable_logging)
    : enable_logging_(enable_logging)
{}

bool HSSCollisionHandler::handle_collision(
    IonState& ion,
    double dt,
    EhssRng& rng,
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
            return false;
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
        p.sigma_eff = (comp.cross_section_m2 > 0.0) ? comp.cross_section_m2 : ion.CCS_m2;
        p.Rn = comp.radius_m > 0.0 ? comp.radius_m : std::sqrt(std::max(p.sigma_eff, 0.0) / M_PI);

        const Vec3 v_neutral = sample_neutral_velocity(p, rng);
        const Vec3 v_rel = ion.vel - v_neutral;
        const double v_rel_mag_actual = norm(v_rel);
        if (v_rel_mag_actual < MIN_RELATIVE_VELOCITY) {
            return false;
        }

        const Vec3 v_post = collide_hs_cpu(ion.vel, v_neutral, p, rng);
        ion.vel = v_post;
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
    
    // Use stored effective cross-section
    const double sigma_eff = ion.CCS_m2;
    
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
    const Vec3 v_neutral = sample_neutral_velocity(p, rng);
    
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
    const Vec3 v_post = collide_hs_cpu(
        ion.vel,      // Ion velocity (lab frame)
        v_neutral,    // Neutral velocity (lab frame)
        p,            // Parameters
        rng           // RNG
    );
    
    // Update ion velocity
    ion.vel = v_post;
    
    // Update statistics
    stats_.total_collisions++;
    
    return true;  // Collision occurred
}

} // namespace ICARION::physics

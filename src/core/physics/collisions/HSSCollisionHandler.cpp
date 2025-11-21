// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "HSSCollisionHandler.h"
#include "collisionHelpers.h"
#include "utils/constants.h"
#include <cmath>

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
    // ===================================================================
    // READ PARAMETERS DIRECTLY FROM ENV (SSOT!)
    // ===================================================================
    const double n = env.particle_density_m_3;
    const double T_K = env.temperature_K;
    const double m_neutral = env.neutral_mass_kg;
    const Vec3 v_gas = env.gas_velocity_m_s;
    
    // Use stored effective cross-section
    const double sigma_eff = ion.CCS_m2;
    
    // Compute relative velocity
    const Vec3 v_rel = ion.vel - v_gas;
    const double v_rel_mag = norm(v_rel);
    
    if (v_rel_mag < 1e-10 || sigma_eff <= 0.0) {
        return false;  // Ion stationary or invalid CCS
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
    // ===================================================================
    
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
    
    // Sample neutral velocity from Maxwell-Boltzmann distribution
    const Vec3 v_neutral = sample_neutral_velocity(p, rng);
    
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

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "EHSSCollisionHandler.h"
#include "collisionHelpers.h"
#include "utils/constants.h"
#include "core/io/logger.h"
#include <cmath>
#include <stdexcept>

namespace ICARION::physics {

EHSSCollisionHandler::EHSSCollisionHandler(
    const GeometryMap& geometry_map,
    bool enable_logging
)
    : geometry_map_(geometry_map)  // Store reference (no copy!)
    , enable_logging_(enable_logging)
{
    if (geometry_map_.empty()) {
        throw std::invalid_argument("EHSSCollisionHandler: geometry_map cannot be empty!");
    }
}

bool EHSSCollisionHandler::handle_collision(
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
    const double neutral_radius = env.neutral_radius_m;
    
    // Compute effective CCS from geometry (or fallback to ion.CCS_m2)
    const double sigma_eff = compute_effective_ccs(ion, neutral_radius);
    
    // Compute relative velocity
    const Vec3 v_rel = ion.vel - v_gas;
    const double v_rel_mag = norm(v_rel);
    
    if (v_rel_mag < 1e-10) {
        return false;  // Ion stationary relative to gas
    }
    
    // Collision probability (exponential distribution of free path)
    // P = 1 - exp(-n·σ·v_rel·dt)
    const double P = 1.0 - std::exp(-n * sigma_eff * v_rel_mag * dt);
    
    // Check if collision occurs
    if (rng.uniform01() >= P) {
        return false;  // No collision
    }
    
    // ===================================================================
    // COLLISION OCCURRED - Apply EHSS momentum transfer
    // ===================================================================
    
    // Sample neutral velocity from Maxwell-Boltzmann distribution
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
    p.Rn = neutral_radius;
    p.sigma_eff = sigma_eff;
    
    const Vec3 v_neutral = sample_neutral_velocity(p, rng);
    
    // Try to find geometry for this species
    auto it = geometry_map_.find(ion.species_id);
    
    Vec3 v_post;
    if (it != geometry_map_.end() && !it->second.first.empty()) {
        // Use geometry-based EHSS collision
        const auto& [centers, radii] = it->second;
        
        v_post = collide_ehss_cpu_geometry_given_neutral(
            ion.vel,      // Ion velocity (lab frame)
            v_neutral,    // Neutral velocity (lab frame)
            p,            // EHSS parameters
            centers,      // Atom positions
            radii,        // Atomic radii
            rng           // RNG
        );
    } else {
        // Fallback to isotropic hard-sphere scattering
        // (if no geometry available or geometry is empty)
        ICARION::io::debug_log(
            "[EHSSCollisionHandler] Warning: No geometry found for species '" + 
            ion.species_id + "', falling back to isotropic HSS collision"
        );
        v_post = collide_hs_cpu(
            ion.vel,      // Ion velocity (lab frame)
            v_neutral,    // Neutral velocity (lab frame)
            p,            // EHSS parameters
            rng           // RNG
        );
    }
    
    // Update ion velocity
    ion.vel = v_post;
    
    // Update statistics
    stats_.total_collisions++;
    
    return true;  // Collision occurred
}

double EHSSCollisionHandler::compute_effective_ccs(
    const IonState& ion,
    double neutral_radius
) const {
    // Try to find geometry for this species
    auto it = geometry_map_.find(ion.species_id);
    
    if (it != geometry_map_.end() && !it->second.first.empty()) {
        // Geometry available - compute CCS from atom-centered spheres
        // 
        // For now, use simple approach: sum of projected areas
        // More accurate: Monte Carlo trajectory method (TODO: Phase 3)
        
        const auto& [centers, radii] = it->second;
        
        double total_area = 0.0;
        for (size_t i = 0; i < radii.size(); ++i) {
            const double r_total = radii[i] + neutral_radius;
            total_area += M_PI * r_total * r_total;
        }
        
        return total_area;  // Effective CCS [m²]
        
    } else {
        // No geometry available - use stored CCS from ion state
        return ion.CCS_m2;
    }
}

} // namespace ICARION::physics

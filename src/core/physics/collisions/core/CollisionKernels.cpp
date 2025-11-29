// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "CollisionKernels.h"
#include "CollisionGeometry.h"
#include "core/utils/mathUtils.h"
#include "utils/constants.h"
#include <cmath>
#include <algorithm>

namespace ICARION::physics::collision_core {

Vec3 CollisionKernels::ehss_collision(
    const Vec3& v_ion_lab,
    const Vec3& v_neutral_lab,
    double ion_mass,
    double neutral_mass,
    double ion_radius,
    const std::vector<Vec3>& atom_centers,
    const std::vector<double>& atom_radii,
    EhssRng& rng,
    int max_attempts,
    double sigma_eff_m2
) {
    // PHYSICS: Exact copy from collisionHelpers.cpp::collide_ehss_cpu_geometry_given_neutral()
    // NO CHANGES to algorithm - only refactored for clarity
    
    // Validate geometry
    const int nat = static_cast<int>(atom_centers.size());
    if (nat == 0 || (int)atom_radii.size() < nat) {
        return v_ion_lab;  // Inconsistent geometry → no collision
    }
    
    // Compute relative velocity (lab frame)
    Vec3 vrel = v_ion_lab - v_neutral_lab;
    double vrel_mag = norm(vrel);
    
    // Handle near-zero relative velocity (assign small random direction)
    if (vrel_mag <= MIN_VELOCITY_MAG) {
        vrel = sample_isotropic_direction(rng) * MIN_VELOCITY_MAG;
        vrel_mag = MIN_VELOCITY_MAG;
    }
    
    Vec3 ehat = vrel / vrel_mag;  // Direction of relative velocity
    
    // Randomly rotate neutral molecule geometry
    double Rm[3][3];
    CollisionGeometry::generate_random_rotation(rng, Rm);
    
    // Rotate all atoms and compute maximum impact parameter b_max
    std::vector<Vec3> rotated_atoms(nat);
    double b_max = 0.0;
    
    for (int j = 0; j < nat; ++j) {
        Vec3 ra = CollisionGeometry::rotate_vector(atom_centers[j], Rm);
        rotated_atoms[j] = ra;
        
        // Project atom onto impact plane (perpendicular to ehat)
        double sra = dot(ehat, ra);
        Vec3 ra_perp = ra - ehat * sra;
        double ra_perp_mag = norm(ra_perp);
        
        // Maximum impact parameter for this atom
        double cur = ra_perp_mag + atom_radii[j];
        if (cur > b_max) {
            b_max = cur;
        }
    }
    
    // Add ion radius
    b_max += ion_radius;
    
    // Optionally expand b_max based on effective cross-section
    if (sigma_eff_m2 > 0.0) {
        double b_sigma = std::sqrt(sigma_eff_m2 / M_PI);
        if (b_sigma > b_max) {
            b_max = b_sigma;
        }
    }
    
    // Construct orthonormal basis for impact plane
    Vec3 t1, t2;
    CollisionGeometry::construct_orthonormal_basis(ehat, t1, t2);
    
    // Sample impact parameters until collision found
    bool hit = false;
    Vec3 n_contact{0, 0, 0};
    
    for (int attempt = 0; attempt < max_attempts && !hit; ++attempt) {
        // Sample impact parameter: b ~ U(0, b_max), φ ~ U(0, 2π)
        double uA = rng.uniform01();
        double uB = rng.uniform01();
        double b = b_max * std::sqrt(uA);
        double phi = 2.0 * M_PI * uB;
        
        // Neutral offset in impact plane
        Vec3 neutral_offset = t1 * (b * std::cos(phi)) + t2 * (b * std::sin(phi));
        
        // Check collision with each atom
        for (int j = 0; j < nat; ++j) {
            const Vec3& ra = rotated_atoms[j];
            double Rsum = atom_radii[j] + ion_radius;
            double Rsum2 = Rsum * Rsum;
            
            // Closest approach distance
            Vec3 rel = neutral_offset - ra;
            double sstar = dot(ehat, rel);
            Vec3 dvec = rel - ehat * sstar;
            double dmin2 = dot(dvec, dvec);
            
            // Collision check
            if (dmin2 <= Rsum2) {
                // Collision detected! Compute contact point
                double h = std::sqrt(std::max(0.0, Rsum2 - dmin2));
                double s_hit = sstar - h;
                Vec3 p_hit_rel = rel - ehat * s_hit;
                
                // Compute collision normal (from atom center to contact point)
                double phr2 = dot(p_hit_rel, p_hit_rel);
                if (phr2 <= MIN_CONTACT_DIST_SQ) {
                    // Degenerate contact → use -ehat as normal
                    n_contact = ehat * -1.0;
                } else {
                    n_contact = p_hit_rel * (1.0 / std::sqrt(phr2));
                    // Ensure normal points toward ion (opposite to relative velocity)
                    if (dot(vrel, n_contact) > 0.0) {
                        n_contact = n_contact * -1.0;
                    }
                }
                
                hit = true;
                break;
            }
        }
        
        // Expand b_max if not finding collisions (helps with low-probability geometries)
        if (!hit && attempt == BMAX_EXPANSION_ATTEMPT) {
            b_max *= BMAX_EXPANSION_FACTOR;
        }
    }
    
    // No collision after max_attempts
    if (!hit) {
        return v_ion_lab;
    }
    
    // Specular reflection: v' = v - 2(v·n)n
    double vdotn = dot(vrel, n_contact);
    Vec3 vrel_reflected = vrel - n_contact * (2.0 * vdotn);
    
    // Transform back to lab frame using COM recombination
    double mt = ion_mass + neutral_mass;
    double inv_mt = 1.0 / mt;
    Vec3 Vcom = (v_ion_lab * ion_mass + v_neutral_lab * neutral_mass) * inv_mt;
    
    return Vcom + vrel_reflected * (neutral_mass * inv_mt);
}

Vec3 CollisionKernels::hss_collision(
    const Vec3& v_ion_lab,
    const Vec3& v_neutral_lab,
    double ion_mass,
    double neutral_mass,
    EhssRng& rng
) {
    // PHYSICS: Exact copy from collisionHelpers.cpp::collide_hs_cpu()
    // NO CHANGES to algorithm
    
    // Compute relative velocity
    Vec3 vrel = v_ion_lab - v_neutral_lab;
    double vrel_mag = norm(vrel);
    
    // Handle near-zero relative velocity
    if (vrel_mag <= MIN_VELOCITY_MAG) {
        vrel = sample_isotropic_direction(rng) * MIN_VELOCITY_MAG;
        vrel_mag = MIN_VELOCITY_MAG;
    }
    
    // Random isotropic scattering: preserve speed, change direction
    Vec3 new_dir = sample_isotropic_direction(rng);
    Vec3 vrel_scattered = new_dir * vrel_mag;
    
    // Transform back to lab frame
    double mt = ion_mass + neutral_mass;
    double inv_mt = 1.0 / mt;
    Vec3 Vcom = (v_ion_lab * ion_mass + v_neutral_lab * neutral_mass) * inv_mt;
    
    return Vcom + vrel_scattered * (neutral_mass * inv_mt);
}

void CollisionKernels::ou_velocity_update(
    IonState& ion_state,
    EhssRng& rng,
    double dt,
    double gamma,
    double temperature_K,
    const Vec3& gas_velocity_m_s,
    bool apply_damping
) {
    // PHYSICS: Exact copy from collisionHelpers.cpp::apply_ou_velocity_kick()
    // NO CHANGES to algorithm
    
    if (gamma <= 0.0 || dt <= 0.0) {
        return;  // Invalid parameters
    }
    
    const double m = ion_state.mass_kg;
    const double kBT_over_m = BOLTZMANN_CONSTANT * temperature_K / m;
    
    const double ux = gas_velocity_m_s.x;
    const double uy = gas_velocity_m_s.y;
    const double uz = gas_velocity_m_s.z;
    
    if (apply_damping) {
        // Full OU process: v(t+dt) = u + (v(t)-u)*exp(-γ*dt) + σ*N(0,1)
        // Use when OU is the ONLY source of damping (no DampingForce)
        const double expfac = std::exp(-gamma * dt);
        const double var_coeff = kBT_over_m * (1.0 - expfac * expfac);
        const double sigma = (var_coeff > 0.0) ? std::sqrt(var_coeff) : 0.0;
        
        ion_state.vel.x = ux + (ion_state.vel.x - ux) * expfac + sigma * rng.normal();
        ion_state.vel.y = uy + (ion_state.vel.y - uy) * expfac + sigma * rng.normal();
        ion_state.vel.z = uz + (ion_state.vel.z - uz) * expfac + sigma * rng.normal();
    } else {
        // Thermal kicks only (no damping): Δv ~ N(0, √(2γkBT/m)*√dt)
        // Use when DampingForce provides continuous friction in RK4
        // This adds diffusion to match correct thermal equilibrium
        const double sigma = std::sqrt(2.0 * gamma * kBT_over_m * dt);
        
        ion_state.vel.x += sigma * rng.normal();
        ion_state.vel.y += sigma * rng.normal();
        ion_state.vel.z += sigma * rng.normal();
    }
}

// ============================================================================
// Private helper methods
// ============================================================================

Vec3 CollisionKernels::sample_isotropic_direction(EhssRng& rng) {
    // Uniform distribution on unit sphere via spherical coordinates
    // cosθ ~ U(-1,1), φ ~ U(0,2π)
    double u1 = rng.uniform01();
    double u2 = rng.uniform01();
    double cosT = 2.0 * u1 - 1.0;
    double sinT = std::sqrt(std::max(0.0, 1.0 - cosT * cosT));
    double phi = 2.0 * M_PI * u2;
    
    return Vec3{
        sinT * std::cos(phi),
        sinT * std::sin(phi),
        cosT
    };
}

Vec3 CollisionKernels::to_com_frame(
    const Vec3& v_ion_lab,
    const Vec3& v_neutral_lab,
    double ion_mass,
    double neutral_mass
) {
    double mt = ion_mass + neutral_mass;
    Vec3 Vcom = (v_ion_lab * ion_mass + v_neutral_lab * neutral_mass) / mt;
    return v_ion_lab - Vcom;
}

Vec3 CollisionKernels::from_com_frame(
    const Vec3& v_ion_com,
    const Vec3& v_neutral_lab,
    double ion_mass,
    double neutral_mass
) {
    double mt = ion_mass + neutral_mass;
    Vec3 Vcom = (v_ion_com * ion_mass + v_neutral_lab * neutral_mass) / mt;
    return Vcom + v_ion_com;
}

} // namespace ICARION::physics::collision_core

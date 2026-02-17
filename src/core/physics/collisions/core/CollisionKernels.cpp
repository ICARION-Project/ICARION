// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
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
    PhysicsRng& rng,
    int max_attempts,
    double sigma_eff_m2
) {
    double Rm[3][3];
    CollisionGeometry::generate_random_rotation(rng, Rm);

    return ehss_collision_with_orientation(
        v_ion_lab,
        v_neutral_lab,
        ion_mass,
        neutral_mass,
        ion_radius,
        atom_centers,
        atom_radii,
        Rm,
        rng,
        max_attempts,
        sigma_eff_m2,
        false
    );
}

Vec3 CollisionKernels::ehss_collision_with_orientation(
    const Vec3& v_ion_lab,
    const Vec3& v_neutral_lab,
    double ion_mass,
    double neutral_mass,
    double ion_radius,
    const std::vector<Vec3>& atom_centers,
    const std::vector<double>& atom_radii,
    const double orientation_axis_frame[3][3],
    PhysicsRng& rng,
    int max_attempts,
    double sigma_eff_m2,
    bool force_hit
) {
    // PHYSICS: Same kernel as ehss_collision, with caller-provided orientation
    // ========================================================================
    // Step 1: Validate geometry
    // ========================================================================
    const int nat = static_cast<int>(atom_centers.size());
    if (nat == 0 || static_cast<int>(atom_radii.size()) < nat) {
        return v_ion_lab;  // Inconsistent geometry → no collision
    }

    // ========================================================================
    // Step 2: Compute relative velocity and collision axis
    // ========================================================================
    Vec3 vrel = v_ion_lab - v_neutral_lab;
    double vrel_mag = norm(vrel);

    // Handle near-zero relative velocity (assign small random direction)
    if (vrel_mag <= MIN_VELOCITY_MAG) {
        vrel = sample_isotropic_direction(rng) * MIN_VELOCITY_MAG;
        vrel_mag = MIN_VELOCITY_MAG;
    }

    Vec3 collision_axis = vrel / vrel_mag;  // Direction of relative velocity

    // ========================================================================
    // Step 3: Rotate molecule with provided orientation in collision-axis frame
    // ========================================================================
    Vec3 t1, t2;
    CollisionGeometry::construct_orthonormal_basis(collision_axis, t1, t2);

    double basis[3][3] = {
        {t1.x, t2.x, collision_axis.x},
        {t1.y, t2.y, collision_axis.y},
        {t1.z, t2.z, collision_axis.z}
    };

    double Rm[3][3];
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            Rm[i][j] = basis[i][0] * orientation_axis_frame[0][j]
                     + basis[i][1] * orientation_axis_frame[1][j]
                     + basis[i][2] * orientation_axis_frame[2][j];
        }
    }

    auto molecule = rotate_and_analyze_molecule(
        atom_centers, atom_radii, collision_axis,
        ion_radius, sigma_eff_m2, Rm
    );

    // ========================================================================
    // Step 4: Sample impact parameters until collision found
    // ========================================================================
    const int attempt_limit = force_hit ? max_attempts * 16 : max_attempts;
    for (int attempt = 0; attempt < attempt_limit; ++attempt) {
        auto impact = sample_impact_geometry(
            collision_axis,
            molecule.max_impact_parameter,
            rng
        );

        auto collision = detect_atom_collision(
            molecule.atom_positions,
            atom_radii,
            impact.neutral_offset,
            collision_axis,
            ion_radius
        );

        if (collision.hit) {
            Vec3 vrel_reflected = compute_reflected_velocity(vrel, collision.contact_normal);
            return to_lab_frame(vrel_reflected, v_ion_lab, v_neutral_lab, ion_mass, neutral_mass);
        }

        if (attempt == BMAX_EXPANSION_ATTEMPT) {
            molecule.max_impact_parameter *= BMAX_EXPANSION_FACTOR;
        }
    }

    return v_ion_lab;
}

Vec3 CollisionKernels::hss_collision(
    const Vec3& v_ion_lab,
    const Vec3& v_neutral_lab,
    double ion_mass,
    double neutral_mass,
    PhysicsRng& rng
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
    PhysicsRng& rng,
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

CollisionKernels::RotatedMolecule CollisionKernels::rotate_and_analyze_molecule(
    const std::vector<Vec3>& atom_centers,
    const std::vector<double>& atom_radii,
    const Vec3& collision_axis,
    double ion_radius,
    double sigma_eff_m2,
    PhysicsRng& rng
) {
    double Rm[3][3];
    CollisionGeometry::generate_random_rotation(rng, Rm);

    return rotate_and_analyze_molecule(
        atom_centers,
        atom_radii,
        collision_axis,
        ion_radius,
        sigma_eff_m2,
        Rm
    );
}

CollisionKernels::RotatedMolecule CollisionKernels::rotate_and_analyze_molecule(
    const std::vector<Vec3>& atom_centers,
    const std::vector<double>& atom_radii,
    const Vec3& collision_axis,
    double ion_radius,
    double sigma_eff_m2,
    const double rotation[3][3]
) {
    const int nat = static_cast<int>(atom_centers.size());

    std::vector<Vec3> rotated_atoms(nat);
    double b_max = 0.0;

    for (int j = 0; j < nat; ++j) {
        Vec3 ra = CollisionGeometry::rotate_vector(atom_centers[j], rotation);
        rotated_atoms[j] = ra;

        double sra = dot(collision_axis, ra);
        Vec3 ra_perp = ra - collision_axis * sra;
        double ra_perp_mag = norm(ra_perp);

        double cur = ra_perp_mag + atom_radii[j];
        if (cur > b_max) {
            b_max = cur;
        }
    }

    b_max += ion_radius;

    if (sigma_eff_m2 > 0.0) {
        double b_sigma = std::sqrt(sigma_eff_m2 / M_PI);
        if (b_sigma > b_max) {
            b_max = b_sigma;
        }
    }

    return RotatedMolecule{std::move(rotated_atoms), b_max};
}

CollisionKernels::ImpactGeometry CollisionKernels::sample_impact_geometry(
    const Vec3& collision_axis,
    double b_max,
    PhysicsRng& rng
) {
    // Construct orthonormal basis for impact plane
    Vec3 t1, t2;
    CollisionGeometry::construct_orthonormal_basis(collision_axis, t1, t2);
    
    // Sample impact parameter: b ~ U(0, b_max), φ ~ U(0, 2π)
    double uA = rng.uniform01();
    double uB = rng.uniform01();
    double b = b_max * std::sqrt(uA);
    double phi = 2.0 * M_PI * uB;
    
    // Neutral offset in impact plane
    Vec3 neutral_offset = t1 * (b * std::cos(phi)) + t2 * (b * std::sin(phi));
    
    return ImpactGeometry{t1, t2, neutral_offset};
}

CollisionKernels::CollisionResult CollisionKernels::detect_atom_collision(
    const std::vector<Vec3>& rotated_atoms,
    const std::vector<double>& atom_radii,
    const Vec3& neutral_offset,
    const Vec3& collision_axis,
    double ion_radius
) {
    const int nat = static_cast<int>(rotated_atoms.size());
    
    // Check collision with each atom
    for (int j = 0; j < nat; ++j) {
        const Vec3& ra = rotated_atoms[j];
        double Rsum = atom_radii[j] + ion_radius;
        double Rsum2 = Rsum * Rsum;
        
        // Closest approach distance
        Vec3 rel = neutral_offset - ra;
        double sstar = dot(collision_axis, rel);
        Vec3 dvec = rel - collision_axis * sstar;
        double dmin2 = dot(dvec, dvec);
        
        // Collision check
        if (dmin2 <= Rsum2) {
            // Collision detected! Compute contact point
            double h = std::sqrt(std::max(0.0, Rsum2 - dmin2));
            double s_hit = sstar - h;
            Vec3 p_hit_rel = rel - collision_axis * s_hit;
            
            // Compute collision normal (from atom center to contact point)
            double phr2 = dot(p_hit_rel, p_hit_rel);
            Vec3 n_contact;
            
            if (phr2 <= MIN_CONTACT_DIST_SQ) {
                // Degenerate contact → use -collision_axis as normal
                n_contact = collision_axis * -1.0;
            } else {
                n_contact = p_hit_rel * (1.0 / std::sqrt(phr2));
                // Ensure normal points toward ion (opposite to relative velocity)
                // Note: relative velocity is along collision_axis
                if (dot(collision_axis, n_contact) > 0.0) {
                    n_contact = n_contact * -1.0;
                }
            }
            
            return CollisionResult{true, n_contact};
        }
    }
    
    // No collision
    return CollisionResult{false, Vec3{0, 0, 0}};
}

Vec3 CollisionKernels::compute_reflected_velocity(
    const Vec3& relative_velocity,
    const Vec3& contact_normal
) {
    // Specular reflection: v' = v - 2(v·n)n
    double vdotn = dot(relative_velocity, contact_normal);
    return relative_velocity - contact_normal * (2.0 * vdotn);
}

Vec3 CollisionKernels::to_lab_frame(
    const Vec3& v_rel_reflected,
    const Vec3& v_ion_lab,
    const Vec3& v_neutral_lab,
    double ion_mass,
    double neutral_mass
) {
    // Compute COM velocity
    double mt = ion_mass + neutral_mass;
    double inv_mt = 1.0 / mt;
    Vec3 Vcom = (v_ion_lab * ion_mass + v_neutral_lab * neutral_mass) * inv_mt;
    
    // Add reflected relative velocity (scaled by mass ratio)
    return Vcom + v_rel_reflected * (neutral_mass * inv_mt);
}

Vec3 CollisionKernels::sample_isotropic_direction(PhysicsRng& rng) {
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

} // namespace ICARION::physics::collision_core

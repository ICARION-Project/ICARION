// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file CollisionKernels.h
 * @brief Low-level collision physics kernels
 * 
 * Pure collision algorithms (NO state, NO I/O, NO configuration).
 * These are the PHYSICS CORE - DO NOT change algorithms without validation!
 * 
 * **Design:**
 * - Stateless static methods (pure functions)
 * - Input: pre-collision velocities + collision parameters
 * - Output: post-collision ion velocity
 * - All quantities in SI units (m/s, kg, K)
 * 
 * **Physics Invariants:**
 * - Momentum conservation in center-of-mass frame
 * - Energy conservation (elastic collisions)
 * - Isotropic scattering distributions (HSS)
 * - Geometry-resolved scattering (EHSS)
 * - Thermalization via Langevin dynamics (OU)
 * 
 * @note This module contains CRITICAL PHYSICS CODE. Any changes must be:
 *       1. Validated against reference data (regression tests)
 *       2. Checked for energy/momentum conservation
 *       3. Verified against thermalization benchmarks
 */

#pragma once

#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "core/physics/collisions/collisionHelpers.h"  // For EhssRng
#include <vector>

namespace ICARION::physics::collision_core {

/**
 * @brief Low-level collision physics kernels
 * 
 * Pure collision algorithms with no dependencies on configuration or I/O.
 * All methods are static and stateless.
 */
class CollisionKernels {
public:
    /**
     * @brief EHSS collision with explicit molecular geometry
     * 
     * **Physics:** Explicit Hard-Sphere Scattering with geometry-resolved atom positions.
     * 
     * Algorithm:
     * 1. Transform to center-of-mass (COM) frame
     * 2. Randomly rotate neutral molecule geometry
     * 3. Sample impact parameter in plane perpendicular to relative velocity
     * 4. Check for collision with any atom in geometry
     * 5. Compute collision normal from contact point
     * 6. Perform specular reflection in COM frame
     * 7. Transform back to lab frame
     * 
     * @param v_ion_lab Ion velocity before collision [m/s] (lab frame)
     * @param v_neutral_lab Neutral molecule velocity [m/s] (lab frame)
     * @param ion_mass Ion mass [kg]
     * @param neutral_mass Neutral molecule mass [kg]
     * @param ion_radius Ion hard-sphere radius [m]
     * @param atom_centers Atom positions [m] in molecule frame
     * @param atom_radii Atomic hard-sphere radii [m]
     * @param rng Random number generator
     * @param max_attempts Maximum impact parameter sampling attempts (default: 256)
     * @param sigma_eff_m2 Effective cross-section [m²] for b_max estimation (optional)
     * 
     * @return Ion velocity after collision [m/s] (lab frame)
     * 
     * @pre atom_centers.size() == atom_radii.size() > 0
     * @pre All radii > 0, masses > 0
     * @pre max_attempts > 0
     * 
     * @post Momentum conserved in COM frame
     * @post Energy conserved (elastic collision)
     * 
     * @note Equivalent to old `collide_ehss_cpu_geometry_given_neutral()` from collisionHelpers.h
     */
    static Vec3 ehss_collision(
        const Vec3& v_ion_lab,
        const Vec3& v_neutral_lab,
        double ion_mass,
        double neutral_mass,
        double ion_radius,
        const std::vector<Vec3>& atom_centers,
        const std::vector<double>& atom_radii,
        EhssRng& rng,
        int max_attempts = 256,
        double sigma_eff_m2 = 0.0
    );
    
    /**
     * @brief HSS collision (isotropic hard-sphere)
     * 
     * **Physics:** Simplified EHSS with single effective sphere (no geometry).
     * 
     * Algorithm:
     * 1. Transform to COM frame
     * 2. Sample random isotropic scattering direction (uniform solid angle)
     * 3. Preserve relative speed, change direction
     * 4. Transform back to lab frame
     * 
     * @param v_ion_lab Ion velocity before collision [m/s] (lab frame)
     * @param v_neutral_lab Neutral molecule velocity [m/s] (lab frame)
     * @param ion_mass Ion mass [kg]
     * @param neutral_mass Neutral molecule mass [kg]
     * @param rng Random number generator
     * 
     * @return Ion velocity after collision [m/s] (lab frame)
     * 
     * @pre masses > 0
     * 
     * @post Momentum conserved in COM frame
     * @post Energy conserved (elastic collision)
     * @post Scattering direction uniformly distributed on unit sphere
     * 
     * @note Equivalent to old `collide_hs_cpu()` from collisionHelpers.h
     */
    static Vec3 hss_collision(
        const Vec3& v_ion_lab,
        const Vec3& v_neutral_lab,
        double ion_mass,
        double neutral_mass,
        EhssRng& rng
    );
    
    /**
     * @brief Ornstein-Uhlenbeck velocity update (Langevin dynamics)
     * 
     * **Physics:** Stochastic thermalization via friction + thermal noise.
     * 
     * Full OU process: dv = -γ(v-u)dt + √(2γkBT/m) dW
     * Integrated solution: v(t+dt) = u + (v(t)-u)*exp(-γ*dt) + σ*N(0,1)
     * where σ² = (kBT/m)*(1 - exp(-2γ*dt))
     * 
     * @param ion_state Ion state (velocity will be modified in-place)
     * @param rng Random number generator
     * @param dt Timestep [s]
     * @param gamma Damping coefficient [1/s]
     * @param temperature_K Gas temperature [K]
     * @param gas_velocity_m_s Bulk gas velocity [m/s] (lab frame)
     * @param apply_damping If true, apply full OU (friction + diffusion).
     *                      If false, apply thermal kicks only (for use with DampingForce)
     * 
     * @pre gamma > 0, dt > 0, temperature_K > 0, ion_state.mass_kg > 0
     * 
     * @post Average velocity → gas_velocity_m_s (thermalization)
     * @post Velocity variance → kBT/m (thermal equilibrium)
     * 
     * @note Equivalent to old `apply_ou_velocity_kick()` from collisionHelpers.h
     */
    static void ou_velocity_update(
        IonState& ion_state,
        EhssRng& rng,
        double dt,
        double gamma,
        double temperature_K,
        const Vec3& gas_velocity_m_s,
        bool apply_damping = true
    );

private:
    // Numerical safety thresholds (from old implementation)
    static constexpr double MIN_VELOCITY_MAG = 1e-12;      // m/s
    static constexpr double MIN_CONTACT_DIST_SQ = 1e-24;  // m²
    
    // EHSS algorithm parameters
    static constexpr int DEFAULT_MAX_ATTEMPTS = 256;
    static constexpr int BMAX_EXPANSION_ATTEMPT = 64;
    static constexpr double BMAX_EXPANSION_FACTOR = 1.5;
    
    /**
     * @brief Sample random isotropic direction on unit sphere
     * 
     * Uses spherical coordinates with uniform sampling:
     * cosθ ~ U(-1,1), φ ~ U(0,2π)
     * 
     * @param rng Random number generator
     * @return Unit vector with isotropic distribution
     */
    static Vec3 sample_isotropic_direction(EhssRng& rng);
    
    /**
     * @brief Transform to center-of-mass frame
     * 
     * Computes COM velocity and returns ion velocity in COM frame.
     * 
     * @param v_ion_lab Ion velocity in lab frame
     * @param v_neutral_lab Neutral velocity in lab frame
     * @param ion_mass Ion mass
     * @param neutral_mass Neutral mass
     * @return Ion velocity in COM frame
     */
    static Vec3 to_com_frame(
        const Vec3& v_ion_lab,
        const Vec3& v_neutral_lab,
        double ion_mass,
        double neutral_mass
    );
    
    /**
     * @brief Transform from center-of-mass frame back to lab frame
     * 
     * @param v_ion_com Ion velocity in COM frame
     * @param v_neutral_lab Neutral velocity in lab frame (needed for COM velocity)
     * @param ion_mass Ion mass
     * @param neutral_mass Neutral mass
     * @return Ion velocity in lab frame
     */
    static Vec3 from_com_frame(
        const Vec3& v_ion_com,
        const Vec3& v_neutral_lab,
        double ion_mass,
        double neutral_mass
    );
};

} // namespace ICARION::physics::collision_core

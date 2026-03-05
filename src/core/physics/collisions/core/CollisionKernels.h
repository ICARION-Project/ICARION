// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file CollisionKernels.h
 * @brief Low-level collision physics kernels
 * 
 * Stateless collision routines (no I/O/config). Inputs are pre-collision
 * velocities/parameters; outputs are post-collision velocities. SI units.
 */

#pragma once

#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "core/types/CollisionTypes.h"  // For PhysicsRng
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
        PhysicsRng& rng,
        int max_attempts = 256,
        double sigma_eff_m2 = 0.0
    );

    /**
     * @brief EHSS collision with fixed molecular orientation
     *
     * Uses a caller-provided orientation (in collision-axis frame) instead of
     * drawing a random rotation internally. When force_hit is enabled, the
     * kernel resamples impact parameters until a collision is detected.
     */
    static Vec3 ehss_collision_with_orientation(
        const Vec3& v_ion_lab,
        const Vec3& v_neutral_lab,
        double ion_mass,
        double neutral_mass,
        double ion_radius,
        const std::vector<Vec3>& atom_centers,
        const std::vector<double>& atom_radii,
        const double orientation_axis_frame[3][3],
        PhysicsRng& rng,
        int max_attempts = 256,
        double sigma_eff_m2 = 0.0,
        bool force_hit = false
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
        PhysicsRng& rng
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
        PhysicsRng& rng,
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
    
    // ========================================================================
    // EHSS Internal Data Structures
    // ========================================================================
    
    /**
     * @brief Result of molecule rotation and impact parameter calculation
     */
    struct RotatedMolecule {
        std::vector<Vec3> atom_positions;  ///< Rotated atom positions
        double max_impact_parameter;        ///< Maximum impact parameter b_max
    };
    
    /**
     * @brief Impact geometry (tangent basis + neutral offset)
     */
    struct ImpactGeometry {
        Vec3 tangent1;          ///< First tangent vector (perpendicular to collision axis)
        Vec3 tangent2;          ///< Second tangent vector
        Vec3 neutral_offset;    ///< Neutral molecule offset in impact plane
    };
    
    /**
     * @brief Result of collision detection
     */
    struct CollisionResult {
        bool hit;               ///< True if collision occurred
        Vec3 contact_normal;    ///< Normal vector at contact point (points toward ion)
    };
    
    // ========================================================================
    // EHSS Helper Functions
    // ========================================================================
    
    /**
     * @brief Rotate molecule randomly and compute max impact parameter
     * 
     * Randomly rotates the molecule geometry and computes the maximum impact
     * parameter b_max by projecting all atoms onto the plane perpendicular
     * to the collision axis.
     * 
     * @param atom_centers Atom positions in molecule frame [m]
     * @param atom_radii Atomic radii [m]
     * @param collision_axis Direction of relative velocity (unit vector)
     * @param ion_radius Ion hard-sphere radius [m]
     * @param sigma_eff_m2 Optional effective cross-section [m²] for b_max expansion
     * @param rng Random number generator
     * @return Rotated molecule with b_max
     */
    static RotatedMolecule rotate_and_analyze_molecule(
        const std::vector<Vec3>& atom_centers,
        const std::vector<double>& atom_radii,
        const Vec3& collision_axis,
        double ion_radius,
        double sigma_eff_m2,
        PhysicsRng& rng
    );

    static RotatedMolecule rotate_and_analyze_molecule(
        const std::vector<Vec3>& atom_centers,
        const std::vector<double>& atom_radii,
        const Vec3& collision_axis,
        double ion_radius,
        double sigma_eff_m2,
        const double rotation[3][3]
    );
    
    /**
     * @brief Sample impact geometry (tangent basis + neutral offset)
     * 
     * Constructs orthonormal tangent basis and samples random impact parameter.
     * 
     * @param collision_axis Direction of relative velocity (unit vector)
     * @param b_max Maximum impact parameter [m]
     * @param rng Random number generator
     * @return Impact geometry
     */
    static ImpactGeometry sample_impact_geometry(
        const Vec3& collision_axis,
        double b_max,
        PhysicsRng& rng
    );
    
    /**
     * @brief Detect collision with rotated atoms
     * 
     * Checks if the ion collides with any atom in the rotated molecule.
     * If collision occurs, computes the contact normal.
     * 
     * @param rotated_atoms Rotated atom positions [m]
     * @param atom_radii Atomic radii [m]
     * @param neutral_offset Neutral molecule offset in impact plane [m]
     * @param collision_axis Direction of relative velocity (unit vector)
     * @param ion_radius Ion hard-sphere radius [m]
     * @return Collision result (hit flag + contact normal)
     */
    static CollisionResult detect_atom_collision(
        const std::vector<Vec3>& rotated_atoms,
        const std::vector<double>& atom_radii,
        const Vec3& neutral_offset,
        const Vec3& collision_axis,
        double ion_radius
    );
    
    /**
     * @brief Compute specular reflection velocity
     * 
     * Reflects relative velocity off contact normal: v' = v - 2(v·n)n
     * 
     * @param relative_velocity Relative velocity before collision [m/s]
     * @param contact_normal Normal vector at contact point (unit vector)
     * @return Reflected relative velocity [m/s]
     */
    static Vec3 compute_reflected_velocity(
        const Vec3& relative_velocity,
        const Vec3& contact_normal
    );
    
    /**
     * @brief Transform reflected relative velocity back to lab frame
     * 
     * Computes COM velocity and adds reflected relative velocity.
     * 
     * @param v_rel_reflected Reflected relative velocity [m/s] (COM frame)
     * @param v_ion_lab Ion velocity before collision [m/s] (lab frame)
     * @param v_neutral_lab Neutral velocity [m/s] (lab frame)
     * @param ion_mass Ion mass [kg]
     * @param neutral_mass Neutral mass [kg]
     * @return Ion velocity after collision [m/s] (lab frame)
     */
    static Vec3 to_lab_frame(
        const Vec3& v_rel_reflected,
        const Vec3& v_ion_lab,
        const Vec3& v_neutral_lab,
        double ion_mass,
        double neutral_mass
    );
    
    /**
     * @brief Sample random isotropic direction on unit sphere
     * 
     * Uses spherical coordinates with uniform sampling:
     * cosθ ~ U(-1,1), φ ~ U(0,2π)
     * 
     * @param rng Random number generator
     * @return Unit vector with isotropic distribution
     */
    static Vec3 sample_isotropic_direction(PhysicsRng& rng);
};

} // namespace ICARION::physics::collision_core

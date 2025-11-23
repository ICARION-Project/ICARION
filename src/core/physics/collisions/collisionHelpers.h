// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        collisionHelpers.h
 *   @brief       Core stochastic collision handling routines.
 *
 *   @details
 *   Implements ion–neutral interactions under various models:
 *   - `collide_hs_cpu`: isotropic hard-sphere scattering
 *   - `collide_ehss_cpu_geometry_given_neutral`: EHSS collisions with rotated geometry
 *   - `sampleNeutralVelocity`: Maxwell–Boltzmann sampling of neutral velocities
 *   - `rand_rotation`: random rotation matrix via uniform quaternion
 *   - `ortho_basis`: orthonormal basis perpendicular to a vector
 *
 *   Functions operate in the lab frame and transform to the center-of-mass frame.
 *   They support both deterministic and stochastic collisions and conserve momentum.
 *
 *   @note
 *   - EHSS collisions use atom-centered spheres and impact-plane sampling.
 *   - Random numbers use `EhssRng` for reproducibility.
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#pragma once

#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <random>
#include <utility>
#include <vector>
#include "core/utils/mathUtils.h"

// -----------------------------
// Parameters for EHSS
// -----------------------------
struct EHSSParams {
    double n, dt, mi, mn, kB, Tn, sigma_eff;
    double ubx, uby, ubz;
    double Rn;
    int    num_atoms    = 0;
    double max_extent   = 0.0;
    int    max_attempts = 256;
};

// -----------------------------
// RNG wrapper
// -----------------------------
class EhssRng {
  public:
    explicit EhssRng(uint64_t seed = 0xDEADBEEFCAFEBABEULL);
    double uniform01();
    double normal();

  private:
    std::mt19937_64                        eng_;
    std::uniform_real_distribution<double> uni_;
    std::normal_distribution<double>       norm_;
};

// -----------------------------
// Collision-routine helper functions
// -----------------------------

/**
 * @brief Sample neutral gas molecule velocity from Maxwell-Boltzmann distribution
 * 
 * @param p EHSS parameters containing temperature (Tn), gas flow velocity (ubx, uby, ubz)
 * @param rng Random number generator
 * 
 * @return Velocity vector [m/s] of neutral molecule in lab frame
 * 
 * Generates thermal velocity from normal distribution with width sqrt(kB*Tn/mn),
 * then adds bulk flow velocity. Used for stochastic collision modeling.
 */
Vec3 sample_neutral_velocity(const EHSSParams& p, EhssRng& rng);

/**
 * @brief Construct orthonormal basis perpendicular to given vector
 * 
 * @param ehat Input unit vector
 * @param t1 Output: first tangent vector (perpendicular to ehat)
 * @param t2 Output: second tangent vector (perpendicular to both ehat and t1)
 * 
 * Uses Gram-Schmidt orthogonalization to construct two unit vectors orthogonal
 * to ehat and each other. Useful for coordinate transformations in collision geometry.
 */
void ortho_basis(const Vec3& ehat, Vec3& t1, Vec3& t2);

/**
 * @brief Generate random rotation matrix from uniform distribution on SO(3)
 * 
 * @param rng Random number generator
 * @param R Output: 3x3 rotation matrix R[i][j]
 * 
 * Uses uniform quaternion sampling to generate isotropic random rotation.
 * Ensures unbiased orientation distribution for collision scattering.
 */
void rand_rotation(EhssRng& rng, double R[3][3]);

// -----------------------------
// Collision-routines
// -----------------------------

/**
 * @brief EHSS collision with explicit molecular geometry
 * 
 * @param v_ion_lab Ion velocity before collision [m/s] in lab frame
 * @param v_neutral_lab Neutral molecule velocity [m/s] in lab frame
 * @param p EHSS parameters (masses, temperature, cross-section, etc.)
 * @param centers Atom positions [m] in molecule frame
 * @param radii Atomic radii [m] for hard-sphere collision detection
 * @param rng Random number generator
 * 
 * @return Ion velocity after collision [m/s] in lab frame
 * 
 * Performs Elastic Hard-Sphere Scattering (EHSS) collision using atom-centered
 * spheres. Samples impact parameter and orientation, detects collision with first
 * atom encountered, then applies hard-sphere momentum transfer in COM frame.
 * Conserves momentum and energy in elastic limit.
 */
Vec3 collide_ehss_cpu_geometry_given_neutral(const Vec3& v_ion_lab, const Vec3& v_neutral_lab,
                                             const EHSSParams& p, const std::vector<Vec3>& centers,
                                             const std::vector<double>& radii, EhssRng& rng);

/**
 * @brief Hard-sphere collision with isotropic scattering
 * 
 * @param v_ion_lab Ion velocity before collision [m/s] in lab frame
 * @param v_neutral_lab Neutral molecule velocity [m/s] in lab frame
 * @param p EHSS parameters (masses, effective cross-section)
 * @param rng Random number generator
 * 
 * @return Ion velocity after collision [m/s] in lab frame
 * 
 * Simple hard-sphere collision model with random isotropic scattering angle.
 * Transforms to center-of-mass frame, applies random deflection, transforms back.
 * Conserves momentum and energy. Faster than EHSS but less physically accurate.
 */
Vec3 collide_hs_cpu(const Vec3& v_ion_lab, const Vec3& v_neutral_lab, const EHSSParams& p,
                    EhssRng& rng);

/**
 * @brief Apply Ornstein-Uhlenbeck thermal velocity kick for friction-based models
 * 
 * @param y Ion state (velocity modified in-place)
 * @param rng Random number generator
 * @param dt Time step [s]
 * @param gamma Friction coefficient [1/s]
 * @param T_K Gas temperature [K]
 * 
 * Applies stochastic thermal kick to maintain Boltzmann distribution in deterministic
 * collision models (Langevin, constant friction). Ensures correct thermal equilibrium
 * without explicit collision events. Based on Ornstein-Uhlenbeck process solution.
 */
void apply_ou_velocity_kick(IonState& y,
                            EhssRng& rng,
                            double dt,
                            double gamma,
                            double T_K);
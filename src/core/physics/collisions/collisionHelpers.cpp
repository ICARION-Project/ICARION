/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        collisionHelpers.cpp
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

#include "core/physics/collisions/collisionHelpers.h"
#include "utils/constants.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>

namespace {
    // Numerical safety thresholds
    constexpr double MIN_VELOCITY_MAG = 1e-12;  ///< Minimum velocity magnitude to avoid division by zero [m/s]
    constexpr double MIN_CONTACT_DIST_SQ = 1e-24;  ///< Minimum squared distance for contact normal calculation [m²]
    
    // EHSS algorithm parameters
    constexpr int DEFAULT_MAX_ATTEMPTS = 256;  ///< Default maximum impact parameter sampling attempts
    constexpr int BMAX_EXPANSION_ATTEMPT = 64;  ///< Attempt number at which to expand b_max by 1.5x
    constexpr double BMAX_EXPANSION_FACTOR = 1.5;  ///< Factor to expand b_max when not finding collisions
}

// -----------------------------
// EhssRng methods
// -----------------------------

/**
 * @brief Construct a reproducible random number generator for EHSS collisions.
 * @param[in] seed Seed value for deterministic behavior.
 */

EhssRng::EhssRng(uint64_t seed) : eng_(seed), uni_(0.0, 1.0), norm_(0.0, 1.0) {
}

/**
 * @brief Sample a uniform random number in [0, 1].
 * @return Uniformly distributed double in [0, 1].
 */

double EhssRng::uniform01() {
    return uni_(eng_);
}

/**
 * @brief Sample a standard normal random number.
 * @return Normally distributed double with mean 0 and stddev 1.
 */

double EhssRng::normal() {
    return norm_(eng_);
}

// -----------------------------
// Collision-routine helper functions
// -----------------------------

/**
 * @brief Sample a neutral gas particle velocity from Maxwell–Boltzmann distribution.
 *
 * Uses Box–Muller transform to generate Gaussian-distributed velocity components
 * with thermal spread determined by temperature and neutral mass.
 *
 * @param[in] p   EHSS parameters (temperature, mass, Boltzmann constant).
 * @param[in] rng Random number generator.
 * @return Velocity vector of a sampled neutral particle.
 */
Vec3 sample_neutral_velocity(const EHSSParams& p, EhssRng& rng) {
    auto normalSample = [&rng](double sigma) -> double {
        double u1 = rng.uniform01();
        double u2 = rng.uniform01();
        return sigma * std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
    };

    // Sample thermal velocity components
    Vec3 vth{normalSample(std::sqrt(p.kB * p.Tn / p.mn)),
             normalSample(std::sqrt(p.kB * p.Tn / p.mn)),
             normalSample(std::sqrt(p.kB * p.Tn / p.mn))};

    // Add bulk gas velocity to get lab-frame neutral velocity
    return Vec3{vth.x + p.ubx, vth.y + p.uby, vth.z + p.ubz};
}

/**
 * @brief Construct an orthonormal basis perpendicular to a given direction.
 *
 * Given a unit vector @p ehat, this function computes two orthonormal vectors
 * @p t1 and @p t2 such that:
 * - @f$ t_1 \perp \hat{e} @f$
 * - @f$ t_2 \perp \hat{e},\ t_2 \perp t_1 @f$
 * - @f$ \{ \hat{e}, t_1, t_2 \} @f$ forms a right-handed orthonormal basis.
 *
 * @param[in]  ehat Unit vector defining the primary direction.
 * @param[out] t1   First orthogonal vector.
 * @param[out] t2   Second orthogonal vector.
 *
 * @note This is useful for defining impact planes or scattering frames.
 */

void ortho_basis(const Vec3& ehat, Vec3& t1, Vec3& t2) {
    // pick a vector not parallel to ehat
    Vec3 a = (std::fabs(ehat.x) < 0.9) ? Vec3{1, 0, 0} : Vec3{0, 1, 0};
    t1     = normalize(cross(ehat, a));
    t2     = normalize(cross(ehat, t1));
}

/**
 * @brief Generate a random 3×3 rotation matrix using a uniform unit quaternion.
 *
 * This function samples a uniformly distributed unit quaternion using three
 * random numbers from @p rng and converts it into a rotation matrix @p R.
 * The resulting matrix is orthogonal and has determinant +1.
 *
 * @param[in,out] rng Random number generator (provides uniform samples).
 * @param[out]    R   3×3 rotation matrix representing a random rotation.
 *
 * @note The quaternion sampling follows the method of Shoemake (1992) for
 *       uniform rotation in 3D space.
 */

void rand_rotation(EhssRng& rng, double R[3][3]) {
    // Uniform unit quaternion
    double u1 = rng.uniform01();
    double u2 = rng.uniform01();
    double u3 = rng.uniform01();

    double q1 = std::sqrt(1.0 - u1) * std::sin(2.0 * M_PI * u2);
    double q2 = std::sqrt(1.0 - u1) * std::cos(2.0 * M_PI * u2);
    double q3 = std::sqrt(u1) * std::sin(2.0 * M_PI * u3);
    double q4 = std::sqrt(u1) * std::cos(2.0 * M_PI * u3);

    double x = q1, y = q2, z = q3, w = q4;
    // Rotation matrix from quaternion (x,y,z,w)
    R[0][0] = 1 - 2 * (y * y + z * z);
    R[0][1] = 2 * (x * y - z * w);
    R[0][2] = 2 * (x * z + y * w);
    R[1][0] = 2 * (x * y + z * w);
    R[1][1] = 1 - 2 * (x * x + z * z);
    R[1][2] = 2 * (y * z - x * w);
    R[2][0] = 2 * (x * z - y * w);
    R[2][1] = 2 * (y * z + x * x);
    R[2][2] = 1 - 2 * (x * x + y * y);
}

// -----------------------------
// Collision-routines
// -----------------------------

/**
 * @brief Simulate an EHSS collision between an ion and a neutral particle using explicit geometry.
 *
 * This function models a hard-sphere scattering event where the neutral particle interacts
 * with a randomly rotated molecular geometry. It samples impact parameters in the plane
 * perpendicular to the relative velocity and checks for contact with any atom in the geometry.
 * If a collision occurs, the ion velocity is updated via specular reflection in the center-of-mass
 * frame.
 *
 * @param[in]  v_ion_lab     Ion velocity in the lab frame [m/s].
 * @param[in]  v_neutral_lab Neutral velocity in the lab frame [m/s].
 * @param[in]  p             EHSS parameters (mass, radius, temperature, etc.).
 * @param[in]  centers       Atom positions in molecule frame [m].
 * @param[in]  radii         Atom radii [m].
 * @param[in,out] rng        Random number generator for sampling directions and impact parameters.
 *
 * @return Updated ion velocity in the lab frame after collision.
 *
 * @note
 * - If no collision is detected after @p p.max_attempts trials, the ion velocity remains unchanged.
 * - The geometry is randomly rotated before each collision attempt using a uniform quaternion.
 * - The collision normal is computed from the contact point and used for specular reflection.
 * - The final velocity is computed in the center-of-mass frame and transformed back to lab frame.
 */

Vec3 collide_ehss_cpu_geometry_given_neutral(const Vec3& v_ion_lab,
                                             const Vec3& v_neutral_lab,
                                             const EHSSParams& p,
                                             const std::vector<Vec3>& centers,
                                             const std::vector<double>& radii,
                                             EhssRng& rng) {
    // quick sanity
    const int nat = (p.num_atoms > 0 ? p.num_atoms : static_cast<int>(centers.size()));
    if ((int)centers.size() < nat || (int)radii.size() < nat) {
        return v_ion_lab; // inconsistent geometry -> no hit
    }

    // relative velocity and guard
    Vec3 vrel = v_ion_lab - v_neutral_lab;
    double vrel_mag = norm(vrel);
    if (vrel_mag <= MIN_VELOCITY_MAG) {
        double u1 = rng.uniform01(), u2 = rng.uniform01();
        double cosT = 2.0*u1 - 1.0;
        double sinT = std::sqrt(std::max(0.0, 1.0 - cosT*cosT));
        double phi = 2.0*M_PI*u2;
        vrel = Vec3{ sinT*std::cos(phi), sinT*std::sin(phi), cosT } * MIN_VELOCITY_MAG;
        vrel_mag = MIN_VELOCITY_MAG;
    }
    Vec3 ehat = vrel / vrel_mag;

    // random rotation
    double Rm[3][3];
    rand_rotation(rng, Rm);

    // rotated buffer (resize to nat to allow index write)
    std::vector<Vec3> rotated_atoms;
    rotated_atoms.resize(nat);

    // compute rotated atoms and b_max in one pass
    double b_max = 0.0;
    for (int j = 0; j < nat; ++j) {
        Vec3 ra = Rmul(Rm, centers[j]);
        rotated_atoms[j] = ra;
        double sra = dot(ehat, ra);
        // perpendicular component squared
        double rx = ra.x - ehat.x * sra;
        double ry = ra.y - ehat.y * sra;
        double rz = ra.z - ehat.z * sra;
        double ra_perp2 = rx*rx + ry*ry + rz*rz;
        double cur = std::sqrt(ra_perp2) + radii[j];
        if (cur > b_max) b_max = cur;
    }
    b_max += p.Rn;
    if (p.sigma_eff > 0.0) {
        double b_sigma = std::sqrt(p.sigma_eff / M_PI);
        if (b_sigma > b_max) b_max = b_sigma;
    }

    // impact-plane basis
    Vec3 t1, t2;
    ortho_basis(ehat, t1, t2);

    // sample impact parameters
    bool hit = false;
    Vec3 n_contact{0,0,0};
    int max_attempts = (p.max_attempts > 0 ? p.max_attempts : DEFAULT_MAX_ATTEMPTS);
    for (int attempt = 0; attempt < max_attempts && !hit; ++attempt) {
        double uA = rng.uniform01(), uB = rng.uniform01();
        double b = b_max * std::sqrt(uA);
        double phi = 2.0*M_PI*uB;
        Vec3 neutral_offset = t1 * (b * std::cos(phi)) + t2 * (b * std::sin(phi));

        for (int j = 0; j < nat; ++j) {
            const Vec3& ra = rotated_atoms[j];
            double Rsum = radii[j] + p.Rn;
            double Rsum2 = Rsum * Rsum;

            Vec3 rel = neutral_offset - ra;
            double sstar = dot(ehat, rel);
            Vec3 dvec = rel - ehat * sstar;
            double dmin2 = dot(dvec, dvec);

            if (dmin2 <= Rsum2) {
                double h = std::sqrt(std::max(0.0, Rsum2 - dmin2));
                double s_hit = sstar - h;
                Vec3 p_hit_rel = rel - ehat * s_hit;

                // avoid normalizing a near-zero vector
                double phr2 = dot(p_hit_rel, p_hit_rel);
                if (phr2 <= MIN_CONTACT_DIST_SQ) {
                    // degenerate contact, choose -ehat as normal
                    n_contact = ehat * -1.0;
                } else {
                    n_contact = p_hit_rel * (1.0 / std::sqrt(phr2)); // faster than normalize()
                    if (dot(vrel, n_contact) > 0.0) n_contact = n_contact * -1.0;
                }
                hit = true;
                break;
            }
        }
        if (!hit && attempt == BMAX_EXPANSION_ATTEMPT) b_max *= BMAX_EXPANSION_FACTOR;
    }

    if (!hit) return v_ion_lab;

    // reflect relative velocity about normal
    double vdotn = dot(vrel, n_contact);
    Vec3 vrel_reflected = vrel - n_contact * (2.0 * vdotn);

    // center-of-mass recombination
    double mt = p.mi + p.mn;
    double inv_mt = 1.0 / mt;
    Vec3 Vcom = (v_ion_lab * p.mi + v_neutral_lab * p.mn) * inv_mt;
    return Vcom + vrel_reflected * (p.mn * inv_mt);
}

/**
 * @brief Simulate a hard-sphere (HS) collision between an ion and a neutral particle.
 *
 * This function models an elastic collision between two spherical particles using
 * isotropic scattering. The relative velocity is reflected specularly about a randomly
 * sampled contact normal, and the post-collision ion velocity is computed in the
 * center-of-mass frame.
 *
 * @param[in]  v_ion_lab     Ion velocity in the lab frame [m/s].
 * @param[in]  v_neutral_lab Neutral velocity in the lab frame [m/s].
 * @param[in]  p             EHSSParams containing masses and radii of ion and neutral particle.
 * @param[in,out] rng        Random number generator for sampling directions.
 *
 * @return Updated ion velocity in the lab frame after collision.
 *
 * @note
 * - If the relative velocity is near zero, a fallback direction with minimal magnitude is assigned.
 * - The contact normal is sampled uniformly over the unit sphere.
 * - The collision is treated as specular reflection in the center-of-mass frame.
 * - This function assumes elastic collisions and does not model energy loss or deformation.
 */

Vec3 collide_hs_cpu(const Vec3& v_ion_lab,
                    const Vec3& v_neutral_lab,
                    const EHSSParams& p,
                    EhssRng& rng) {
    // Relative velocity
    Vec3 vrel = v_ion_lab - v_neutral_lab;
    double vrel_mag = norm(vrel);
    if (vrel_mag <= MIN_VELOCITY_MAG) {
        // assign tiny random direction
        double u1 = rng.uniform01();
        double u2 = rng.uniform01();
        double cosT = 2.0*u1 - 1.0;
        double sinT = std::sqrt(std::max(0.0, 1.0 - cosT*cosT));
        double phi  = 2.0*M_PI*u2;
        vrel = Vec3{ sinT*std::cos(phi), sinT*std::sin(phi), cosT } * MIN_VELOCITY_MAG;
        vrel_mag = MIN_VELOCITY_MAG;
    }

    // Random isotropic scattering direction
    double u1 = rng.uniform01();
    double u2 = rng.uniform01();
    double cosT = 2.0*u1 - 1.0;
    double sinT = std::sqrt(std::max(0.0, 1.0 - cosT*cosT));
    double phi  = 2.0*M_PI*u2;
    Vec3 new_dir{ sinT*std::cos(phi), sinT*std::sin(phi), cosT };

    Vec3 vrel_scattered = new_dir * vrel_mag;

    // Center-of-mass recombination
    double mt = p.mi + p.mn;
    double inv_mt = 1.0 / mt;
    Vec3 Vcom = (v_ion_lab * p.mi + v_neutral_lab * p.mn) * inv_mt;

    // Ion velocity after collision
    return Vcom + vrel_scattered * (p.mn * inv_mt);
}

void apply_ou_velocity_kick(IonState& y,
                            EhssRng& rng,
                            double dt,
                            double gamma,
                            double T_K,
                            const Vec3& gas_velocity_m_s,
                            bool apply_damping)
{
    if (gamma <= 0.0 || dt <= 0.0) return;

    const double m = y.mass_kg;
    const double kBT_over_m = BOLTZMANN_CONSTANT * T_K / m;
    
    // gas (mean) velocity passed from environment (usually [0,0,0])
    const double ux = gas_velocity_m_s.x;
    const double uy = gas_velocity_m_s.y;
    const double uz = gas_velocity_m_s.z;

    if (apply_damping) {
        // Full OU process: v(t+dt) = u + (v(t)-u)*exp(-γ*dt) + σ*N(0,1)
        // Use when OU is the ONLY source of damping (no DampingForce)
        const double expfac = std::exp(-gamma * dt);
        const double var_coeff = kBT_over_m * (1.0 - expfac * expfac);
        const double sigma = (var_coeff > 0.0) ? std::sqrt(var_coeff) : 0.0;

        y.vel.x = ux + (y.vel.x - ux) * expfac + sigma * rng.normal();
        y.vel.y = uy + (y.vel.y - uy) * expfac + sigma * rng.normal();
        y.vel.z = uz + (y.vel.z - uz) * expfac + sigma * rng.normal();
    } else {
        // Thermal kicks only (no damping): Δv ~ N(0, √(2γkBT/m)*√dt)
        // Use when DampingForce provides continuous friction in RK4
        // This adds diffusion to match correct thermal equilibrium
        const double sigma = std::sqrt(2.0 * gamma * kBT_over_m * dt);

        y.vel.x += sigma * rng.normal();
        y.vel.y += sigma * rng.normal();
        y.vel.z += sigma * rng.normal();
    }
}
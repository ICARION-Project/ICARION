#ifndef COLLISION_HS_GPU_CUH
#define COLLISION_HS_GPU_CUH

#include "core/types/Vec3.h"
#include <curand_kernel.h>
#include <cmath>

/**
 * @brief Parameters for hard-sphere collision models (GPU version).
 * 
 * This struct mirrors EHSSParams from CPU code but is suitable for GPU kernels.
 */
struct EHSSParamsGPU {
    double n;          ///< Neutral number density [m^-3]
    double dt;         ///< Time step [s]
    double mi;         ///< Ion mass [kg]
    double mn;         ///< Neutral mass [kg]
    double kB;         ///< Boltzmann constant [J/K]
    double Tn;         ///< Neutral temperature [K]
    double sigma_eff;  ///< Effective collision cross section [m^2]
    double ubx;        ///< Bulk neutral velocity X [m/s]
    double uby;        ///< Bulk neutral velocity Y [m/s]
    double ubz;        ///< Bulk neutral velocity Z [m/s]
};

/**
 * @brief Sample a neutral gas velocity from Maxwell-Boltzmann distribution (GPU).
 *
 * Uses Box-Muller transform to generate three independent Gaussian samples
 * for the neutral velocity components in the lab frame.
 *
 * @param[in] p    Collision parameters (neutral mass, temperature, bulk velocity).
 * @param[in,out] rng_state Per-ion cuRAND state.
 * @return Vec3 Neutral velocity in lab frame [m/s].
 */
__device__ inline Vec3 sample_neutral_velocity_gpu(const EHSSParamsGPU& p, curandState* rng_state) {
    // Box-Muller for Gaussian samples
    double u1 = curand_uniform_double(rng_state);
    double u2 = curand_uniform_double(rng_state);
    double u3 = curand_uniform_double(rng_state);
    double u4 = curand_uniform_double(rng_state);
    double u5 = curand_uniform_double(rng_state);
    double u6 = curand_uniform_double(rng_state);
    
    double R1 = sqrt(-2.0 * log(u1 + 1e-30));
    double R2 = sqrt(-2.0 * log(u3 + 1e-30));
    double R3 = sqrt(-2.0 * log(u5 + 1e-30));
    
    double theta1 = 2.0 * M_PI * u2;
    double theta2 = 2.0 * M_PI * u4;
    double theta3 = 2.0 * M_PI * u6;
    
    // Standard deviation for thermal velocity
    double sigma_v = sqrt(p.kB * p.Tn / p.mn);
    
    Vec3 v_thermal;
    v_thermal.x = R1 * cos(theta1) * sigma_v;
    v_thermal.y = R2 * cos(theta2) * sigma_v;
    v_thermal.z = R3 * cos(theta3) * sigma_v;
    
    // Add bulk velocity
    Vec3 v_neutral;
    v_neutral.x = v_thermal.x + p.ubx;
    v_neutral.y = v_thermal.y + p.uby;
    v_neutral.z = v_thermal.z + p.ubz;
    
    return v_neutral;
}

/**
 * @brief GPU device function for isotropic hard-sphere collision.
 *
 * Implements the same algorithm as collide_hs_cpu but uses cuRAND for random numbers.
 * Algorithm:
 * 1. Compute relative velocity (ion - neutral)
 * 2. Sample isotropic scattering direction using 2 uniform random numbers
 * 3. Apply center-of-mass transformation to get new ion velocity in lab frame
 *
 * @param[in] v_ion_lab     Ion velocity before collision [m/s]
 * @param[in] v_neutral_lab Neutral velocity [m/s]
 * @param[in] p             Collision parameters (masses)
 * @param[in,out] rng_state Per-ion cuRAND state
 * @return Vec3 Ion velocity after collision [m/s]
 */
__device__ inline Vec3 collide_hs_gpu(const Vec3& v_ion_lab, const Vec3& v_neutral_lab,
                                const EHSSParamsGPU& p, curandState* rng_state) {
    // 1. Relative velocity
    Vec3 vrel;
    vrel.x = v_ion_lab.x - v_neutral_lab.x;
    vrel.y = v_ion_lab.y - v_neutral_lab.y;
    vrel.z = v_ion_lab.z - v_neutral_lab.z;
    
    double vrel_mag = sqrt(vrel.x * vrel.x + vrel.y * vrel.y + vrel.z * vrel.z);
    
    // 2. Sample isotropic scattering direction (spherical coordinates)
    double cosT = 2.0 * curand_uniform_double(rng_state) - 1.0;
    double sinT = sqrt(1.0 - cosT * cosT);
    double phi = 2.0 * M_PI * curand_uniform_double(rng_state);
    
    Vec3 new_dir;
    new_dir.x = sinT * cos(phi);
    new_dir.y = sinT * sin(phi);
    new_dir.z = cosT;
    
    // Scattered relative velocity
    Vec3 vrel_scattered;
    vrel_scattered.x = new_dir.x * vrel_mag;
    vrel_scattered.y = new_dir.y * vrel_mag;
    vrel_scattered.z = new_dir.z * vrel_mag;
    
    // 3. Center-of-mass transformation
    double m_total = p.mi + p.mn;
    Vec3 Vcom;
    Vcom.x = (v_ion_lab.x * p.mi + v_neutral_lab.x * p.mn) / m_total;
    Vcom.y = (v_ion_lab.y * p.mi + v_neutral_lab.y * p.mn) / m_total;
    Vcom.z = (v_ion_lab.z * p.mi + v_neutral_lab.z * p.mn) / m_total;
    
    Vec3 v_ion_new;
    v_ion_new.x = Vcom.x + vrel_scattered.x * (p.mn / m_total);
    v_ion_new.y = Vcom.y + vrel_scattered.y * (p.mn / m_total);
    v_ion_new.z = Vcom.z + vrel_scattered.z * (p.mn / m_total);
    
    return v_ion_new;
}

/**
 * @brief Handle collision event on GPU (hard-sphere Monte Carlo).
 *
 * This function:
 * 1. Samples a neutral velocity from Maxwell-Boltzmann distribution
 * 2. Computes relative velocity and collision probability
 * 3. If collision occurs, updates ion velocity using isotropic hard-sphere scattering
 *
 * @param[in,out] v_ion Ion velocity (updated if collision occurs) [m/s]
 * @param[in] p Collision parameters
 * @param[in,out] rng_state Per-ion cuRAND state
 */
__device__ inline void handle_collision_hs_gpu(Vec3& v_ion, const EHSSParamsGPU& p, curandState* rng_state) {
    // 1. Sample neutral velocity
    Vec3 v_neutral = sample_neutral_velocity_gpu(p, rng_state);
    
    // 2. Compute collision probability
    Vec3 vrel;
    vrel.x = v_ion.x - v_neutral.x;
    vrel.y = v_ion.y - v_neutral.y;
    vrel.z = v_ion.z - v_neutral.z;
    double v_rel = sqrt(vrel.x * vrel.x + vrel.y * vrel.y + vrel.z * vrel.z);
    
    double P = 1.0 - exp(-p.n * p.sigma_eff * v_rel * p.dt);
    
    // 3. Check if collision occurs
    if (curand_uniform_double(rng_state) < P) {
        v_ion = collide_hs_gpu(v_ion, v_neutral, p, rng_state);
    }
}

#endif // COLLISION_HS_GPU_CUH

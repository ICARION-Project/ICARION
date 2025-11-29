#ifndef COLLISION_EHSS_GPU_CUH
#define COLLISION_EHSS_GPU_CUH

#include "core/types/Vec3.h"
#include <curand_kernel.h>
#include <cmath>

/**
 * @brief Construct an orthonormal basis perpendicular to a given direction (GPU version).
 *
 * Given a unit vector ehat, computes two orthonormal vectors t1 and t2 such that:
 * - t1 ⊥ ehat
 * - t2 ⊥ ehat, t2 ⊥ t1
 * - {ehat, t1, t2} forms a right-handed orthonormal basis
 *
 * @param[in] ehat Unit vector defining the primary direction
 * @param[out] t1 First orthogonal vector
 * @param[out] t2 Second orthogonal vector
 */
__device__ inline void ortho_basis_gpu(const Vec3& ehat, Vec3& t1, Vec3& t2) {
    // Pick a vector not parallel to ehat
    Vec3 a;
    if (fabs(ehat.x) < 0.9) {
        a.x = 1.0; a.y = 0.0; a.z = 0.0;
    } else {
        a.x = 0.0; a.y = 1.0; a.z = 0.0;
    }
    
    // t1 = normalize(cross(ehat, a))
    Vec3 cross1;
    cross1.x = ehat.y * a.z - ehat.z * a.y;
    cross1.y = ehat.z * a.x - ehat.x * a.z;
    cross1.z = ehat.x * a.y - ehat.y * a.x;
    
    double norm1 = sqrt(cross1.x * cross1.x + cross1.y * cross1.y + cross1.z * cross1.z);
    t1.x = cross1.x / norm1;
    t1.y = cross1.y / norm1;
    t1.z = cross1.z / norm1;
    
    // t2 = normalize(cross(ehat, t1))
    Vec3 cross2;
    cross2.x = ehat.y * t1.z - ehat.z * t1.y;
    cross2.y = ehat.z * t1.x - ehat.x * t1.z;
    cross2.z = ehat.x * t1.y - ehat.y * t1.x;
    
    double norm2 = sqrt(cross2.x * cross2.x + cross2.y * cross2.y + cross2.z * cross2.z);
    t2.x = cross2.x / norm2;
    t2.y = cross2.y / norm2;
    t2.z = cross2.z / norm2;
}

/**
 * @brief Generate a random 3×3 rotation matrix using a uniform unit quaternion (GPU version).
 *
 * Samples a uniformly distributed unit quaternion using three random numbers
 * and converts it into a rotation matrix R.
 *
 * @param[in,out] rng_state cuRAND state for random number generation
 * @param[out] R 3×3 rotation matrix (row-major)
 */
__device__ inline void rand_rotation_gpu(curandState* rng_state, double R[3][3]) {
    // Uniform unit quaternion (Shoemake method)
    double u1 = curand_uniform_double(rng_state);
    double u2 = curand_uniform_double(rng_state);
    double u3 = curand_uniform_double(rng_state);
    
    double q1 = sqrt(1.0 - u1) * sin(2.0 * M_PI * u2);
    double q2 = sqrt(1.0 - u1) * cos(2.0 * M_PI * u2);
    double q3 = sqrt(u1) * sin(2.0 * M_PI * u3);
    double q4 = sqrt(u1) * cos(2.0 * M_PI * u3);
    
    double x = q1, y = q2, z = q3, w = q4;
    
    // Rotation matrix from quaternion (x, y, z, w)
    R[0][0] = 1.0 - 2.0 * (y * y + z * z);
    R[0][1] = 2.0 * (x * y - z * w);
    R[0][2] = 2.0 * (x * z + y * w);
    
    R[1][0] = 2.0 * (x * y + z * w);
    R[1][1] = 1.0 - 2.0 * (x * x + z * z);
    R[1][2] = 2.0 * (y * z - x * w);
    
    R[2][0] = 2.0 * (x * z - y * w);
    R[2][1] = 2.0 * (y * z + x * w);
    R[2][2] = 1.0 - 2.0 * (x * x + y * y);
}

/**
 * @brief Multiply a 3×3 rotation matrix by a Vec3 (GPU version).
 *
 * @param[in] R 3×3 rotation matrix (row-major)
 * @param[in] v Input vector
 * @return Rotated vector R * v
 */
__device__ inline Vec3 Rmul_gpu(const double R[3][3], const Vec3& v) {
    Vec3 result;
    result.x = R[0][0] * v.x + R[0][1] * v.y + R[0][2] * v.z;
    result.y = R[1][0] * v.x + R[1][1] * v.y + R[1][2] * v.z;
    result.z = R[2][0] * v.x + R[2][1] * v.y + R[2][2] * v.z;
    return result;
}

/**
 * @brief Parameters for EHSS geometry-based collision (GPU version).
 *
 * Extended version of EHSSParamsGPU with geometry-specific fields.
 */
struct EHSSGeometryParamsGPU {
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
    double Rn;         ///< Neutral radius [m]
    int num_atoms;     ///< Number of atoms in molecule
    int max_attempts;  ///< Maximum collision attempts (default 256)
};

/**
 * @brief GPU device function for EHSS collision with molecular geometry.
 *
 * This is the GPU version of collide_ehss_cpu_geometry_given_neutral.
 * Models hard-sphere scattering with explicit atom-centered spheres.
 *
 * Algorithm:
 * 1. Compute relative velocity and direction
 * 2. Apply random rotation to molecule geometry
 * 3. Compute maximum impact parameter b_max
 * 4. Sample impact plane using random b and phi
 * 5. Check collision with each atom's sphere
 * 6. If hit, compute contact normal and reflect velocity
 * 7. Transform back to lab frame
 *
 * @param[in] v_ion_lab Ion velocity before collision [m/s]
 * @param[in] v_neutral_lab Neutral velocity [m/s]
 * @param[in] p Collision parameters with geometry info
 * @param[in] centers Atom positions in molecule frame [m] (device memory)
 * @param[in] radii Atom radii [m] (device memory)
 * @param[in,out] rng_state Per-ion cuRAND state
 * @return Vec3 Ion velocity after collision [m/s]
 */
__device__ Vec3 collide_ehss_gpu_geometry(const Vec3& v_ion_lab, const Vec3& v_neutral_lab,
                                          const EHSSGeometryParamsGPU& p,
                                          const Vec3* centers, const double* radii,
                                          curandState* rng_state) {
    const int nat = p.num_atoms;
    if (nat <= 0) return v_ion_lab;  // No geometry
    
    // 1. Relative velocity
    Vec3 vrel;
    vrel.x = v_ion_lab.x - v_neutral_lab.x;
    vrel.y = v_ion_lab.y - v_neutral_lab.y;
    vrel.z = v_ion_lab.z - v_neutral_lab.z;
    
    double vrel_mag = sqrt(vrel.x * vrel.x + vrel.y * vrel.y + vrel.z * vrel.z);
    
    // Handle near-zero relative velocity
    if (vrel_mag <= 1e-12) {
        double u1 = curand_uniform_double(rng_state);
        double u2 = curand_uniform_double(rng_state);
        double cosT = 2.0 * u1 - 1.0;
        double sinT = sqrt(fmax(0.0, 1.0 - cosT * cosT));
        double phi = 2.0 * M_PI * u2;
        vrel.x = sinT * cos(phi) * 1e-12;
        vrel.y = sinT * sin(phi) * 1e-12;
        vrel.z = cosT * 1e-12;
        vrel_mag = 1e-12;
    }
    
    // Unit vector along relative velocity
    Vec3 ehat;
    ehat.x = vrel.x / vrel_mag;
    ehat.y = vrel.y / vrel_mag;
    ehat.z = vrel.z / vrel_mag;
    
    // 2. Random rotation matrix
    double Rm[3][3];
    rand_rotation_gpu(rng_state, Rm);
    
    // 3. Rotate atoms and compute b_max
    // Note: We'll use local arrays (max 32 atoms typical for molecules)
    // For larger molecules, this could be optimized with shared memory
    Vec3 rotated_atoms[32];  // Stack allocation for small molecules
    double b_max = 0.0;
    
    for (int j = 0; j < nat && j < 32; ++j) {
        Vec3 ra = Rmul_gpu(Rm, centers[j]);
        rotated_atoms[j] = ra;
        
        // Dot product: ehat · ra
        double sra = ehat.x * ra.x + ehat.y * ra.y + ehat.z * ra.z;
        
        // Perpendicular component
        double rx = ra.x - ehat.x * sra;
        double ry = ra.y - ehat.y * sra;
        double rz = ra.z - ehat.z * sra;
        double ra_perp2 = rx * rx + ry * ry + rz * rz;
        
        double cur = sqrt(ra_perp2) + radii[j];
        if (cur > b_max) b_max = cur;
    }
    
    b_max += p.Rn;
    
    // Enlarge b_max if sigma_eff is larger
    if (p.sigma_eff > 0.0) {
        double b_sigma = sqrt(p.sigma_eff / M_PI);
        if (b_sigma > b_max) b_max = b_sigma;
    }
    
    // 4. Impact-plane basis
    Vec3 t1, t2;
    ortho_basis_gpu(ehat, t1, t2);
    
    // 5. Sample impact parameters and check for collision
    bool hit = false;
    Vec3 n_contact;
    n_contact.x = 0.0; n_contact.y = 0.0; n_contact.z = 0.0;
    
    int max_attempts = (p.max_attempts > 0) ? p.max_attempts : 256;
    
    for (int attempt = 0; attempt < max_attempts && !hit; ++attempt) {
        // Sample impact parameter
        double uA = curand_uniform_double(rng_state);
        double uB = curand_uniform_double(rng_state);
        double b = b_max * sqrt(uA);
        double phi = 2.0 * M_PI * uB;
        
        // Neutral offset in impact plane
        Vec3 neutral_offset;
        neutral_offset.x = t1.x * (b * cos(phi)) + t2.x * (b * sin(phi));
        neutral_offset.y = t1.y * (b * cos(phi)) + t2.y * (b * sin(phi));
        neutral_offset.z = t1.z * (b * cos(phi)) + t2.z * (b * sin(phi));
        
        // Check collision with each atom
        for (int j = 0; j < nat && j < 32; ++j) {
            const Vec3& ra = rotated_atoms[j];
            double Rsum = radii[j] + p.Rn;
            double Rsum2 = Rsum * Rsum;
            
            // Vector from atom center to neutral offset
            Vec3 rel;
            rel.x = neutral_offset.x - ra.x;
            rel.y = neutral_offset.y - ra.y;
            rel.z = neutral_offset.z - ra.z;
            
            // Project onto ehat direction
            double sstar = rel.x * ehat.x + rel.y * ehat.y + rel.z * ehat.z;
            
            // Perpendicular distance vector
            Vec3 dvec;
            dvec.x = rel.x - ehat.x * sstar;
            dvec.y = rel.y - ehat.y * sstar;
            dvec.z = rel.z - ehat.z * sstar;
            
            double dmin2 = dvec.x * dvec.x + dvec.y * dvec.y + dvec.z * dvec.z;
            
            if (dmin2 <= Rsum2) {
                // Collision detected!
                double h = sqrt(fmax(0.0, Rsum2 - dmin2));
                double s_hit = sstar - h;
                
                // Contact point relative to atom center
                Vec3 p_hit_rel;
                p_hit_rel.x = rel.x - ehat.x * s_hit;
                p_hit_rel.y = rel.y - ehat.y * s_hit;
                p_hit_rel.z = rel.z - ehat.z * s_hit;
                
                double phr2 = p_hit_rel.x * p_hit_rel.x + 
                             p_hit_rel.y * p_hit_rel.y + 
                             p_hit_rel.z * p_hit_rel.z;
                
                if (phr2 <= 1e-24) {
                    // Degenerate contact - use -ehat as normal
                    n_contact.x = -ehat.x;
                    n_contact.y = -ehat.y;
                    n_contact.z = -ehat.z;
                } else {
                    // Normalize to get contact normal
                    double inv_phr = 1.0 / sqrt(phr2);
                    n_contact.x = p_hit_rel.x * inv_phr;
                    n_contact.y = p_hit_rel.y * inv_phr;
                    n_contact.z = p_hit_rel.z * inv_phr;
                    
                    // Ensure normal points toward ion
                    double vdotn = vrel.x * n_contact.x + vrel.y * n_contact.y + vrel.z * n_contact.z;
                    if (vdotn > 0.0) {
                        n_contact.x = -n_contact.x;
                        n_contact.y = -n_contact.y;
                        n_contact.z = -n_contact.z;
                    }
                }
                hit = true;
                break;  // Exit atom loop
            }
        }
        
        // Expand search radius if no hit after many attempts
        if (!hit && attempt == 64) {
            b_max *= 1.5;
        }
    }
    
    // If no hit, return original velocity
    if (!hit) return v_ion_lab;
    
    // 6. Reflect relative velocity about contact normal (specular reflection)
    double vdotn = vrel.x * n_contact.x + vrel.y * n_contact.y + vrel.z * n_contact.z;
    Vec3 vrel_reflected;
    vrel_reflected.x = vrel.x - n_contact.x * (2.0 * vdotn);
    vrel_reflected.y = vrel.y - n_contact.y * (2.0 * vdotn);
    vrel_reflected.z = vrel.z - n_contact.z * (2.0 * vdotn);
    
    // 7. Center-of-mass transformation
    double mt = p.mi + p.mn;
    double inv_mt = 1.0 / mt;
    
    Vec3 Vcom;
    Vcom.x = (v_ion_lab.x * p.mi + v_neutral_lab.x * p.mn) * inv_mt;
    Vcom.y = (v_ion_lab.y * p.mi + v_neutral_lab.y * p.mn) * inv_mt;
    Vcom.z = (v_ion_lab.z * p.mi + v_neutral_lab.z * p.mn) * inv_mt;
    
    Vec3 v_ion_new;
    v_ion_new.x = Vcom.x + vrel_reflected.x * (p.mn * inv_mt);
    v_ion_new.y = Vcom.y + vrel_reflected.y * (p.mn * inv_mt);
    v_ion_new.z = Vcom.z + vrel_reflected.z * (p.mn * inv_mt);
    
    return v_ion_new;
}

#endif // COLLISION_EHSS_GPU_CUH

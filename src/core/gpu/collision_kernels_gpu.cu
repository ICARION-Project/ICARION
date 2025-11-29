// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file collision_kernels_gpu.cu
 * @brief GPU collision kernel implementations
 */

#include "collision_kernels_gpu.cuh"
#include <cmath>

namespace ICARION {
namespace gpu {

// Physical constants
constexpr double BOLTZMANN_CONSTANT = 1.380649e-23; // J/K
constexpr double PI = 3.14159265358979323846;
constexpr double TWO_PI = 2.0 * PI;

// Numerical safety thresholds
constexpr double MIN_VELOCITY_MAG = 1e-12;  // m/s
constexpr double MIN_CONTACT_DIST_SQ = 1e-24;  // m²

// EHSS parameters
constexpr int MAX_COLLISION_ATTEMPTS = 256;
constexpr double BMAX_EXPANSION_FACTOR = 1.5;

// ============================================================================
// Device Helper Functions
// ============================================================================

/**
 * @brief Sample velocity component from Maxwell-Boltzmann distribution
 * 
 * Uses Box-Muller transform: v ~ N(0, sqrt(kT/m))
 */
__device__ inline double sample_maxwell_boltzmann_component(
    curandState* state,
    double temperature_K,
    double mass_kg
) {
    double sigma = sqrt(BOLTZMANN_CONSTANT * temperature_K / mass_kg);
    double u1 = curand_uniform_double(state);
    double u2 = curand_uniform_double(state);
    return sigma * sqrt(-2.0 * log(u1)) * cos(TWO_PI * u2);
}

/**
 * @brief Sample isotropic unit vector (uniform on sphere)
 * 
 * Marsaglia (1972) method: rejection sampling on unit sphere
 */
__device__ inline void sample_isotropic_direction(
    curandState* state,
    double& x,
    double& y,
    double& z
) {
    double s;
    do {
        x = 2.0 * curand_uniform_double(state) - 1.0;
        y = 2.0 * curand_uniform_double(state) - 1.0;
        s = x * x + y * y;
    } while (s >= 1.0 || s < 1e-10);
    
    double factor = sqrt((1.0 - s) / s);
    x *= factor;
    y *= factor;
    z = 1.0 - 2.0 * s;
}

/**
 * @brief Normalize vector in-place
 */
__device__ inline void normalize(double& x, double& y, double& z) {
    double mag = sqrt(x * x + y * y + z * z);
    if (mag > MIN_VELOCITY_MAG) {
        x /= mag;
        y /= mag;
        z /= mag;
    }
}

/**
 * @brief Dot product
 */
__device__ inline double dot(double x1, double y1, double z1,
                              double x2, double y2, double z2) {
    return x1 * x2 + y1 * y2 + z1 * z2;
}

/**
 * @brief Cross product: (x3, y3, z3) = (x1, y1, z1) × (x2, y2, z2)
 */
__device__ inline void cross(double x1, double y1, double z1,
                              double x2, double y2, double z2,
                              double& x3, double& y3, double& z3) {
    x3 = y1 * z2 - z1 * y2;
    y3 = z1 * x2 - x1 * z2;
    z3 = x1 * y2 - y1 * x2;
}

/**
 * @brief Rotate vector by Euler angles (Z-Y-Z convention)
 */
__device__ inline void rotate_euler_zyz(
    double x, double y, double z,
    double alpha, double beta, double gamma,
    double& x_rot, double& y_rot, double& z_rot
) {
    // Rotation matrix elements (Z-Y-Z Euler angles)
    double ca = cos(alpha), sa = sin(alpha);
    double cb = cos(beta),  sb = sin(beta);
    double cg = cos(gamma), sg = sin(gamma);
    
    double r11 = ca * cb * cg - sa * sg;
    double r12 = -ca * cb * sg - sa * cg;
    double r13 = ca * sb;
    double r21 = sa * cb * cg + ca * sg;
    double r22 = -sa * cb * sg + ca * cg;
    double r23 = sa * sb;
    double r31 = -sb * cg;
    double r32 = sb * sg;
    double r33 = cb;
    
    x_rot = r11 * x + r12 * y + r13 * z;
    y_rot = r21 * x + r22 * y + r23 * z;
    z_rot = r31 * x + r32 * y + r33 * z;
}

// ============================================================================
// Kernel Implementations
// ============================================================================

__global__ void init_curand_states(
    curandState* states,
    unsigned long long seed,
    int n_threads
) {
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    if (tid < n_threads) {
        curand_init(seed, tid, 0, &states[tid]);
    }
}

__global__ void hss_collision_kernel(
    double* vx_inout,
    double* vy_inout,
    double* vz_inout,
    const double* mass,
    const double* ccs,
    const bool* active,
    curandState* curand_states,
    const EnvironmentParams_GPU env,
    double dt,
    int n_ions
) {
    // Grid-stride loop
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;
    
    // Get thread-local RNG state
    curandState local_state = curand_states[tid];
    
    for (int i = tid; i < n_ions; i += stride) {
        // Skip inactive ions
        if (!active[i]) continue;
        
        // Read ion state
        double vx = vx_inout[i];
        double vy = vy_inout[i];
        double vz = vz_inout[i];
        double m_ion = mass[i];
        double sigma = ccs[i];
        
        // Sample neutral velocity from Maxwell-Boltzmann
        double vn_x = env.gas_velocity_x + sample_maxwell_boltzmann_component(
            &local_state, env.temperature_K, env.neutral_mass_kg);
        double vn_y = env.gas_velocity_y + sample_maxwell_boltzmann_component(
            &local_state, env.temperature_K, env.neutral_mass_kg);
        double vn_z = env.gas_velocity_z + sample_maxwell_boltzmann_component(
            &local_state, env.temperature_K, env.neutral_mass_kg);
        
        // Relative velocity
        double vrel_x = vx - vn_x;
        double vrel_y = vy - vn_y;
        double vrel_z = vz - vn_z;
        double vrel_mag = sqrt(vrel_x * vrel_x + vrel_y * vrel_y + vrel_z * vrel_z);
        
        if (vrel_mag < MIN_VELOCITY_MAG) continue;
        
        // Collision probability: P = 1 - exp(-v_rel * σ * n * dt / 4)
        // where n = P / (kT) is number density
        double number_density = env.pressure_Pa / (BOLTZMANN_CONSTANT * env.temperature_K);
        double collision_rate = 0.25 * vrel_mag * sigma * number_density;
        double collision_prob = 1.0 - exp(-collision_rate * dt);
        
        // Check if collision occurs
        double u = curand_uniform_double(&local_state);
        if (u >= collision_prob) continue;
        
        // ====== COLLISION OCCURS ======
        
        // Transform to center-of-mass frame
        double m_neutral = env.neutral_mass_kg;
        double m_total = m_ion + m_neutral;
        double vcm_x = (m_ion * vx + m_neutral * vn_x) / m_total;
        double vcm_y = (m_ion * vy + m_neutral * vn_y) / m_total;
        double vcm_z = (m_ion * vz + m_neutral * vn_z) / m_total;
        
        double vion_cm_x = vx - vcm_x;
        double vion_cm_y = vy - vcm_y;
        double vion_cm_z = vz - vcm_z;
        
        // Sample isotropic scattering direction (uniform sphere)
        double nx, ny, nz;
        sample_isotropic_direction(&local_state, nx, ny, nz);
        
        // Preserve relative speed magnitude, change direction
        double vion_cm_mag = sqrt(vion_cm_x * vion_cm_x + 
                                   vion_cm_y * vion_cm_y + 
                                   vion_cm_z * vion_cm_z);
        
        vion_cm_x = nx * vion_cm_mag;
        vion_cm_y = ny * vion_cm_mag;
        vion_cm_z = nz * vion_cm_mag;
        
        // Transform back to lab frame
        vx = vion_cm_x + vcm_x;
        vy = vion_cm_y + vcm_y;
        vz = vion_cm_z + vcm_z;
        
        // Write back
        vx_inout[i] = vx;
        vy_inout[i] = vy;
        vz_inout[i] = vz;
    }
    
    // Save RNG state
    curand_states[tid] = local_state;
}

__global__ void ehss_collision_kernel(
    double* vx_inout,
    double* vy_inout,
    double* vz_inout,
    const double* mass,
    const double* ccs,
    const int* species_indices,
    const bool* active,
    curandState* curand_states,
    const EnvironmentParams_GPU env,
    const GeometryData_GPU geometry,
    double dt,
    int n_ions
) {
    // Grid-stride loop
    int tid = blockIdx.x * blockDim.x + threadIdx.x;
    int stride = gridDim.x * blockDim.x;
    
    curandState local_state = curand_states[tid];
    
    for (int i = tid; i < n_ions; i += stride) {
        if (!active[i]) continue;
        
        // Read ion state
        double vx = vx_inout[i];
        double vy = vy_inout[i];
        double vz = vz_inout[i];
        double m_ion = mass[i];
        double sigma = ccs[i];
        int species_idx = species_indices[i];
        
        // Sample neutral velocity
        double vn_x = env.gas_velocity_x + sample_maxwell_boltzmann_component(
            &local_state, env.temperature_K, env.neutral_mass_kg);
        double vn_y = env.gas_velocity_y + sample_maxwell_boltzmann_component(
            &local_state, env.temperature_K, env.neutral_mass_kg);
        double vn_z = env.gas_velocity_z + sample_maxwell_boltzmann_component(
            &local_state, env.temperature_K, env.neutral_mass_kg);
        
        // Relative velocity
        double vrel_x = vx - vn_x;
        double vrel_y = vy - vn_y;
        double vrel_z = vz - vn_z;
        double vrel_mag = sqrt(vrel_x * vrel_x + vrel_y * vrel_y + vrel_z * vrel_z);
        
        if (vrel_mag < MIN_VELOCITY_MAG) continue;
        
        // Collision probability
        double number_density = env.pressure_Pa / (BOLTZMANN_CONSTANT * env.temperature_K);
        double collision_rate = 0.25 * vrel_mag * sigma * number_density;
        double collision_prob = 1.0 - exp(-collision_rate * dt);
        
        double u = curand_uniform_double(&local_state);
        if (u >= collision_prob) continue;
        
        // ====== COLLISION WITH GEOMETRY ======
        
        // Look up geometry for this species
        if (species_idx < 0 || species_idx >= geometry.num_species) {
            // Fallback to HSS if no geometry
            double m_neutral = env.neutral_mass_kg;
            double m_total = m_ion + m_neutral;
            double vcm_x = (m_ion * vx + m_neutral * vn_x) / m_total;
            double vcm_y = (m_ion * vy + m_neutral * vn_y) / m_total;
            double vcm_z = (m_ion * vz + m_neutral * vn_z) / m_total;
            
            double vion_cm_x = vx - vcm_x;
            double vion_cm_y = vy - vcm_y;
            double vion_cm_z = vz - vcm_z;
            
            double nx, ny, nz;
            sample_isotropic_direction(&local_state, nx, ny, nz);
            
            double vion_cm_mag = sqrt(vion_cm_x * vion_cm_x + 
                                       vion_cm_y * vion_cm_y + 
                                       vion_cm_z * vion_cm_z);
            
            vx = nx * vion_cm_mag + vcm_x;
            vy = ny * vion_cm_mag + vcm_y;
            vz = nz * vion_cm_mag + vcm_z;
            
            vx_inout[i] = vx;
            vy_inout[i] = vy;
            vz_inout[i] = vz;
            continue;
        }
        
        int atom_offset = geometry.atom_offsets[species_idx];
        int n_atoms = geometry.atom_counts[species_idx];
        
        // Sample random molecular orientation (3 Euler angles)
        double alpha = TWO_PI * curand_uniform_double(&local_state);
        double beta = acos(2.0 * curand_uniform_double(&local_state) - 1.0);
        double gamma = TWO_PI * curand_uniform_double(&local_state);
        
        // Rotate neutral molecule geometry
        // (Store rotated positions in registers for small molecules, shared memory for large)
        
        // Build impact parameter sampling basis (perpendicular to v_rel)
        double vrel_unit_x = vrel_x / vrel_mag;
        double vrel_unit_y = vrel_y / vrel_mag;
        double vrel_unit_z = vrel_z / vrel_mag;
        
        // Tangent vectors (Gram-Schmidt)
        double t1_x, t1_y, t1_z;
        if (fabs(vrel_unit_x) < 0.9) {
            t1_x = 1.0; t1_y = 0.0; t1_z = 0.0;
        } else {
            t1_x = 0.0; t1_y = 1.0; t1_z = 0.0;
        }
        // t1 = t1 - (t1·v_rel)*v_rel
        double proj = dot(t1_x, t1_y, t1_z, vrel_unit_x, vrel_unit_y, vrel_unit_z);
        t1_x -= proj * vrel_unit_x;
        t1_y -= proj * vrel_unit_y;
        t1_z -= proj * vrel_unit_z;
        normalize(t1_x, t1_y, t1_z);
        
        double t2_x, t2_y, t2_z;
        cross(vrel_unit_x, vrel_unit_y, vrel_unit_z, t1_x, t1_y, t1_z, t2_x, t2_y, t2_z);
        
        // Estimate max impact parameter from CCS
        double bmax = sqrt(sigma / PI);
        
        // Sample impact parameter (uniform disk)
        bool collision_detected = false;
        double contact_x, contact_y, contact_z;
        
        for (int attempt = 0; attempt < MAX_COLLISION_ATTEMPTS && !collision_detected; ++attempt) {
            // Random point in disk of radius bmax
            double r = bmax * sqrt(curand_uniform_double(&local_state));
            double theta = TWO_PI * curand_uniform_double(&local_state);
            double b_x = r * cos(theta);
            double b_y = r * sin(theta);
            
            // Impact point = b_x*t1 + b_y*t2
            double impact_x = b_x * t1_x + b_y * t2_x;
            double impact_y = b_x * t1_y + b_y * t2_y;
            double impact_z = b_x * t1_z + b_y * t2_z;
            
            // Ray-trace through rotated atoms
            for (int j = 0; j < n_atoms; ++j) {
                int idx = atom_offset + j;
                double ax = geometry.atom_x[idx];
                double ay = geometry.atom_y[idx];
                double az = geometry.atom_z[idx];
                double ar = geometry.atom_radii[idx];
                
                // Rotate atom position
                double ax_rot, ay_rot, az_rot;
                rotate_euler_zyz(ax, ay, az, alpha, beta, gamma, ax_rot, ay_rot, az_rot);
                
                // Check collision: distance from ion trajectory to atom center < sum of radii
                double dx = ax_rot - impact_x;
                double dy = ay_rot - impact_y;
                double dz = az_rot - impact_z;
                
                double dist_sq = dx * dx + dy * dy + dz * dz;
                double sum_radii = ar + env.neutral_radius_m;  // Assume ion radius ~ neutral radius
                
                if (dist_sq < sum_radii * sum_radii) {
                    collision_detected = true;
                    contact_x = ax_rot;
                    contact_y = ay_rot;
                    contact_z = az_rot;
                    break;
                }
            }
            
            // Expand bmax if no collision found (adaptive sampling)
            if (!collision_detected && attempt % 64 == 63) {
                bmax *= BMAX_EXPANSION_FACTOR;
            }
        }
        
        // If collision detected, perform specular reflection
        if (collision_detected) {
            // Collision normal: from contact point to ion trajectory
            double normal_x = impact_x - contact_x;  // Simplified: should be from atom center
            double normal_y = impact_y - contact_y;
            double normal_z = impact_z - contact_z;
            normalize(normal_x, normal_y, normal_z);
            
            // Transform to COM frame
            double m_neutral = env.neutral_mass_kg;
            double m_total = m_ion + m_neutral;
            double vcm_x = (m_ion * vx + m_neutral * vn_x) / m_total;
            double vcm_y = (m_ion * vy + m_neutral * vn_y) / m_total;
            double vcm_z = (m_ion * vz + m_neutral * vn_z) / m_total;
            
            double vion_cm_x = vx - vcm_x;
            double vion_cm_y = vy - vcm_y;
            double vion_cm_z = vz - vcm_z;
            
            // Specular reflection: v' = v - 2(v·n)n
            double vdotn = dot(vion_cm_x, vion_cm_y, vion_cm_z, normal_x, normal_y, normal_z);
            vion_cm_x -= 2.0 * vdotn * normal_x;
            vion_cm_y -= 2.0 * vdotn * normal_y;
            vion_cm_z -= 2.0 * vdotn * normal_z;
            
            // Transform back to lab frame
            vx = vion_cm_x + vcm_x;
            vy = vion_cm_y + vcm_y;
            vz = vion_cm_z + vcm_z;
        } else {
            // Fallback to isotropic scattering if geometry collision fails
            double m_neutral = env.neutral_mass_kg;
            double m_total = m_ion + m_neutral;
            double vcm_x = (m_ion * vx + m_neutral * vn_x) / m_total;
            double vcm_y = (m_ion * vy + m_neutral * vn_y) / m_total;
            double vcm_z = (m_ion * vz + m_neutral * vn_z) / m_total;
            
            double vion_cm_x = vx - vcm_x;
            double vion_cm_y = vy - vcm_y;
            double vion_cm_z = vz - vcm_z;
            
            double nx, ny, nz;
            sample_isotropic_direction(&local_state, nx, ny, nz);
            
            double vion_cm_mag = sqrt(vion_cm_x * vion_cm_x + 
                                       vion_cm_y * vion_cm_y + 
                                       vion_cm_z * vion_cm_z);
            
            vx = nx * vion_cm_mag + vcm_x;
            vy = ny * vion_cm_mag + vcm_y;
            vz = nz * vion_cm_mag + vcm_z;
        }
        
        // Write back
        vx_inout[i] = vx;
        vy_inout[i] = vy;
        vz_inout[i] = vz;
    }
    
    curand_states[tid] = local_state;
}

// ============================================================================
// Host-side Launch Wrappers
// ============================================================================

void launch_hss_collision_batch(
    double* vx_inout,
    double* vy_inout,
    double* vz_inout,
    const double* mass,
    const double* ccs,
    const bool* active,
    curandState* curand_states,
    const EnvironmentParams_GPU& env,
    double dt,
    int n_ions,
    cudaStream_t stream
) {
    // Launch configuration
    constexpr int THREADS_PER_BLOCK = 256;
    int blocks = (n_ions + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    blocks = std::min(blocks, 2048);  // Limit max blocks for scheduler efficiency
    
    hss_collision_kernel<<<blocks, THREADS_PER_BLOCK, 0, stream>>>(
        vx_inout, vy_inout, vz_inout,
        mass, ccs, active,
        curand_states, env, dt, n_ions
    );
}

void launch_ehss_collision_batch(
    double* vx_inout,
    double* vy_inout,
    double* vz_inout,
    const double* mass,
    const double* ccs,
    const int* species_indices,
    const bool* active,
    curandState* curand_states,
    const EnvironmentParams_GPU& env,
    const GeometryData_GPU& geometry,
    double dt,
    int n_ions,
    cudaStream_t stream
) {
    // Launch configuration
    constexpr int THREADS_PER_BLOCK = 256;
    int blocks = (n_ions + THREADS_PER_BLOCK - 1) / THREADS_PER_BLOCK;
    blocks = std::min(blocks, 2048);
    
    ehss_collision_kernel<<<blocks, THREADS_PER_BLOCK, 0, stream>>>(
        vx_inout, vy_inout, vz_inout,
        mass, ccs, species_indices, active,
        curand_states, env, geometry, dt, n_ions
    );
}

} // namespace gpu
} // namespace ICARION

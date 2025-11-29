// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file collision_kernels_gpu.cuh
 * @brief GPU collision kernels (EHSS/HSS)
 * 
 * CUDA kernels for batch collision processing.
 * Implements EHSS and HSS collision models on GPU with cuRAND for stochastic sampling.
 * 
 * **Design:**
 * - Batch processing: one thread per ion
 * - Grid-stride loop for arbitrary ion counts
 * - cuRAND states pre-initialized, one per thread
 * - Collision probability evaluated per ion
 * - Only colliding ions perform full collision calculation
 * 
 * **Performance:**
 * - Expected 5-20× speedup vs CPU for N > 5000 ions
 * - Memory bandwidth: ~50 GB/s (ion state + geometry data)
 * - Compute: dominated by RNG and transcendental functions
 * 
 * **Phase 11 Implementation Plan:**
 * 1. HSS kernel (simpler, no geometry)
 * 2. EHSS kernel (geometry-resolved)
 * 3. cuRAND integration
 * 4. Batch collision handler class
 * 5. SimulationEngine integration
 */

#pragma once

#include <cstdint>  // For uint8_t

#ifdef __CUDACC__
#include <cuda_runtime.h>
#include <curand_kernel.h>
#else
// Host-side includes for header parsing
#include <cstddef>
#endif

namespace ICARION {
namespace gpu {

/**
 * @brief GPU-friendly environment parameters
 */
struct EnvironmentParams_GPU {
    double temperature_K;
    double pressure_Pa;
    double neutral_mass_kg;
    double neutral_radius_m;
    double gas_velocity_x;
    double gas_velocity_y;
    double gas_velocity_z;
    double mean_free_path_m;
    double thermal_velocity_m_s;
};

/**
 * @brief GPU geometry data for EHSS
 */
struct GeometryData_GPU {
    double* atom_x = nullptr;
    double* atom_y = nullptr;
    double* atom_z = nullptr;
    double* atom_radii = nullptr;
    int* atom_counts = nullptr;
    int* atom_offsets = nullptr;
    int num_species = 0;
};

}  // namespace gpu
}  // namespace ICARION

/**
 * @brief Initialize cuRAND states for collision kernels
 * 
 * Call once before collision processing. Each thread gets independent RNG state.
 * 
 * @param states Output: cuRAND states (one per thread)
 * @param seed Random seed
 * @param n_threads Total number of threads (grid_size * block_size)
 */
__global__ void init_curand_states(
    curandState* states,
    unsigned long long seed,
    int n_threads
);

/**
 * @brief HSS collision kernel (batch processing)
 * 
 * Isotropic hard-sphere scattering without geometry.
 * 
 * **Algorithm:**
 * 1. Compute collision probability: P = 1 - exp(-v_rel * σ * n * dt / 4)
 * 2. Sample uniform random: if (u < P) → collision occurs
 * 3. Sample neutral velocity from Maxwell-Boltzmann distribution
 * 4. Transform to COM frame
 * 5. Sample isotropic scattering direction (uniform sphere)
 * 6. Transform back to lab frame
 * 7. Update ion velocity
 * 
 * @param vx_inout Ion x-velocity [m/s] (modified in-place)
 * @param vy_inout Ion y-velocity [m/s]
 * @param vz_inout Ion z-velocity [m/s]
 * @param mass Ion masses [kg]
 * @param ccs Ion collision cross-sections [m²]
 * @param active Ion active flags (skip inactive ions)
 * @param curand_states cuRAND states (one per thread)
 * @param env Environment parameters (constant memory)
 * @param dt Timestep [s]
 * @param n_ions Total number of ions
 */
__global__ void hss_collision_kernel(
    double* vx_inout,
    double* vy_inout,
    double* vz_inout,
    const double* mass,
    const double* ccs,
    const uint8_t* active,
    curandState* curand_states,
    const EnvironmentParams_GPU env,
    double dt,
    int n_ions
);

/**
 * @brief EHSS collision kernel (geometry-resolved)
 * 
 * Explicit hard-sphere scattering with molecular geometry.
 * 
 * **Algorithm:**
 * 1. Compute collision probability (same as HSS)
 * 2. Sample neutral velocity from Maxwell-Boltzmann
 * 3. Look up geometry for ion species
 * 4. Randomly rotate molecule (3 Euler angles)
 * 5. Sample impact parameter in plane perpendicular to v_rel
 * 6. Ray-trace through rotated atoms, find first contact
 * 7. Compute collision normal from contact point
 * 8. Perform specular reflection in COM frame
 * 9. Transform back to lab frame
 * 10. Update ion velocity
 * 
 * @param vx_inout Ion x-velocity [m/s] (modified in-place)
 * @param vy_inout Ion y-velocity [m/s]
 * @param vz_inout Ion z-velocity [m/s]
 * @param mass Ion masses [kg]
 * @param ccs Ion collision cross-sections [m²]
 * @param species_indices Ion species indices (for geometry lookup)
 * @param active Ion active flags
 * @param curand_states cuRAND states (one per thread)
 * @param env Environment parameters (constant memory)
 * @param geometry Geometry data (global memory)
 * @param dt Timestep [s]
 * @param n_ions Total number of ions
 */
__global__ void ehss_collision_kernel(
    double* vx_inout,
    double* vy_inout,
    double* vz_inout,
    const double* mass,
    const double* ccs,
    const int* species_indices,
    const uint8_t* active,
    curandState* curand_states,
    const EnvironmentParams_GPU env,
    const GeometryData_GPU geometry,
    double dt,
    int n_ions
);

/**
 * @brief Host-side wrapper for HSS collision batch
 * 
 * Launches HSS kernel with optimal grid configuration.
 * 
 * @param vx_inout Device pointer: ion x-velocities
 * @param vy_inout Device pointer: ion y-velocities
 * @param vz_inout Device pointer: ion z-velocities
 * @param mass Device pointer: ion masses [kg]
 * @param ccs Device pointer: ion CCS [m²]
 * @param active Device pointer: ion active flags
 * @param curand_states Device pointer: cuRAND states
 * @param env Environment parameters (copied to constant memory)
 * @param dt Timestep [s]
 * @param n_ions Total number of ions
 * @param stream CUDA stream for async execution
 */
void launch_hss_collision_batch(
    double* vx_inout,
    double* vy_inout,
    double* vz_inout,
    const double* mass,
    const double* ccs,
    const uint8_t* active,
    curandState* curand_states,
    const EnvironmentParams_GPU& env,
    double dt,
    int n_ions,
    cudaStream_t stream
);

/**
 * @brief Host-side wrapper for EHSS collision batch
 * 
 * Launches EHSS kernel with optimal grid configuration.
 * 
 * @param vx_inout Device pointer: ion x-velocities
 * @param vy_inout Device pointer: ion y-velocities
 * @param vz_inout Device pointer: ion z-velocities
 * @param mass Device pointer: ion masses [kg]
 * @param ccs Device pointer: ion CCS [m²]
 * @param species_indices Device pointer: ion species indices
 * @param active Device pointer: ion active flags
 * @param curand_states Device pointer: cuRAND states
 * @param env Environment parameters (copied to constant memory)
 * @param geometry Geometry data (device pointers)
 * @param dt Timestep [s]
 * @param n_ions Total number of ions
 * @param stream CUDA stream for async execution
 */
void launch_ehss_collision_batch(
    double* vx_inout,
    double* vy_inout,
    double* vz_inout,
    const double* mass,
    const double* ccs,
    const int* species_indices,
    const uint8_t* active,
    curandState* curand_states,
    const EnvironmentParams_GPU& env,
    const GeometryData_GPU& geometry,
    double dt,
    int n_ions,
    cudaStream_t stream
);

} // namespace gpu
} // namespace ICARION

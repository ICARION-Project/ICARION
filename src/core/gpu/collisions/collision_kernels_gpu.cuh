// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file collision_kernels_gpu.cuh
 * @brief GPU collision kernels (EHSS/HSS)
 * 
 * CUDA kernels for batch collision processing. Implements experimental GPU
 * variants of the HSS (isotropic) and EHSS (geometry-resolved) models using
 * cuRAND for stochastic sampling. These kernels are currently only exercised
 * via `GPUCollisionHelper`; species-to-geometry mapping is still rudimentary
 * (single-species assumption) and results have not been validated against the
 * CPU path.
 *
 * **Design:**
 * - Batch processing: one thread per ion, grid-stride loop for arbitrary N
 * - cuRAND states must be pre-initialized (one per launched thread)
 * - Collision probability evaluated per ion; non-colliding ions do minimal work
 *
 * **Caveats:**
 * - Environment fields `mean_free_path_m` and `thermal_velocity_m_s` are
 *   unused by the kernels (kept for parity with host structs).
 * - EHSS uses a simple specular reflection with neutral radius as ion radius
 *   and falls back to isotropic scattering if no geometry hit is found.
 * - Geometry upload and species index mapping must be handled by the caller;
 *   invalid indices fall back to HSS behaviour.
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

namespace icarion {
namespace gpu {

/**
 * @brief Maximum mixture components supported in GPU path
 *
 * GPU kernels keep component data in fixed-size arrays for cache friendliness.
 * Raise this limit if future instruments routinely exceed it.
 */
constexpr int MAX_GPU_GAS_COMPONENTS = 8;

/**
 * @brief GPU-friendly environment parameters
 */
struct EnvironmentParams_GPU {
    double temperature_K = 0.0;
    double pressure_Pa = 0.0;
    double neutral_mass_kg = 0.0;
    double neutral_radius_m = 0.0;
    double gas_velocity_x = 0.0;
    double gas_velocity_y = 0.0;
    double gas_velocity_z = 0.0;
    double mean_free_path_m = 0.0;
    double thermal_velocity_m_s = 0.0;
    int num_components = 0;
    double component_density_m3[MAX_GPU_GAS_COMPONENTS] = {0.0};
    double component_mass_kg[MAX_GPU_GAS_COMPONENTS] = {0.0};
    double component_radius_m[MAX_GPU_GAS_COMPONENTS] = {0.0};
    double component_cross_section_m2[MAX_GPU_GAS_COMPONENTS] = {0.0};
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
}  // namespace icarion

// CUDA kernel declarations (must be outside namespace)

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
 * 1. Compute collision probability: P = 1 - exp(-v_rel * σ * n * dt)
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
 * @param env Environment parameters (passed by value)
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
    const icarion::gpu::EnvironmentParams_GPU env,
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
 * 7. Compute collision normal from contact point (approximate)
 * 8. Perform specular reflection in COM frame (uses neutral radius as ion radius)
 * 9. Transform back to lab frame; if no geometry hit, fall back to isotropic scattering
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
 * @param env Environment parameters (passed by value)
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
    const icarion::gpu::EnvironmentParams_GPU env,
    const icarion::gpu::GeometryData_GPU geometry,
    double dt,
    int n_ions
);

// Host-side launch wrappers (back in namespace for C++ linkage)
namespace icarion {
namespace gpu {

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
} // namespace icarion

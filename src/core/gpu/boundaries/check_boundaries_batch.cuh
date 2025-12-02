// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file check_boundaries_batch.cuh
 * @brief GPU batch boundary checking for cylindrical geometry
 * 
 * Validates ion positions against cylindrical domain boundaries:
 * - Axial bounds: -ε ≤ z < length_m
 * - Radial bound: r ≤ radius_m + ε
 *
 * Parallelizes boundary checks across all active ions in batch. Only
 * cylindrical geometry is supported on GPU; Orbitrap boundaries are currently
 * evaluated on the CPU path. Domain transitions are handled by the caller via
 * the `is_last_domain` flag.
 */

#pragma once

#include <cuda_runtime.h>
#include <cstdint>

namespace ICARION {
namespace gpu {

/**
 * @brief Geometry parameters for cylindrical domain
 */
struct CylindricalGeometry {
    double length_m;     ///< Axial length [m]
    double radius_m;     ///< Radial radius [m]
    double epsilon;      ///< Floating-point tolerance [m] (1e-12 m = 1 pm)
    bool is_last_domain; ///< True if last domain (strict z < length check)
};

/**
 * @brief Check boundaries for batch of ions (GPU kernel)
 * 
 * @param pos_x X positions [m] (device array, length N)
 * @param pos_y Y positions [m] (device array, length N)
 * @param pos_z Z positions [m] (device array, length N)
 * @param active Active flags (device array, length N, modified in-place)
 * @param geom Geometry parameters (copied into registers for each launch)
 * @param N Number of ions in batch
 * 
 * For each ion i:
 * - If position is outside bounds → active[i] = false
 * - If position is inside → active[i] unchanged
 * 
 * Boundary checks (cylindrical):
 * 1. Axial lower: z ≥ -ε (ions starting exactly at z=0 may have roundoff)
 * 2. Axial upper: z < length (last domain) or z ≤ length + ε (transition domain)
 * 3. Radial: sqrt(x² + y²) ≤ radius + ε
 */
__global__ void check_boundaries_batch_kernel(
    const double* __restrict__ pos_x,
    const double* __restrict__ pos_y,
    const double* __restrict__ pos_z,
    bool* __restrict__ active,
    CylindricalGeometry geom,
    int N
);

/**
 * @brief Launch boundary checking kernel (host wrapper)
 * 
 * @param pos_x X positions [m] (device memory)
 * @param pos_y Y positions [m] (device memory)
 * @param pos_z Z positions [m] (device memory)
 * @param active Active flags (device memory, modified)
 * @param length_m Domain length [m]
 * @param radius_m Domain radius [m]
 * @param is_last_domain True if last domain
 * @param N Number of ions
 * @param stream CUDA stream (nullptr = default stream)
 * 
 * @return cudaError_t Error code (cudaSuccess if no error)
 */
cudaError_t check_boundaries_batch(
    const double* pos_x,
    const double* pos_y,
    const double* pos_z,
    bool* active,
    double length_m,
    double radius_m,
    bool is_last_domain,
    int N,
    cudaStream_t stream = nullptr
);

}  // namespace gpu
}  // namespace ICARION

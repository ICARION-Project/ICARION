// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file check_boundaries_batch.cu
 * @brief GPU batch boundary checking implementation
 * 
 * Phase 11: GPU Acceleration - Boundary Handling Module
 */

#include "check_boundaries_batch.cuh"
#include <cstdio>

namespace ICARION {
namespace gpu {

// =============================================================================
// CUDA Kernel: Boundary Checking
// =============================================================================

/**
 * @brief Check cylindrical boundaries for batch of ions
 * 
 * Thread mapping: 1 thread = 1 ion
 * Grid size: (N + 255) / 256 blocks of 256 threads
 * 
 * Boundary logic (cylindrical geometry):
 * 1. Axial lower bound: z ≥ -ε
 *    - Allows slight negative z due to floating-point roundoff
 *    - Prevents false positives for ions initialized at z=0
 * 
 * 2. Axial upper bound:
 *    - Last domain: z < length (strict, ions exit simulation)
 *    - Transition domain: z ≤ length + ε (allow transition to next domain)
 * 
 * 3. Radial bound: r ≤ radius + ε
 *    - r = sqrt(x² + y²)
 *    - Small epsilon for roundoff tolerance
 * 
 * Performance notes:
 * - sqrt() is expensive (~20 cycles), but unavoidable
 * - Early exit on axial checks reduces sqrt() calls
 * - Coalesced memory access (stride-1 on pos_x, pos_y, pos_z, active)
 */
__global__ void check_boundaries_batch_kernel(
    const double* __restrict__ pos_x,
    const double* __restrict__ pos_y,
    const double* __restrict__ pos_z,
    bool* __restrict__ active,
    CylindricalGeometry geom,
    int N
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= N) return;
    
    // Skip if already inactive
    if (!active[idx]) return;
    
    // Load position
    double x = pos_x[idx];
    double y = pos_y[idx];
    double z = pos_z[idx];
    
    // Axial lower bound: z ≥ -ε
    if (z < -geom.epsilon) {
        active[idx] = false;
        return;
    }
    
    // Axial upper bound (depends on domain type)
    bool z_valid;
    if (geom.is_last_domain) {
        z_valid = (z < geom.length_m);  // Strict check
    } else {
        z_valid = (z <= geom.length_m + geom.epsilon);  // Transition check
    }
    
    if (!z_valid) {
        active[idx] = false;
        return;
    }
    
    // Radial bound: sqrt(x² + y²) ≤ radius + ε
    double r = sqrt(x*x + y*y);
    if (r > geom.radius_m + geom.epsilon) {
        active[idx] = false;
        return;
    }
    
    // Ion is inside boundaries, active flag unchanged
}

// =============================================================================
// Host Wrapper
// =============================================================================

cudaError_t check_boundaries_batch(
    const double* pos_x,
    const double* pos_y,
    const double* pos_z,
    bool* active,
    double length_m,
    double radius_m,
    bool is_last_domain,
    int N,
    cudaStream_t stream
) {
    if (N <= 0) return cudaSuccess;
    
    // Configure geometry
    CylindricalGeometry geom;
    geom.length_m = length_m;
    geom.radius_m = radius_m;
    geom.epsilon = 1e-12;  // 1 pm tolerance (matches CPU DOMAIN_BOUNDARY_EPSILON)
    geom.is_last_domain = is_last_domain;
    
    // Configure kernel launch
    constexpr int BLOCK_SIZE = 256;
    int grid_size = (N + BLOCK_SIZE - 1) / BLOCK_SIZE;
    
    // Launch kernel
    check_boundaries_batch_kernel<<<grid_size, BLOCK_SIZE, 0, stream>>>(
        pos_x, pos_y, pos_z, active, geom, N
    );
    
    // Check for kernel launch errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        fprintf(stderr, "[GPU Boundary Check] Kernel launch failed: %s\n", 
                cudaGetErrorString(err));
        return err;
    }
    
    // Optionally synchronize (only if stream == nullptr for debugging)
    if (stream == nullptr) {
        err = cudaDeviceSynchronize();
        if (err != cudaSuccess) {
            fprintf(stderr, "[GPU Boundary Check] Kernel execution failed: %s\n",
                    cudaGetErrorString(err));
            return err;
        }
    }
    
    return cudaSuccess;
}

}  // namespace gpu
}  // namespace ICARION

// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file FieldArrayGPU_kernels.cuh
 * @brief Device-side field interpolation functions
 * 
 * This header contains __device__ functions that must be visible to
 * CUDA kernels that use field interpolation.
 */

#ifndef ICARION_FIELD_ARRAY_GPU_KERNELS_CUH
#define ICARION_FIELD_ARRAY_GPU_KERNELS_CUH

#include "FieldArrayGPU.h"
#include "core/types/Vec3.h"
#include <cuda_runtime.h>

namespace icarion {
namespace gpu {

/**
 * @brief Device function: Interpolate electric field at position
 * 
 * Uses hardware trilinear interpolation via texture memory.
 * Returns zero field if position is outside grid bounds or if no E-field is available.
 * 
 * @param fields FieldArrayGPU with E-field textures
 * @param pos Position to sample [m]
 * @return Interpolated E-field [V/m]
 */
__device__ inline Vec3 interpolate_E_field(
    const FieldArrayGPU& fields,
    const Vec3& pos
) {
    if (!fields.has_E_field) {
        return Vec3(0.0, 0.0, 0.0);
    }
    
    // Convert world position to grid coordinates
    double gx = (pos.x - fields.origin.x) * fields.inv_spacing.x;
    double gy = (pos.y - fields.origin.y) * fields.inv_spacing.y;
    double gz = (pos.z - fields.origin.z) * fields.inv_spacing.z;
    
    // Check bounds (with small margin for numerical precision)
    if (gx < -0.5 || gx > fields.nx - 0.5 ||
        gy < -0.5 || gy > fields.ny - 0.5 ||
        gz < -0.5 || gz > fields.nz - 0.5) {
        return Vec3(0.0, 0.0, 0.0);
    }
    
    // Sample from textures (hardware trilinear interpolation)
    // Note: texture coordinates are in voxel space (0.5 = center of first voxel)
    float Ex = tex3D<float>(fields.Ex_tex, gx + 0.5f, gy + 0.5f, gz + 0.5f);
    float Ey = tex3D<float>(fields.Ey_tex, gx + 0.5f, gy + 0.5f, gz + 0.5f);
    float Ez = tex3D<float>(fields.Ez_tex, gx + 0.5f, gy + 0.5f, gz + 0.5f);
    
    return Vec3(static_cast<double>(Ex), static_cast<double>(Ey), static_cast<double>(Ez));
}

/**
 * @brief Device function: Interpolate magnetic field at position
 * 
 * Uses hardware trilinear interpolation via texture memory.
 * Returns zero field if position is outside grid bounds or if no B-field is available.
 * 
 * @param fields FieldArrayGPU with B-field textures
 * @param pos Position to sample [m]
 * @return Interpolated B-field [T]
 */
__device__ inline Vec3 interpolate_B_field(
    const FieldArrayGPU& fields,
    const Vec3& pos
) {
    if (!fields.has_B_field) {
        return Vec3(0.0, 0.0, 0.0);
    }
    
    // Convert world position to grid coordinates
    double gx = (pos.x - fields.origin.x) * fields.inv_spacing.x;
    double gy = (pos.y - fields.origin.y) * fields.inv_spacing.y;
    double gz = (pos.z - fields.origin.z) * fields.inv_spacing.z;
    
    // Check bounds
    if (gx < -0.5 || gx > fields.nx - 0.5 ||
        gy < -0.5 || gy > fields.ny - 0.5 ||
        gz < -0.5 || gz > fields.nz - 0.5) {
        return Vec3(0.0, 0.0, 0.0);
    }
    
    // Sample from textures
    float Bx = tex3D<float>(fields.Bx_tex, gx + 0.5f, gy + 0.5f, gz + 0.5f);
    float By = tex3D<float>(fields.By_tex, gx + 0.5f, gy + 0.5f, gz + 0.5f);
    float Bz = tex3D<float>(fields.Bz_tex, gx + 0.5f, gy + 0.5f, gz + 0.5f);
    
    return Vec3(static_cast<double>(Bx), static_cast<double>(By), static_cast<double>(Bz));
}

} // namespace gpu
} // namespace icarion

#endif // ICARION_FIELD_ARRAY_GPU_KERNELS_CUH

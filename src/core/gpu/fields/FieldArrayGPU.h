// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_FIELD_ARRAY_GPU_H
#define ICARION_FIELD_ARRAY_GPU_H

#include "core/types/Vec3.h"
#include <cuda_runtime.h>

namespace icarion {
namespace gpu {

/**
 * @brief GPU-side field array with texture memory (single precision)
 * 
 * Stores electric and magnetic field components in CUDA texture memory for
 * hardware trilinear interpolation. Converts host doubles to float; intended
 * for the experimental GPUIntegrationHelper path.
 * 
 * Memory layout:
 * - 3D texture objects for each field component (Ex, Ey, Ez, Bx, By, Bz)
 * - Grid metadata (origin, spacing, dimensions)
 * - CUDA arrays backing the texture objects
 * 
 * @note Uses single-precision textures; precision-sensitive setups should
 *       prefer CPU interpolation.
 */
struct FieldArrayGPU {
    // Texture objects for electric field components
    cudaTextureObject_t Ex_tex{0};
    cudaTextureObject_t Ey_tex{0};
    cudaTextureObject_t Ez_tex{0};
    
    // Texture objects for magnetic field components
    cudaTextureObject_t Bx_tex{0};
    cudaTextureObject_t By_tex{0};
    cudaTextureObject_t Bz_tex{0};
    
    // Grid geometry
    Vec3 origin{0.0, 0.0, 0.0};    // Grid origin [m]
    Vec3 spacing{1.0, 1.0, 1.0};   // Grid spacing [m]
    Vec3 inv_spacing{1.0, 1.0, 1.0}; // 1/spacing for fast division
    
    // Grid dimensions
    int nx{0}, ny{0}, nz{0};
    
    // Field availability flags
    bool has_E_field{false};
    bool has_B_field{false};
    
    // Backing CUDA arrays (for cleanup)
    cudaArray_t Ex_array{nullptr};
    cudaArray_t Ey_array{nullptr};
    cudaArray_t Ez_array{nullptr};
    cudaArray_t Bx_array{nullptr};
    cudaArray_t By_array{nullptr};
    cudaArray_t Bz_array{nullptr};
};

/**
 * @brief Upload electric field to GPU texture memory
 * 
 * @param Ex CPU array of Ex component [V/m] (size: nx*ny*nz)
 * @param Ey CPU array of Ey component [V/m] (size: nx*ny*nz)
 * @param Ez CPU array of Ez component [V/m] (size: nx*ny*nz)
 * @param nx Number of grid points in x
 * @param ny Number of grid points in y
 * @param nz Number of grid points in z
 * @param origin Grid origin [m]
 * @param spacing Grid spacing [m]
 * @param fields Output: FieldArrayGPU structure with initialized textures
 * 
 * @note The caller must call free_field_array_gpu() to release resources.
 */
void upload_E_field(
    const double* Ex, const double* Ey, const double* Ez,
    int nx, int ny, int nz,
    const Vec3& origin,
    const Vec3& spacing,
    FieldArrayGPU& fields
);

/**
 * @brief Upload magnetic field to GPU texture memory
 * 
 * @param Bx CPU array of Bx component [T] (size: nx*ny*nz)
 * @param By CPU array of By component [T] (size: nx*ny*nz)
 * @param Bz CPU array of Bz component [T] (size: nx*ny*nz)
 * @param nx Number of grid points in x
 * @param ny Number of grid points in y
 * @param nz Number of grid points in z
 * @param origin Grid origin [m]
 * @param spacing Grid spacing [m]
 * @param fields Output: FieldArrayGPU structure with initialized textures
 * 
 * @note The caller must call free_field_array_gpu() to release resources.
 */
void upload_B_field(
    const double* Bx, const double* By, const double* Bz,
    int nx, int ny, int nz,
    const Vec3& origin,
    const Vec3& spacing,
    FieldArrayGPU& fields
);

/**
 * @brief Free GPU field array resources
 * 
 * Destroys texture objects and frees CUDA arrays.
 * 
 * @param fields FieldArrayGPU to free
 */
void free_field_array_gpu(FieldArrayGPU& fields);

} // namespace gpu
} // namespace icarion

// Include device-side interpolation functions for CUDA kernels
// These must be inline __device__ functions visible to all kernels
#ifdef __CUDACC__
#include "FieldArrayGPU_kernels.cuh"
#endif

#endif // ICARION_FIELD_ARRAY_GPU_H

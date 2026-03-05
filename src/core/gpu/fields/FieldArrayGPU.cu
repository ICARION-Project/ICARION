// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file FieldArrayGPU.cu
 * @brief Implementation of GPU field interpolation (single-precision textures)
 * @details Converts host double grids to float textures; experimental path
 *          used only by GPUIntegrationHelper.
 */

#include "FieldArrayGPU.h"
#include <cuda_runtime.h>
#include <stdexcept>
#include <cmath>

namespace icarion {
namespace gpu {

// ============================================================================
// Helper: Create 3D texture from CPU data
// ============================================================================

static void create_3d_texture(
    const double* cpu_data,
    int nx, int ny, int nz,
    cudaArray_t& cuda_array,
    cudaTextureObject_t& tex_obj
) {
    // Create 3D CUDA array
    cudaChannelFormatDesc channel_desc = cudaCreateChannelDesc<float>();
    cudaExtent extent = make_cudaExtent(nx, ny, nz);
    
    cudaError_t err = cudaMalloc3DArray(&cuda_array, &channel_desc, extent);
    if (err != cudaSuccess) {
        throw std::runtime_error("Failed to allocate CUDA 3D array: " + 
                                std::string(cudaGetErrorString(err)));
    }
    
    // Convert double to float and copy to device
    size_t total_size = nx * ny * nz;
    float* float_data = new float[total_size];
    for (size_t i = 0; i < total_size; ++i) {
        float_data[i] = static_cast<float>(cpu_data[i]);
    }
    
    // Copy data to 3D array
    cudaMemcpy3DParms copy_params{};
    copy_params.srcPtr = make_cudaPitchedPtr(
        (void*)float_data, nx * sizeof(float), nx, ny
    );
    copy_params.dstArray = cuda_array;
    copy_params.extent = extent;
    copy_params.kind = cudaMemcpyHostToDevice;
    
    err = cudaMemcpy3D(&copy_params);
    delete[] float_data;  // Free temporary buffer
    if (err != cudaSuccess) {
        cudaFreeArray(cuda_array);
        throw std::runtime_error("Failed to copy data to CUDA array: " + 
                                std::string(cudaGetErrorString(err)));
    }
    
    // Create texture object
    cudaResourceDesc res_desc;
    memset(&res_desc, 0, sizeof(res_desc));
    res_desc.resType = cudaResourceTypeArray;
    res_desc.res.array.array = cuda_array;
    
    cudaTextureDesc tex_desc;
    memset(&tex_desc, 0, sizeof(tex_desc));
    tex_desc.addressMode[0] = cudaAddressModeClamp;
    tex_desc.addressMode[1] = cudaAddressModeClamp;
    tex_desc.addressMode[2] = cudaAddressModeClamp;
    tex_desc.filterMode = cudaFilterModeLinear;  // Hardware trilinear interpolation
    tex_desc.readMode = cudaReadModeElementType;
    tex_desc.normalizedCoords = 0;  // Use unnormalized coordinates
    
    err = cudaCreateTextureObject(&tex_obj, &res_desc, &tex_desc, nullptr);
    if (err != cudaSuccess) {
        cudaFreeArray(cuda_array);
        throw std::runtime_error("Failed to create texture object: " + 
                                std::string(cudaGetErrorString(err)));
    }
}

// ============================================================================
// Public API: Upload fields
// ============================================================================

void upload_E_field(
    const double* Ex, const double* Ey, const double* Ez,
    int nx, int ny, int nz,
    const Vec3& origin,
    const Vec3& spacing,
    FieldArrayGPU& fields
) {
    // Store grid geometry
    fields.nx = nx;
    fields.ny = ny;
    fields.nz = nz;
    fields.origin = origin;
    fields.spacing = spacing;
    fields.inv_spacing = Vec3(1.0 / spacing.x, 1.0 / spacing.y, 1.0 / spacing.z);
    
    // Create textures for each component
    create_3d_texture(Ex, nx, ny, nz, fields.Ex_array, fields.Ex_tex);
    create_3d_texture(Ey, nx, ny, nz, fields.Ey_array, fields.Ey_tex);
    create_3d_texture(Ez, nx, ny, nz, fields.Ez_array, fields.Ez_tex);
    
    fields.has_E_field = true;
}

void upload_B_field(
    const double* Bx, const double* By, const double* Bz,
    int nx, int ny, int nz,
    const Vec3& origin,
    const Vec3& spacing,
    FieldArrayGPU& fields
) {
    // Store grid geometry (if not already set by E-field)
    if (fields.nx == 0) {
        fields.nx = nx;
        fields.ny = ny;
        fields.nz = nz;
        fields.origin = origin;
        fields.spacing = spacing;
        fields.inv_spacing = Vec3(1.0 / spacing.x, 1.0 / spacing.y, 1.0 / spacing.z);
    }
    
    // Create textures for each component
    create_3d_texture(Bx, nx, ny, nz, fields.Bx_array, fields.Bx_tex);
    create_3d_texture(By, nx, ny, nz, fields.By_array, fields.By_tex);
    create_3d_texture(Bz, nx, ny, nz, fields.Bz_array, fields.Bz_tex);
    
    fields.has_B_field = true;
}

void free_field_array_gpu(FieldArrayGPU& fields) {
    // Destroy texture objects
    if (fields.Ex_tex) cudaDestroyTextureObject(fields.Ex_tex);
    if (fields.Ey_tex) cudaDestroyTextureObject(fields.Ey_tex);
    if (fields.Ez_tex) cudaDestroyTextureObject(fields.Ez_tex);
    if (fields.Bx_tex) cudaDestroyTextureObject(fields.Bx_tex);
    if (fields.By_tex) cudaDestroyTextureObject(fields.By_tex);
    if (fields.Bz_tex) cudaDestroyTextureObject(fields.Bz_tex);
    
    // Free CUDA arrays
    if (fields.Ex_array) cudaFreeArray(fields.Ex_array);
    if (fields.Ey_array) cudaFreeArray(fields.Ey_array);
    if (fields.Ez_array) cudaFreeArray(fields.Ez_array);
    if (fields.Bx_array) cudaFreeArray(fields.Bx_array);
    if (fields.By_array) cudaFreeArray(fields.By_array);
    if (fields.Bz_array) cudaFreeArray(fields.Bz_array);
    
    // Reset structure
    fields = FieldArrayGPU{};
}

// Device functions are now in FieldArrayGPU_kernels.cuh

} // namespace gpu
} // namespace icarion

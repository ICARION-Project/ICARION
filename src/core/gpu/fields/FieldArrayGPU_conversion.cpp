// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file FieldArrayGPU_conversion.cpp
 * @brief Implementation of CPU-GPU field array conversion
 */

#include "FieldArrayGPU_conversion.h"
#include <stdexcept>
#include <cmath>

namespace icarion {
namespace gpu {

bool is_field_valid_for_gpu(const FieldArray& cpu_field) {
    // Check basic validity
    if (!cpu_field.is_valid()) {
        return false;
    }
    
    // Check that we have field data
    size_t expected_size = cpu_field.nx * cpu_field.ny * cpu_field.nz;
    if (cpu_field.Ex.size() != expected_size ||
        cpu_field.Ey.size() != expected_size ||
        cpu_field.Ez.size() != expected_size) {
        return false;
    }
    
    // Check grid spacing consistency
    if (cpu_field.nx < 2 || cpu_field.ny < 2 || cpu_field.nz < 2) {
        return false;  // Need at least 2 points per dimension for interpolation
    }
    
    return true;
}

void upload_field_array_to_gpu(
    const FieldArray& cpu_field,
    FieldArrayGPU& gpu_field
) {
    // Validate input
    if (!is_field_valid_for_gpu(cpu_field)) {
        throw std::runtime_error(
            "FieldArray is invalid for GPU upload. "
            "Check dimensions and field data consistency."
        );
    }
    
    // Extract grid parameters
    int nx = static_cast<int>(cpu_field.nx);
    int ny = static_cast<int>(cpu_field.ny);
    int nz = static_cast<int>(cpu_field.nz);
    
    // Compute grid origin and spacing
    // Origin is at first grid point
    Vec3 origin(
        cpu_field.xs.front(),
        cpu_field.ys.front(),
        cpu_field.zs.front()
    );
    
    // Spacing between grid points (assume uniform spacing)
    Vec3 spacing;
    if (nx > 1) {
        spacing.x = (cpu_field.xs.back() - cpu_field.xs.front()) / (nx - 1);
    } else {
        spacing.x = 1.0;  // Fallback (shouldn't happen due to validation)
    }
    
    if (ny > 1) {
        spacing.y = (cpu_field.ys.back() - cpu_field.ys.front()) / (ny - 1);
    } else {
        spacing.y = 1.0;
    }
    
    if (nz > 1) {
        spacing.z = (cpu_field.zs.back() - cpu_field.zs.front()) / (nz - 1);
    } else {
        spacing.z = 1.0;
    }
    
    // Upload electric field to GPU
    upload_E_field(
        cpu_field.Ex.data(),
        cpu_field.Ey.data(),
        cpu_field.Ez.data(),
        nx, ny, nz,
        origin,
        spacing,
        gpu_field
    );
    
    // Note: Magnetic field not uploaded (FieldArray doesn't store Bx, By, Bz).
    // Only the E-field path is exercised by the experimental GPUIntegrationHelper.
}

} // namespace gpu
} // namespace icarion

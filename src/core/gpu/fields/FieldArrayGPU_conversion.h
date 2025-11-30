// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file FieldArrayGPU_conversion.h
 * @brief Conversion between CPU FieldArray and GPU FieldArrayGPU
 * 
 * Provides bridge between ICARION's CPU field representation (FieldArray)
 * and GPU texture-based field storage (FieldArrayGPU).
 */

#ifndef ICARION_FIELD_ARRAY_GPU_CONVERSION_H
#define ICARION_FIELD_ARRAY_GPU_CONVERSION_H

#include "core/gpu/fields/FieldArrayGPU.h"
#include "core/io/fieldArrayLoader.h"

namespace icarion {
namespace gpu {

/**
 * @brief Upload CPU FieldArray to GPU as FieldArrayGPU
 * 
 * Converts ICARION's standard FieldArray format (vectors of Ex, Ey, Ez)
 * to GPU texture memory format for hardware-accelerated interpolation.
 * 
 * @param cpu_field CPU field array (from load_field_array() or similar)
 * @param gpu_field Output GPU field array with initialized textures
 * 
 * @throws std::runtime_error if field is invalid or GPU upload fails
 * 
 * Usage:
 * @code
 * FieldArray cpu_field = load_field_array("fields.h5");
 * FieldArrayGPU gpu_field;
 * upload_field_array_to_gpu(cpu_field, gpu_field);
 * // ... use gpu_field in kernels ...
 * free_field_array_gpu(gpu_field);
 * @endcode
 * 
 * @note Only electric field (Ex, Ey, Ez) is uploaded. Magnetic field
 *       support can be added if FieldArray is extended to include Bx, By, Bz.
 */
void upload_field_array_to_gpu(
    const FieldArray& cpu_field,
    FieldArrayGPU& gpu_field
);

/**
 * @brief Check if a FieldArray can be uploaded to GPU
 * 
 * Validates that the field has:
 * - Valid dimensions (nx, ny, nz > 0)
 * - Consistent array sizes
 * - Non-empty field data
 * 
 * @param cpu_field Field array to check
 * @return true if field is valid for GPU upload
 */
bool is_field_valid_for_gpu(const FieldArray& cpu_field);

} // namespace gpu
} // namespace icarion

#endif // ICARION_FIELD_ARRAY_GPU_CONVERSION_H

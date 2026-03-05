// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_INTEGRATE_RK4_BATCH_CUH
#define ICARION_INTEGRATE_RK4_BATCH_CUH

#include "core/types/Vec3.h"
#include "core/types/gpu/IonState_GPU.h"
#include "core/gpu/damping/DeviceDamping.h"
#include <cuda_runtime.h>

namespace icarion {
namespace gpu {

// Forward declaration
struct FieldArrayGPU;

/**
 * @brief Batch RK4 integration on GPU with constant fields (experimental)
 * 
 * Integrates all ions from t to t+dt using 4th-order Runge-Kutta. Used only
 * by GPUIntegrationHelper (and GPUIntegrationStrategy when eligible).
 * 
 * @param ions_in Input ion states at time t
 * @param ions_out Output ion states at time t+dt (must be pre-allocated)
 * @param E_field Electric field vector [V/m] (constant)
 * @param B_field Magnetic field vector [T] (constant)
 * @param dt Timestep [s]
 * @param stream CUDA stream for async execution (default: 0)
 * 
 * @note Inactive ions are copied unchanged to output.
 */
void integrate_rk4_batch(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const Vec3& E_field,
    const Vec3& B_field,
    double dt,
    const DeviceDamping& damping,
    cudaStream_t stream = 0
);

/**
 * @brief Batch RK4 integration on GPU with field interpolation (experimental)
 * 
 * Integrates all ions from t to t+dt using 4th-order Runge-Kutta with
 * position-dependent field evaluation from texture memory. Field textures are
 * single-precision; CPU parity not guaranteed.
 * 
 * @param ions_in Input ion states at time t
 * @param ions_out Output ion states at time t+dt (must be pre-allocated)
 * @param field_array GPU field array with texture objects
 * @param dt Timestep [s]
 * @param stream CUDA stream for async execution (default: 0)
 * 
 * @note Field interpolation uses hardware trilinear interpolation (fast!)
 * @note Inactive ions are copied unchanged to output.
 */
void integrate_rk4_batch_with_fields(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const FieldArrayGPU& field_array,
    double dt,
    const DeviceDamping& damping,
    cudaStream_t stream = 0
);

} // namespace gpu
} // namespace icarion

#endif // ICARION_INTEGRATE_RK4_BATCH_CUH

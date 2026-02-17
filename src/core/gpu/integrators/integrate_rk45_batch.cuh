// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#ifndef ICARION_INTEGRATE_RK45_BATCH_CUH
#define ICARION_INTEGRATE_RK45_BATCH_CUH

#include "core/types/Vec3.h"
#include "core/types/gpu/IonState_GPU.h"
#include <cuda_runtime.h>

namespace icarion {
namespace gpu {

// Forward declaration
struct FieldArrayGPU;

/**
 * @brief RK45 (Dormand-Prince) adaptive timestep parameters
 */
struct RK45Params {
    double atol = 1e-12;              ///< Absolute tolerance [m, m/s]
    double rtol = 1e-6;               ///< Relative tolerance
    double safety_factor = 0.9;       ///< Safety factor for step size
    double min_step_factor = 0.2;     ///< Minimum dt reduction per step
    double max_step_factor = 5.0;     ///< Maximum dt increase per step
    int max_substeps = 1000;          ///< Max substeps per main step
};

/**
 * @brief Batch RK45 integration on GPU with constant fields (experimental)
 * 
 * Integrates all ions from t to t+dt using adaptive 4th/5th-order
 * Dormand-Prince method with local error control. Used by GPUIntegrationHelper
 * (and GPUIntegrationStrategy when eligible). Experimental and not fully
 * validated against the CPU RK45.
 * 
 * Each ion adapts its substep independently based on local error.
 * 
 * @param ions_in Input ion states at time t
 * @param ions_out Output ion states at time t+dt (must be pre-allocated)
 * @param E_field Electric field vector [V/m] (constant)
 * @param B_field Magnetic field vector [T] (constant)
 * @param dt Maximum timestep [s]
 * @param params RK45 adaptive parameters
 * @param stream CUDA stream for async execution (default: 0)
 * 
 * @note Each ion may take different number of substeps
 * @note Inactive ions are copied unchanged to output
 */
void integrate_rk45_batch(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const Vec3& E_field,
    const Vec3& B_field,
    double dt,
    const RK45Params& params,
    cudaStream_t stream = 0
);

/**
 * @brief Batch RK45 integration on GPU with field interpolation (experimental)
 * 
 * Integrates all ions from t to t+dt using adaptive Dormand-Prince with
 * position-dependent field evaluation from texture memory. Field textures are
 * single-precision; CPU parity not guaranteed.
 * 
 * @param ions_in Input ion states at time t
 * @param ions_out Output ion states at time t+dt (must be pre-allocated)
 * @param field_array GPU field array with texture objects
 * @param dt Maximum timestep [s]
 * @param params RK45 adaptive parameters
 * @param stream CUDA stream for async execution (default: 0)
 * 
 * @note Field interpolation uses hardware trilinear interpolation
 * @note Each ion adapts substep size independently
 * @note Inactive ions are copied unchanged to output
 */
void integrate_rk45_batch_with_fields(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const FieldArrayGPU& field_array,
    double dt,
    const RK45Params& params,
    cudaStream_t stream = 0
);

} // namespace gpu
} // namespace icarion

#endif // ICARION_INTEGRATE_RK45_BATCH_CUH

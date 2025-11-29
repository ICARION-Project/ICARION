#ifndef ICARION_INTEGRATE_RK4_BATCH_CUH
#define ICARION_INTEGRATE_RK4_BATCH_CUH

#include "core/types/Vec3.h"
#include "utils/IonState_GPU.h"
#include <cuda_runtime.h>

namespace icarion {
namespace gpu {

/**
 * @brief Batch RK4 integration on GPU
 * 
 * Integrates all ions from t to t+dt using 4th-order Runge-Kutta.
 * Uses grid-stride loop pattern for optimal performance across all N.
 * 
 * @param ions_in Input ion states at time t
 * @param ions_out Output ion states at time t+dt (must be pre-allocated)
 * @param E_field Electric field vector [V/m] (constant for now)
 * @param B_field Magnetic field vector [T] (constant for now)
 * @param dt Timestep [s]
 * @param stream CUDA stream for async execution (default: 0)
 * 
 * @note Currently uses constant fields. Phase 3 will add field interpolation.
 * @note Inactive ions are copied unchanged to output.
 */
void integrate_rk4_batch(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const Vec3& E_field,
    const Vec3& B_field,
    double dt,
    cudaStream_t stream = 0
);

} // namespace gpu
} // namespace icarion

#endif // ICARION_INTEGRATE_RK4_BATCH_CUH

#ifndef ICARION_INTEGRATE_BORIS_BATCH_CUH
#define ICARION_INTEGRATE_BORIS_BATCH_CUH

#include "core/types/Vec3.h"
#include "utils/IonState_GPU.h"
#include <cuda_runtime.h>

namespace icarion {
namespace gpu {

// Forward declaration
struct FieldArrayGPU;

/**
 * @brief Batch Boris pusher integration on GPU with constant fields
 * 
 * Integrates all ions from t to t+dt using symplectic Boris algorithm.
 * Optimal for magnetic field-dominated systems (ICR, Orbitrap, Penning).
 * 
 * Algorithm:
 *   1. Half-step electric kick: v^- = v^n + (q/m)*E*dt/2
 *   2. Magnetic rotation: v^+ = rotate(v^-, (q/m)*B*dt)
 *   3. Half-step electric kick: v^(n+1) = v^+ + (q/m)*E*dt/2
 *   4. Position update: x^(n+1) = x^n + v^(n+1)*dt
 * 
 * @param ions_in Input ion states at time t
 * @param ions_out Output ion states at time t+dt (must be pre-allocated)
 * @param E_field Electric field vector [V/m] (constant)
 * @param B_field Magnetic field vector [T] (constant)
 * @param dt Timestep [s]
 * @param stream CUDA stream for async execution (default: 0)
 * 
 * @note Symplectic → conserves energy in pure B-field
 * @note Inactive ions are copied unchanged to output
 */
void integrate_boris_batch(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const Vec3& E_field,
    const Vec3& B_field,
    double dt,
    cudaStream_t stream = 0
);

/**
 * @brief Batch Boris pusher integration on GPU with field interpolation
 * 
 * Integrates all ions from t to t+dt using Boris algorithm with
 * position-dependent field evaluation from texture memory.
 * 
 * Fields are evaluated at:
 *   - E(x^n) for both electric kicks
 *   - B(x^n) for magnetic rotation
 * 
 * @param ions_in Input ion states at time t
 * @param ions_out Output ion states at time t+dt (must be pre-allocated)
 * @param field_array GPU field array with texture objects
 * @param dt Timestep [s]
 * @param stream CUDA stream for async execution (default: 0)
 * 
 * @note Field interpolation uses hardware trilinear interpolation
 * @note Symplectic properties preserved for smooth fields
 * @note Inactive ions are copied unchanged to output
 */
void integrate_boris_batch_with_fields(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const FieldArrayGPU& field_array,
    double dt,
    cudaStream_t stream = 0
);

} // namespace gpu
} // namespace icarion

#endif // ICARION_INTEGRATE_BORIS_BATCH_CUH

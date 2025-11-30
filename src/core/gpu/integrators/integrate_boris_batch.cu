// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file integrate_boris_batch.cu
 * @brief GPU batch integration using Boris pusher algorithm
 * 
 * Symplectic method for charged particle motion in electromagnetic fields.
 * Optimal for magnetic field-dominated systems (ICR, Orbitrap, Penning trap).
 * 
 * Algorithm:
 *   1. Half-step electric kick: v^- = v^n + (q/m)*E*dt/2
 *   2. Magnetic rotation: v^+ = rotate(v^-, (q/m)*B*dt)
 *   3. Half-step electric kick: v^(n+1) = v^+ + (q/m)*E*dt/2
 *   4. Position update: x^(n+1) = x^n + v^(n+1)*dt
 */

#include "integrate_boris_batch.cuh"
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/fields/FieldArrayGPU.h"
#include "core/gpu/fields/FieldArrayGPU_kernels.cuh"
#include <cuda_runtime.h>
#include <cmath>

namespace icarion {
namespace gpu {

/**
 * @brief Boris magnetic rotation
 * 
 * Rotates velocity vector by angle θ = |ω_B| * dt where ω_B = (q/m)*B
 * 
 * Algorithm (Boris 1970):
 *   t = (q/m) * B * (dt/2)
 *   s = 2*t / (1 + t·t)
 *   v' = v + v × t
 *   v_out = v + v' × s
 * 
 * This is exact for constant B field and second-order accurate for varying B.
 */
__device__ Vec3 boris_rotation(
    const Vec3& v_in,
    const Vec3& B_field,
    double charge,
    double mass,
    double dt
) {
    // t = (q/m) * B * (dt/2)
    double qm_half = 0.5 * charge / mass;
    Vec3 t = B_field * (qm_half * dt);
    
    // t² = t·t
    double t_squared = t.x * t.x + t.y * t.y + t.z * t.z;
    
    // s = 2*t / (1 + t²)
    Vec3 s = t * (2.0 / (1.0 + t_squared));
    
    // v' = v + v × t
    Vec3 v_cross_t = {
        v_in.y * t.z - v_in.z * t.y,
        v_in.z * t.x - v_in.x * t.z,
        v_in.x * t.y - v_in.y * t.x
    };
    Vec3 v_prime = v_in + v_cross_t;
    
    // v_out = v + v' × s
    Vec3 vprime_cross_s = {
        v_prime.y * s.z - v_prime.z * s.y,
        v_prime.z * s.x - v_prime.x * s.z,
        v_prime.x * s.y - v_prime.y * s.x
    };
    Vec3 v_out = v_in + vprime_cross_s;
    
    return v_out;
}

/**
 * @brief Boris pusher integration kernel
 * 
 * Each thread processes multiple ions using grid-stride pattern.
 */
__global__ void integrate_boris_batch_kernel(
    // Input state (t)
    const double* __restrict__ x_in,
    const double* __restrict__ y_in,
    const double* __restrict__ z_in,
    const double* __restrict__ vx_in,
    const double* __restrict__ vy_in,
    const double* __restrict__ vz_in,
    const double* __restrict__ mass,
    const double* __restrict__ charge,
    const bool* __restrict__ active_in,
    
    // Output state (t + dt)
    double* __restrict__ x_out,
    double* __restrict__ y_out,
    double* __restrict__ z_out,
    double* __restrict__ vx_out,
    double* __restrict__ vy_out,
    double* __restrict__ vz_out,
    bool* __restrict__ active_out,
    
    // Fields
    const FieldArrayGPU* field_array,
    Vec3 E_const,
    Vec3 B_const,
    
    // Integration parameters
    double dt,
    int N
) {
    // Grid-stride loop
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += gridDim.x * blockDim.x) {
        active_out[i] = active_in[i];
        
        if (!active_in[i]) {
            // Copy unchanged
            x_out[i] = x_in[i];
            y_out[i] = y_in[i];
            z_out[i] = z_in[i];
            vx_out[i] = vx_in[i];
            vy_out[i] = vy_in[i];
            vz_out[i] = vz_in[i];
            continue;
        }
        
        // Load ion state
        Vec3 pos = {x_in[i], y_in[i], z_in[i]};
        Vec3 vel = {vx_in[i], vy_in[i], vz_in[i]};
        double m = mass[i];
        double q = charge[i];
        
        // Evaluate fields at current position
        Vec3 E_field, B_field;
        if (field_array != nullptr) {
            E_field = interpolate_E_field(*field_array, pos);
            B_field = interpolate_B_field(*field_array, pos);
        } else {
            E_field = E_const;
            B_field = B_const;
        }
        
        // Boris algorithm:
        
        // 1. Half-step electric kick: v^- = v^n + (q/m)*E*dt/2
        double qm = q / m;
        Vec3 v_minus = vel + E_field * (qm * dt * 0.5);
        
        // 2. Magnetic rotation: v^+ = rotate(v^-, (q/m)*B*dt)
        Vec3 v_plus = boris_rotation(v_minus, B_field, q, m, dt);
        
        // 3. Half-step electric kick: v^(n+1) = v^+ + (q/m)*E*dt/2
        Vec3 vel_new = v_plus + E_field * (qm * dt * 0.5);
        
        // 4. Position update: x^(n+1) = x^n + v^(n+1)*dt
        Vec3 pos_new = pos + vel_new * dt;
        
        // Store result
        x_out[i] = pos_new.x;
        y_out[i] = pos_new.y;
        z_out[i] = pos_new.z;
        vx_out[i] = vel_new.x;
        vy_out[i] = vel_new.y;
        vz_out[i] = vel_new.z;
    }
}

// ============================================================================
// Host-side API
// ============================================================================

void integrate_boris_batch(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const Vec3& E_field,
    const Vec3& B_field,
    double dt,
    cudaStream_t stream
) {
    int N = ions_in.count;
    if (N == 0) return;
    
    // Launch configuration (256 threads per block)
    int threads = 256;
    int blocks = (N + threads - 1) / threads;
    
    integrate_boris_batch_kernel<<<blocks, threads, 0, stream>>>(
        ions_in.x, ions_in.y, ions_in.z,
        ions_in.vx, ions_in.vy, ions_in.vz,
        ions_in.mass, ions_in.charge, ions_in.active,
        ions_out.x, ions_out.y, ions_out.z,
        ions_out.vx, ions_out.vy, ions_out.vz,
        ions_out.active,
        nullptr,  // No field array
        E_field, B_field,
        dt, N
    );
}

void integrate_boris_batch_with_fields(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const FieldArrayGPU& field_array,
    double dt,
    cudaStream_t stream
) {
    int N = ions_in.count;
    if (N == 0) return;
    
    int threads = 256;
    int blocks = (N + threads - 1) / threads;
    
    integrate_boris_batch_kernel<<<blocks, threads, 0, stream>>>(
        ions_in.x, ions_in.y, ions_in.z,
        ions_in.vx, ions_in.vy, ions_in.vz,
        ions_in.mass, ions_in.charge, ions_in.active,
        ions_out.x, ions_out.y, ions_out.z,
        ions_out.vx, ions_out.vy, ions_out.vz,
        ions_out.active,
        &field_array,
        Vec3{0, 0, 0}, Vec3{0, 0, 0},  // Not used with field array
        dt, N
    );
}

} // namespace gpu
} // namespace icarion

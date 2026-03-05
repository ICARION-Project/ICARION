// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file integrate_rk45_batch.cu
 * @brief GPU batch integration using adaptive RK45 (Dormand-Prince)
 * 
 * Implements embedded 4th/5th order Runge-Kutta with automatic
 * step size control based on local error estimation.
 * 
 * Each ion adapts its substep size independently → optimal efficiency.
 */

#include "integrate_rk45_batch.cuh"
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/fields/FieldArrayGPU.h"
#include "core/gpu/fields/FieldArrayGPU_kernels.cuh"
#include <cuda_runtime.h>
#include <cmath>

namespace icarion {
namespace gpu {

/**
 * @brief Compute acceleration for a single ion (shared with RK4)
 */
__device__ Vec3 compute_acceleration_rk45(
    const Vec3& pos,
    const Vec3& vel,
    double mass,
    double charge,
    const FieldArrayGPU* fields,
    const Vec3& E_const,
    const Vec3& B_const
) {
    Vec3 E_field, B_field;
    if (fields != nullptr) {
        E_field = interpolate_E_field(*fields, pos);
        B_field = interpolate_B_field(*fields, pos);
    } else {
        E_field = E_const;
        B_field = B_const;
    }
    
    // F = q*E + q*(v × B)
    Vec3 F_electric = E_field * charge;
    Vec3 v_cross_B = {
        vel.y * B_field.z - vel.z * B_field.y,
        vel.z * B_field.x - vel.x * B_field.z,
        vel.x * B_field.y - vel.y * B_field.x
    };
    Vec3 F_magnetic = v_cross_B * charge;
    return (F_electric + F_magnetic) / mass;
}

/**
 * @brief Dormand-Prince RK45 coefficients (embedded 4/5 pair)
 * 
 * Butcher tableau:
 *    0   |
 *   1/5  | 1/5
 *   3/10 | 3/40      9/40
 *   4/5  | 44/45    -56/15     32/9
 *   8/9  | 19372/6561  -25360/2187  64448/6561  -212/729
 *    1   | 9017/3168   -355/33   46732/5247   49/176  -5103/18656
 *    1   | 35/384       0        500/1113     125/192  -2187/6784   11/84  (4th order)
 *        | 5179/57600  0        7571/16695   393/640  -92097/339200  187/2100  1/40  (5th order)
 */
__device__ void rk45_step(
    const Vec3& pos0,
    const Vec3& vel0,
    double mass,
    double charge,
    double dt,
    const FieldArrayGPU* fields,
    const Vec3& E_const,
    const Vec3& B_const,
    Vec3& pos_out4,  // 4th order solution
    Vec3& vel_out4,
    Vec3& pos_out5,  // 5th order solution
    Vec3& vel_out5
) {
    // Stage 1: k1 = f(t, y)
    Vec3 k1_vel = vel0;
    Vec3 k1_acc = compute_acceleration_rk45(pos0, vel0, mass, charge, fields, E_const, B_const);
    
    // Stage 2: k2 = f(t + dt/5, y + k1*dt/5)
    Vec3 pos2 = pos0 + k1_vel * (dt / 5.0);
    Vec3 vel2 = vel0 + k1_acc * (dt / 5.0);
    Vec3 k2_vel = vel2;
    Vec3 k2_acc = compute_acceleration_rk45(pos2, vel2, mass, charge, fields, E_const, B_const);
    
    // Stage 3: k3 = f(t + 3*dt/10, y + 3*k1*dt/40 + 9*k2*dt/40)
    Vec3 pos3 = pos0 + k1_vel * (3.0 * dt / 40.0) + k2_vel * (9.0 * dt / 40.0);
    Vec3 vel3 = vel0 + k1_acc * (3.0 * dt / 40.0) + k2_acc * (9.0 * dt / 40.0);
    Vec3 k3_vel = vel3;
    Vec3 k3_acc = compute_acceleration_rk45(pos3, vel3, mass, charge, fields, E_const, B_const);
    
    // Stage 4: k4 = f(t + 4*dt/5, y + 44*k1*dt/45 - 56*k2*dt/15 + 32*k3*dt/9)
    Vec3 pos4 = pos0 + k1_vel * (44.0 * dt / 45.0) + k2_vel * (-56.0 * dt / 15.0) + k3_vel * (32.0 * dt / 9.0);
    Vec3 vel4 = vel0 + k1_acc * (44.0 * dt / 45.0) + k2_acc * (-56.0 * dt / 15.0) + k3_acc * (32.0 * dt / 9.0);
    Vec3 k4_vel = vel4;
    Vec3 k4_acc = compute_acceleration_rk45(pos4, vel4, mass, charge, fields, E_const, B_const);
    
    // Stage 5: k5 = f(t + 8*dt/9, ...)
    Vec3 pos5 = pos0 + k1_vel * (19372.0 * dt / 6561.0) + k2_vel * (-25360.0 * dt / 2187.0) 
              + k3_vel * (64448.0 * dt / 6561.0) + k4_vel * (-212.0 * dt / 729.0);
    Vec3 vel5 = vel0 + k1_acc * (19372.0 * dt / 6561.0) + k2_acc * (-25360.0 * dt / 2187.0) 
              + k3_acc * (64448.0 * dt / 6561.0) + k4_acc * (-212.0 * dt / 729.0);
    Vec3 k5_vel = vel5;
    Vec3 k5_acc = compute_acceleration_rk45(pos5, vel5, mass, charge, fields, E_const, B_const);
    
    // Stage 6: k6 = f(t + dt, ...)
    Vec3 pos6 = pos0 + k1_vel * (9017.0 * dt / 3168.0) + k2_vel * (-355.0 * dt / 33.0) 
              + k3_vel * (46732.0 * dt / 5247.0) + k4_vel * (49.0 * dt / 176.0) 
              + k5_vel * (-5103.0 * dt / 18656.0);
    Vec3 vel6 = vel0 + k1_acc * (9017.0 * dt / 3168.0) + k2_acc * (-355.0 * dt / 33.0) 
              + k3_acc * (46732.0 * dt / 5247.0) + k4_acc * (49.0 * dt / 176.0) 
              + k5_acc * (-5103.0 * dt / 18656.0);
    Vec3 k6_vel = vel6;
    Vec3 k6_acc = compute_acceleration_rk45(pos6, vel6, mass, charge, fields, E_const, B_const);
    
    // 4th order solution (for propagation)
    pos_out4 = pos0 + k1_vel * (35.0 * dt / 384.0) + k3_vel * (500.0 * dt / 1113.0)
             + k4_vel * (125.0 * dt / 192.0) + k5_vel * (-2187.0 * dt / 6784.0) 
             + k6_vel * (11.0 * dt / 84.0);
    vel_out4 = vel0 + k1_acc * (35.0 * dt / 384.0) + k3_acc * (500.0 * dt / 1113.0)
             + k4_acc * (125.0 * dt / 192.0) + k5_acc * (-2187.0 * dt / 6784.0) 
             + k6_acc * (11.0 * dt / 84.0);
    
    // 5th order solution (for error estimation)
    pos_out5 = pos0 + k1_vel * (5179.0 * dt / 57600.0) + k3_vel * (7571.0 * dt / 16695.0)
             + k4_vel * (393.0 * dt / 640.0) + k5_vel * (-92097.0 * dt / 339200.0) 
             + k6_vel * (187.0 * dt / 2100.0) + k6_vel * (dt / 40.0);  // k7 = k6 (FSAL)
    vel_out5 = vel0 + k1_acc * (5179.0 * dt / 57600.0) + k3_acc * (7571.0 * dt / 16695.0)
             + k4_acc * (393.0 * dt / 640.0) + k5_acc * (-92097.0 * dt / 339200.0) 
             + k6_acc * (187.0 * dt / 2100.0) + k6_acc * (dt / 40.0);
}

/**
 * @brief Compute error estimate and new step size
 */
__device__ double compute_error_and_step(
    const Vec3& pos4,
    const Vec3& vel4,
    const Vec3& pos5,
    const Vec3& vel5,
    const RK45Params& params,
    double dt_current
) {
    // Error = max(|pos5 - pos4|, |vel5 - vel4|) scaled by tolerance
    double err_pos_x = fabs(pos5.x - pos4.x) / (params.atol + params.rtol * fabs(pos4.x));
    double err_pos_y = fabs(pos5.y - pos4.y) / (params.atol + params.rtol * fabs(pos4.y));
    double err_pos_z = fabs(pos5.z - pos4.z) / (params.atol + params.rtol * fabs(pos4.z));
    double err_vel_x = fabs(vel5.x - vel4.x) / (params.atol + params.rtol * fabs(vel4.x));
    double err_vel_y = fabs(vel5.y - vel4.y) / (params.atol + params.rtol * fabs(vel4.y));
    double err_vel_z = fabs(vel5.z - vel4.z) / (params.atol + params.rtol * fabs(vel4.z));
    
    double error = fmax(err_pos_x, fmax(err_pos_y, fmax(err_pos_z, 
                   fmax(err_vel_x, fmax(err_vel_y, err_vel_z)))));
    
    // Avoid division by zero
    error = fmax(error, 1e-15);
    
    // New step size: dt_new = dt * safety * (tol/error)^(1/5)
    double factor = params.safety_factor * pow(1.0 / error, 0.2);  // 1/5 = 0.2
    factor = fmin(params.max_step_factor, fmax(params.min_step_factor, factor));
    
    return dt_current * factor;
}

/**
 * @brief RK45 adaptive integration kernel
 * 
 * Each thread integrates one ion from t to t+dt using adaptive substeps.
 */
__global__ void integrate_rk45_batch_kernel(
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
    double dt_max,
    RK45Params params,
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
        
        // Load initial state
        Vec3 pos = {x_in[i], y_in[i], z_in[i]};
        Vec3 vel = {vx_in[i], vy_in[i], vz_in[i]};
        double m = mass[i];
        double q = charge[i];
        
        // Adaptive integration loop
        double t_current = 0.0;
        double dt_sub = dt_max * 0.1;  // Start with small step
        int substeps = 0;
        
        while (t_current < dt_max && substeps < params.max_substeps) {
            // Don't overshoot
            if (t_current + dt_sub > dt_max) {
                dt_sub = dt_max - t_current;
            }
            
            // Take RK45 step
            Vec3 pos4, vel4, pos5, vel5;
            rk45_step(pos, vel, m, q, dt_sub, field_array, E_const, B_const,
                     pos4, vel4, pos5, vel5);
            
            // Compute error and new step size
            double dt_new = compute_error_and_step(pos4, vel4, pos5, vel5, params, dt_sub);
            
            // Check if error is acceptable (error < 1.0)
            double err = 0.0;
            double err_pos_x = fabs(pos5.x - pos4.x) / (params.atol + params.rtol * fabs(pos4.x));
            double err_pos_y = fabs(pos5.y - pos4.y) / (params.atol + params.rtol * fabs(pos4.y));
            double err_pos_z = fabs(pos5.z - pos4.z) / (params.atol + params.rtol * fabs(pos4.z));
            err = fmax(err_pos_x, fmax(err_pos_y, err_pos_z));
            
            if (err <= 1.0) {
                // Accept step (use 4th order solution)
                pos = pos4;
                vel = vel4;
                t_current += dt_sub;
                substeps++;
            }
            
            // Update step size for next iteration
            dt_sub = dt_new;
        }
        
        // Store final state
        x_out[i] = pos.x;
        y_out[i] = pos.y;
        z_out[i] = pos.z;
        vx_out[i] = vel.x;
        vy_out[i] = vel.y;
        vz_out[i] = vel.z;
    }
}

// ============================================================================
// Host-side API
// ============================================================================

void integrate_rk45_batch(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const Vec3& E_field,
    const Vec3& B_field,
    double dt,
    const RK45Params& params,
    cudaStream_t stream
) {
    int N = ions_in.count;
    if (N == 0) return;
    
    // Launch configuration (256 threads per block, enough blocks for all ions)
    int threads = 256;
    int blocks = (N + threads - 1) / threads;
    
    integrate_rk45_batch_kernel<<<blocks, threads, 0, stream>>>(
        ions_in.x, ions_in.y, ions_in.z,
        ions_in.vx, ions_in.vy, ions_in.vz,
        ions_in.mass, ions_in.charge, ions_in.active,
        ions_out.x, ions_out.y, ions_out.z,
        ions_out.vx, ions_out.vy, ions_out.vz,
        ions_out.active,
        nullptr,  // No field array
        E_field, B_field,
        dt, params, N
    );
}

void integrate_rk45_batch_with_fields(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const FieldArrayGPU& field_array,
    double dt,
    const RK45Params& params,
    cudaStream_t stream
) {
    int N = ions_in.count;
    if (N == 0) return;
    
    int threads = 256;
    int blocks = (N + threads - 1) / threads;
    
    integrate_rk45_batch_kernel<<<blocks, threads, 0, stream>>>(
        ions_in.x, ions_in.y, ions_in.z,
        ions_in.vx, ions_in.vy, ions_in.vz,
        ions_in.mass, ions_in.charge, ions_in.active,
        ions_out.x, ions_out.y, ions_out.z,
        ions_out.vx, ions_out.vy, ions_out.vz,
        ions_out.active,
        &field_array,
        Vec3{0, 0, 0}, Vec3{0, 0, 0},  // Not used with field array
        dt, params, N
    );
}

} // namespace gpu
} // namespace icarion

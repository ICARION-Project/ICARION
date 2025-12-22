// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "integrate_rk4_batch.cuh"
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/fields/FieldArrayGPU.h"
#include "core/gpu/fields/FieldArrayGPU_kernels.cuh"
#include <cuda_runtime.h>

namespace icarion {
namespace gpu {

/**
 * @brief Compute acceleration for a single ion on GPU (internal helper)
 * 
 * Uses field interpolation if field array is provided, otherwise uses constant fields.
 * 
 * @param pos Position [m]
 * @param vel Velocity [m/s]
 * @param mass Mass [kg]
 * @param charge Charge [C]
 * @param fields Field array for interpolation (nullptr for constant fields)
 * @param E_const Constant electric field [V/m] (used if fields == nullptr)
 * @param B_const Constant magnetic field [T] (used if fields == nullptr)
 * @return Acceleration [m/s²]
 */
__device__ Vec3 compute_acceleration(
    const Vec3& pos,
    const Vec3& vel,
    double mass,
    double charge,
    const FieldArrayGPU* fields,
    const Vec3& E_const,
    const Vec3& B_const
) {
    // Evaluate fields at position
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
    
    // v × B (cross product)
    Vec3 v_cross_B = {
        vel.y * B_field.z - vel.z * B_field.y,
        vel.z * B_field.x - vel.x * B_field.z,
        vel.x * B_field.y - vel.y * B_field.x
    };
    Vec3 F_magnetic = v_cross_B * charge;
    
    Vec3 F_total = F_electric + F_magnetic;
    
    // a = F/m
    return F_total / mass;
}

/**
 * @brief RK4 integration kernel (grid-stride loop for any N)
 * 
 * Each thread processes multiple ions using grid-stride pattern:
 * for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += gridDim.x * blockDim.x)
 * 
 * This allows optimal occupancy regardless of N.
 */
__global__ void integrate_rk4_batch_kernel(
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
    
    // Fields (interpolated or constant)
    const FieldArrayGPU* field_array,
    Vec3 E_const,
    Vec3 B_const,
    
    // Integration parameters
    double dt,
    int N,
    DeviceDamping damping
) {
    // Grid-stride loop
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; i < N; i += gridDim.x * blockDim.x) {
        // Always copy active flag (ions don't change active state during integration)
        active_out[i] = active_in[i];
        
        // Skip inactive ions
        if (!active_in[i]) {
            // Copy unchanged state
            x_out[i] = x_in[i];
            y_out[i] = y_in[i];
            z_out[i] = z_in[i];
            vx_out[i] = vx_in[i];
            vy_out[i] = vy_in[i];
            vz_out[i] = vz_in[i];
            continue;
        }
        
        // Load ion state
        Vec3 pos0 = {x_in[i], y_in[i], z_in[i]};
        Vec3 vel0 = {vx_in[i], vy_in[i], vz_in[i]};
        double m = mass[i];
        double q = charge[i];
        
        // RK4 Stage 1: k1 = f(t, y)
        Vec3 k1_vel = vel0;
        Vec3 k1_acc = compute_acceleration(pos0, vel0, m, q, field_array, E_const, B_const);
        
        // RK4 Stage 2: k2 = f(t + dt/2, y + k1*dt/2)
        Vec3 pos2 = pos0 + k1_vel * (dt * 0.5);
        Vec3 vel2 = vel0 + k1_acc * (dt * 0.5);
        Vec3 k2_vel = vel2;
        Vec3 k2_acc = compute_acceleration(pos2, vel2, m, q, field_array, E_const, B_const);
        
        // RK4 Stage 3: k3 = f(t + dt/2, y + k2*dt/2)
        Vec3 pos3 = pos0 + k2_vel * (dt * 0.5);
        Vec3 vel3 = vel0 + k2_acc * (dt * 0.5);
        Vec3 k3_vel = vel3;
        Vec3 k3_acc = compute_acceleration(pos3, vel3, m, q, field_array, E_const, B_const);
        
        // RK4 Stage 4: k4 = f(t + dt, y + k3*dt)
        Vec3 pos4 = pos0 + k3_vel * dt;
        Vec3 vel4 = vel0 + k3_acc * dt;
        Vec3 k4_vel = vel4;
        Vec3 k4_acc = compute_acceleration(pos4, vel4, m, q, field_array, E_const, B_const);
        
        // RK4 Final: y(t+dt) = y(t) + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        Vec3 pos_new = pos0 + (k1_vel + k2_vel * 2.0 + k3_vel * 2.0 + k4_vel) * (dt / 6.0);
        Vec3 vel_new = vel0 + (k1_acc + k2_acc * 2.0 + k3_acc * 2.0 + k4_acc) * (dt / 6.0);
        
        // Apply linear damping (v <- v * exp(-nu*dt)) if enabled
        if (damping.enabled) {
            float nu = damping.nu_per_ion ? damping.nu_per_ion[i] : damping.nu_const;
            if (nu > 0.0f) {
                double factor = exp(-static_cast<double>(nu) * dt);
                vel_new = vel_new * factor;
            }
        }

        // Store result
        x_out[i] = pos_new.x;
        y_out[i] = pos_new.y;
        z_out[i] = pos_new.z;
        vx_out[i] = vel_new.x;
        vy_out[i] = vel_new.y;
        vz_out[i] = vel_new.z;
    }
}

// Host function to launch kernel with constant fields
void integrate_rk4_batch(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const Vec3& E_field,
    const Vec3& B_field,
    double dt,
    const DeviceDamping& damping,
    cudaStream_t stream
) {
    if (ions_in.count == 0) {
        return;
    }
    
    int N = ions_in.count;
    
    // Compute optimal grid dimensions
    int threads_per_block = 256;  // Good occupancy for most GPUs
    int blocks = (N + threads_per_block - 1) / threads_per_block;
    
    // Limit blocks to avoid excessive scheduling overhead
    // Grid-stride loop handles large N efficiently
    blocks = min(blocks, 2048);
    
    // Launch kernel with nullptr field array (use constant fields)
    integrate_rk4_batch_kernel<<<blocks, threads_per_block, 0, stream>>>(
        // Input
        ions_in.x, ions_in.y, ions_in.z,
        ions_in.vx, ions_in.vy, ions_in.vz,
        ions_in.mass, ions_in.charge, ions_in.active,
        
        // Output
        ions_out.x, ions_out.y, ions_out.z,
        ions_out.vx, ions_out.vy, ions_out.vz,
        ions_out.active,
        
        // Fields (constant)
        nullptr, E_field, B_field,
        
        // Parameters
        dt, N, damping
    );
    
    // Check for kernel launch errors
    cudaError_t err = cudaGetLastError();
    if (err != cudaSuccess) {
        throw std::runtime_error(
            std::string("RK4 kernel launch failed: ") + cudaGetErrorString(err)
        );
    }
}

// Host function to launch kernel with field interpolation
void integrate_rk4_batch_with_fields(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    const FieldArrayGPU& field_array,
    double dt,
    const DeviceDamping& damping,
    cudaStream_t stream
) {
    if (ions_in.count == 0) {
        return;
    }
    
    int N = ions_in.count;
    
    // Compute optimal grid dimensions
    int threads_per_block = 256;
    int blocks = (N + threads_per_block - 1) / threads_per_block;
    blocks = min(blocks, 2048);
    
    // Copy field array to device memory
    FieldArrayGPU* d_field_array;
    cudaError_t err = cudaMalloc(&d_field_array, sizeof(FieldArrayGPU));
    if (err != cudaSuccess) {
        throw std::runtime_error("Failed to allocate device memory for field array");
    }
    
    err = cudaMemcpyAsync(d_field_array, &field_array, sizeof(FieldArrayGPU),
                          cudaMemcpyHostToDevice, stream);
    if (err != cudaSuccess) {
        cudaFree(d_field_array);
        throw std::runtime_error("Failed to copy field array to device");
    }
    
    // Launch kernel with field interpolation
    Vec3 zero = {0.0, 0.0, 0.0};
    integrate_rk4_batch_kernel<<<blocks, threads_per_block, 0, stream>>>(
        // Input
        ions_in.x, ions_in.y, ions_in.z,
        ions_in.vx, ions_in.vy, ions_in.vz,
        ions_in.mass, ions_in.charge, ions_in.active,
        
        // Output
        ions_out.x, ions_out.y, ions_out.z,
        ions_out.vx, ions_out.vy, ions_out.vz,
        ions_out.active,
        
        // Fields (interpolated)
        d_field_array, zero, zero,
        
        // Parameters
        dt, N, damping
    );
    
    // Synchronize stream before freeing device memory
    cudaStreamSynchronize(stream);
    cudaFree(d_field_array);
    
    // Check for kernel errors
    err = cudaGetLastError();
    if (err != cudaSuccess) {
        throw std::runtime_error(
            std::string("RK4 kernel with field interpolation failed: ") + cudaGetErrorString(err)
        );
    }
}

} // namespace gpu
} // namespace icarion

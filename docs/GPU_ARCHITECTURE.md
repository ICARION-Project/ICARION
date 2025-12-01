# GPU Architecture Details

**Purpose:** Detailed GPU implementation documentation for ICARION

This document contains the detailed GPU architecture information that was previously in ARCHITECTURE.md. For high-level overview, see [ARCHITECTURE_OVERVIEW.md](ARCHITECTURE_OVERVIEW.md).

## GPU Acceleration Implementation

**Status:** Core features implemented (v1.0)

### Hybrid CPU/GPU Architecture

ICARION uses automatic dispatch based on ion count (threshold: 5000 active ions):

```cpp
// In SimulationEngine::integrate_timestep()
bool use_gpu = (ions.size() >= gpu_threshold_) && gpu_helper_ && gpu_context_;

if (use_gpu) {
    gpu_helper_->integrate_batch(ions, dt, t);
} else {
    // CPU integration
    for (auto& ion : ions) {
        integration_strategy_->integrate_single_step(ion, dt, t, force_context_);
    }
}
```

### GPU Memory Management

**Structure-of-Arrays (SoA) Layout:**
```cpp
// CPU: Array-of-Structures (AoS)
struct IonState {
    Vec3 pos, vel;
    double mass_kg, ion_charge_C;
    bool is_active;
};
std::vector<IonState> ions;

// GPU: Structure-of-Arrays (SoA) 
struct IonStateGPU {
    float* pos_x, *pos_y, *pos_z;
    float* vel_x, *vel_y, *vel_z;
    float* mass_kg, *charge_C;
    bool* is_active;
};
```

**Benefits:**
- Coalesced memory access (all threads access same field)
- SIMD-friendly operations
- Reduced memory bandwidth requirements

### GPU Integration Kernels

**RK4 Integration Kernel:**
```cuda
__global__ void integrate_rk4_kernel(
    IonStateGPU ions_in,
    IonStateGPU ions_out,
    int n_ions,
    double dt,
    double t,
    FieldArrayGPU field_data
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_ions || !ions_in.is_active[idx]) return;
    
    // Load ion state
    Vec3 pos = {ions_in.pos_x[idx], ions_in.pos_y[idx], ions_in.pos_z[idx]};
    Vec3 vel = {ions_in.vel_x[idx], ions_in.vel_y[idx], ions_in.vel_z[idx]};
    double mass = ions_in.mass_kg[idx];
    double charge = ions_in.charge_C[idx];
    
    // RK4 integration
    Vec3 k1_pos = vel;
    Vec3 k1_vel = compute_acceleration_gpu(pos, vel, mass, charge, t, field_data);
    
    Vec3 k2_pos = vel + 0.5 * dt * k1_vel;
    Vec3 k2_vel = compute_acceleration_gpu(pos + 0.5 * dt * k1_pos, vel + 0.5 * dt * k1_vel, 
                                         mass, charge, t + 0.5 * dt, field_data);
    
    // k3, k4 calculations...
    
    // Update ion state
    pos += (dt/6.0) * (k1_pos + 2*k2_pos + 2*k3_pos + k4_pos);
    vel += (dt/6.0) * (k1_vel + 2*k2_vel + 2*k3_vel + k4_vel);
    
    // Store results
    ions_out.pos_x[idx] = pos.x;
    ions_out.pos_y[idx] = pos.y;
    ions_out.pos_z[idx] = pos.z;
    ions_out.vel_x[idx] = vel.x;
    ions_out.vel_y[idx] = vel.y;
    ions_out.vel_z[idx] = vel.z;
}
```

### GPU Collision Implementation

**HSS Collision Kernel:**
```cuda
__global__ void hss_collision_kernel(
    IonStateGPU ions,
    int n_ions,
    double dt,
    curandState* rng_states,
    CollisionParams params
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_ions || !ions.is_active[idx]) return;
    
    curandState* rng = &rng_states[idx];
    
    // Load ion state
    Vec3 vel = {ions.vel_x[idx], ions.vel_y[idx], ions.vel_z[idx]};
    double mass = ions.mass_kg[idx];
    
    // Collision probability
    double v_rel = length(vel - params.gas_velocity);
    double collision_freq = params.cross_section_m2 * params.gas_density * v_rel;
    double collision_prob = 1.0 - exp(-collision_freq * dt);
    
    if (curand_uniform(rng) < collision_prob) {
        // Sample post-collision velocity from Maxwell-Boltzmann
        Vec3 v_th = sample_thermal_velocity_gpu(params.temperature_K, params.gas_mass_kg, rng);
        
        // Update ion velocity
        ions.vel_x[idx] = v_th.x + params.gas_velocity.x;
        ions.vel_y[idx] = v_th.y + params.gas_velocity.y;
        ions.vel_z[idx] = v_th.z + params.gas_velocity.z;
    }
}
```

### Field Array GPU Implementation

**Texture Memory for Field Interpolation:**
```cuda
// Bind field arrays to texture memory for fast interpolation
texture<float, 3, cudaReadModeElementType> tex_field_Ex;
texture<float, 3, cudaReadModeElementType> tex_field_Ey;
texture<float, 3, cudaReadModeElementType> tex_field_Ez;

__device__ Vec3 interpolate_field_gpu(Vec3 pos, FieldArrayGPU field_data) {
    // Convert world coordinates to texture coordinates
    float tx = (pos.x - field_data.origin.x) / field_data.spacing.x;
    float ty = (pos.y - field_data.origin.y) / field_data.spacing.y;
    float tz = (pos.z - field_data.origin.z) / field_data.spacing.z;
    
    // Trilinear interpolation via texture hardware
    float Ex = tex3D(tex_field_Ex, tx, ty, tz) * field_data.scale_factor;
    float Ey = tex3D(tex_field_Ey, tx, ty, tz) * field_data.scale_factor;
    float Ez = tex3D(tex_field_Ez, tx, ty, tz) * field_data.scale_factor;
    
    return {Ex, Ey, Ez};
}
```

### GPU Error Handling

**Automatic Fallback Strategy:**
```cpp
bool GPUIntegrationHelper::integrate_batch(std::vector<IonState>& ions, double dt, double t) {
    try {
        // 1. Copy data to GPU
        copy_ions_to_gpu(ions);
        
        // 2. Launch integration kernel
        launch_integration_kernel(ions.size(), dt, t);
        
        // 3. Check for errors
        CUDA_CHECK(cudaDeviceSynchronize());
        
        // 4. Copy results back
        copy_ions_from_gpu(ions);
        
        return true;
        
    } catch (const CudaException& e) {
        spdlog::warn("GPU integration failed: {}", e.what());
        spdlog::info("Falling back to CPU integration");
        return false;  // Caller will use CPU path
    }
}
```

### Performance Characteristics

**GPU Benefits:**
- Parallel processing of thousands of ions simultaneously
- Coalesced memory access patterns reduce bandwidth requirements
- Texture memory acceleration for field array interpolation

**When GPU Helps:**
- Large ion populations (threshold-based dispatch)
- Compute-bound simulations (complex force calculations)
- Long simulation times where setup overhead amortizes

**When GPU May Not Help:**
- Small ion populations (overhead dominates)
- Memory-bound simulations (frequent data transfer)
- Simple analytical fields (CPU already fast)

**Performance varies significantly with:**
- GPU hardware (compute capability, memory bandwidth)
- Simulation complexity (force models, field types)
- Problem size (ion count, integration steps)
- System configuration (PCIe bandwidth, CPU performance)

Benchmark your specific use case to determine optimal settings.

### GPU Memory Layout

**Optimization Strategies:**
1. **Coalesced Access**: All threads in warp access consecutive memory
2. **Bank Conflict Avoidance**: Shared memory layout optimized
3. **Register Usage**: Minimize register spilling to local memory
4. **Occupancy**: Balance threads/block vs. shared memory usage

### CUDA Kernel Launch Configuration

**Block/Grid Sizing:**
```cpp
void launch_integration_kernel(int n_ions, double dt, double t) {
    // Optimal block size (determined empirically)
    constexpr int BLOCK_SIZE = 256;
    
    // Grid size to cover all ions
    int grid_size = (n_ions + BLOCK_SIZE - 1) / BLOCK_SIZE;
    
    // Launch kernel
    integrate_rk4_kernel<<<grid_size, BLOCK_SIZE>>>(
        gpu_ions_in_, gpu_ions_out_, n_ions, dt, t, gpu_field_data_
    );
    
    // Check for launch errors
    CUDA_CHECK(cudaGetLastError());
}
```

### GPU Debugging and Profiling

**CUDA Error Checking:**
```cpp
#define CUDA_CHECK(call) do { \
    cudaError_t error = call; \
    if (error != cudaSuccess) { \
        throw CudaException(cudaGetErrorString(error), __FILE__, __LINE__); \
    } \
} while(0)
```

**Profiling Tools:**
- `nvprof`: Legacy CUDA profiler
- `nsys`: Nsight Systems timeline profiler
- `ncu`: Nsight Compute kernel profiler
- Custom timing: `cudaEventRecord()` for kernel timing

### Testing Strategy

**GPU Validation Tests:**
- **Correctness**: CPU/GPU result parity within tolerance
- **Performance**: Speedup measurement vs. CPU baseline
- **Robustness**: Error handling and fallback mechanisms
- **Memory**: Leak detection and proper cleanup

**Example Test:**
```cpp
TEST(GPUIntegrationTest, RK4_Parity) {
    std::vector<IonState> ions_cpu = create_test_ions(10000);
    std::vector<IonState> ions_gpu = ions_cpu;
    
    // CPU integration
    cpu_integrator->integrate_batch(ions_cpu, dt, t);
    
    // GPU integration
    gpu_helper->integrate_batch(ions_gpu, dt, t);
    
    // Compare results
    for (size_t i = 0; i < ions_cpu.size(); ++i) {
        EXPECT_NEAR(ions_cpu[i].pos.x, ions_gpu[i].pos.x, 1e-6);
        EXPECT_NEAR(ions_cpu[i].pos.y, ions_gpu[i].pos.y, 1e-6);
        EXPECT_NEAR(ions_cpu[i].pos.z, ions_gpu[i].pos.z, 1e-6);
    }
}
```

---

For more GPU implementation details, see:
- `src/core/gpu/` - GPU implementation files
- `tests/gpu/` - GPU-specific tests
- `docs/DEVELOPERS_GUIDE.md` - GPU development guide
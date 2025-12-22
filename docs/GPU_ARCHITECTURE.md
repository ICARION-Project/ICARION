# GPU Architecture Details

**Purpose:** Detailed GPU implementation documentation for ICARION

This document contains the detailed GPU architecture information for ICARION. For high-level overview, see [ARCHITECTURE.md](ARCHITECTURE.md).

## GPU Acceleration Implementation

**Status:** Core features implemented (v1.0) with CPU fallback; space-charge and GPU boundary helpers exist but are not called from the main loop yet.

### Hybrid CPU/GPU Architecture

ICARION uses automatic dispatch inside `SimulationEngine::process_timestep()`:

- Integration: `try_gpu_integration` uses total ion count and a cached integrator type. Default threshold is 5000 ions, halved for Boris (cheap kernel). If below threshold, unknown integrator, or GPU unavailable, CPU path runs.
- Collisions: `try_gpu_collisions` uses **active** ion count with a 5000 threshold. Falls back to CPU on failure or below threshold.
- Space charge: P³M helper exists (`try_gpu_space_charge`) but is not invoked from the main loop yet.
- Boundary checks: Helper exists for absorption/cylindrical domains but is not wired into the main loop.

### GPU Memory Management

**AoS Upload (current implementation):**
- CPU holds ions as `std::vector<IonState>` (AoS). GPU helpers upload AoS buffers internally; SoA conversion is not used for GPU paths yet.
- `try_gpu_boundary_check` converts SoA (`IonEnsemble`) back to AoS before dispatching to the helper.

**Implications:**
- Coalescing is limited by AoS layout; future SoA upload would improve bandwidth.
- Single-domain field provider (domain 0) is passed to GPU integration helper when available.

### GPU Integration Kernels

**Integration Helpers:**
- Batch kernels exist for RK4, RK45 (with atol/rtol), and Boris. Dispatch is selected without repeated dynamic casts (type cached).
- Kernel signatures are helper-specific; code above is schematic. Refer to `GPUIntegrationHelper` for exact interfaces.

### GPU Collision Implementation

**Collision Helper (HSS/EHSS):**
- Uses active ion count threshold (default 5000).
- EHSS geometry upload is noted as TODO; helper currently supports HSS/EHSS dispatch with CPU fallback.

### Field Array GPU Implementation

- GPU integration helper accepts a field provider pointer (when present in domain 0’s ForceRegistry). Field upload details are encapsulated in the helper; texture binding is implementation-specific and not exposed at the interface.

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
// Kernel launch sizing is helper-specific; 256 threads/block is a common default.
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

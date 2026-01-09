# GPU Architecture Details

**Purpose:** Detailed GPU implementation documentation for ICARION

This document contains the detailed GPU architecture information for ICARION. For high-level overview, see [ARCHITECTURE.md](ARCHITECTURE.md).

## GPU Acceleration Implementation

**Status:** Core GPU paths are implemented with CPU fallback. Space charge can use `SpaceChargeGPUModel` when enabled; the boundary-check helper exists but is not wired into the main loop.

### Hybrid CPU/GPU Architecture

ICARION uses automatic dispatch inside `SimulationEngine::process_timestep()` via the strategy/handler interfaces:

- Integration: `IntegrationStrategyFactory` wraps the selected CPU strategy in `GPUIntegrationStrategy` when `simulation.enable_gpu` is true and ICARION is built with CUDA. Batch dispatch happens inside `SimulationEngine::perform_integration()` when all active ions share the same `dt`. The GPU path requires a single domain with exactly one `ElectricFieldForce` and a field provider; it falls back on space charge, magnetic forces, damping unless GPU damping is enabled, unsupported force mixes, or missing field providers. `GridFieldProvider` + `FieldArray` is required for non-zero fields; non-grid or snapshot-based providers run with zero fields, so avoid GPU dispatch in those cases. GPU damping must be enabled explicitly via `GPUIntegrationStrategy` setters (not wired from config).
- Collisions: `CollisionHandlerFactory` can return a `GPUCollisionHandler` for HSS/EHSS when `simulation.enable_gpu` is true. `SimulationEngine::perform_collisions()` groups ions per domain and calls `handle_batch()` when the handler advertises `supports_batch()` and the domain has uniform `dt`. The handler falls back to CPU below its threshold or on GPU errors.
- Reactions: `ReactionHandlerFactory` can return a `GPUReactionHandler` (wrapper around `GPUReactionBackend`) when GPU is enabled; `SimulationEngine::perform_reactions()` uses the batch hook when available and `dt` is uniform across active ions.
- Space charge: `SpaceChargeModelFactory` may construct `SpaceChargeGPUModel` when `physics.enable_space_charge_gpu` is set and CUDA is available; `SimulationEngine::update_space_charge_models()` calls `update_fields()` on the chosen model.
- Boundary checks: A helper exists for absorption/cylindrical domains but is not currently dispatched in the main loop.

### GPU Memory Management

**Host-side layout (current implementation):**
- The simulation state lives in `IonEnsemble` (SoA). GPU integration/collision paths gather active ions into temporary `std::vector<IonState>` (AoS) before upload. GPU reactions operate directly on the SoA ensemble.
- Device buffers use `IonStateGPU` (SoA). There is no direct SoA upload from `IonEnsemble` yet.

**Implications:**
- Coalescing is limited by AoS layout; future SoA upload would improve bandwidth.
- The GPU integration helper is passed the selected domain’s field provider; only `GridFieldProvider` with a `FieldArray` is supported. Non-grid providers are treated as unsupported and result in zero-field integration, so GPU dispatch should be avoided in those cases.

### GPU Integration Kernels

**Integration Helpers:**
- Batch kernels exist for RK4, RK45 (adaptive), and Boris. Dispatch is selected by the `GPUIntegrationStrategy::Kind` set at construction.
- Kernel signatures are helper-specific; refer to `GPUIntegrationHelper` for exact interfaces and parameter defaults.

### GPU Collision Implementation

**Collision Helper (HSS/EHSS):**
- Uses active ion count threshold (default 5000).
- EHSS geometry upload is supported when a geometry map is provided; the helper falls back to CPU on errors or missing prerequisites.

### Field Array GPU Implementation

- GPU integration helper accepts a field provider pointer from the selected domain. Only `GridFieldProvider` + `FieldArray` are uploaded; providers using snapshots are skipped and non-grid providers are treated as unsupported (zero fields).

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
        
    } catch (const std::exception& e) {
        // Log and fall back to CPU
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
using Catch::Matchers::WithinAbs;

TEST_CASE("RK4: CPU/GPU parity - Free particle", "[gpu][rk4][parity]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }

    std::vector<IonState> ions_cpu = create_test_ions(10000);
    std::vector<IonState> ions_gpu = ions_cpu;

    // CPU integration (per-ion)
    RK4Strategy cpu_integrator;
    auto ensemble_cpu = core::IonEnsemble::from_legacy(ions_cpu);
    for (size_t i = 0; i < ensemble_cpu.size(); ++i) {
        cpu_integrator.step(ensemble_cpu, i, 0.0, dt, force_registry);
    }

    // GPU integration
    auto context = GPUContext::create(0);
    auto helper = GPUIntegrationHelper::create(*context, 100);
    REQUIRE(helper);
    helper->integrate_batch_rk4(ions_gpu, dt, 0.0);

    // Compare results
    for (size_t i = 0; i < ions_cpu.size(); ++i) {
        REQUIRE_THAT(ions_gpu[i].pos.x, WithinAbs(ions_cpu[i].pos.x, 1e-6));
    }
}
```

---

For more GPU implementation details, see:
- `src/core/gpu/` - GPU implementation files
- `tests/gpu/` - GPU-specific tests
- `docs/DEVELOPERS_GUIDE.md` - GPU development guide

# GPU Acceleration - Remaining Work

**Date:** 2025-11-29  
**Current Status:** Phase 4+5 Complete (Basic RK4 Integration)  
**Branch:** feature/gpu-acceleration

---

## ✅ Completed (Phase 1-5)

### Phase 1: GPU Infrastructure ✅
- ✅ GPUContext (CUDA device management, RAII)
- ✅ GPUMemoryPool (buffer allocation/reuse)
- ✅ IonState_GPU (SoA layout)
- ✅ AoS ↔ SoA conversion (upload/download)

### Phase 2: GPU Batch Integration ✅
- ✅ `integrate_rk4_batch.cu` - RK4 kernel (grid-stride loop)
- ✅ GPUIntegrationHelper (batch API)
- ✅ Zero-field integration working

### Phase 4: SimulationEngine Integration ✅
- ✅ `try_gpu_integration()` method (5 lines)
- ✅ Automatic CPU/GPU dispatch (N >= 5000 threshold)
- ✅ Graceful fallback to CPU
- ✅ Statistics logging

### Phase 5: Testing & Validation ✅
- ✅ test_gpu_integration.cpp (6 tests passing)
- ✅ test_cpu_gpu_parity.cpp (6 tests passing)
- ✅ 53/53 total tests passing
- ✅ Kernel bugs fixed (active flag, buffer initialization)

**Lines of Code:** ~1250 lines GPU infrastructure

---

## 🚧 Remaining Work (Phase 6-9)

### Phase 6: GPU Field Interpolation 🔴 HIGH PRIORITY
**Effort:** ~2-3 days  
**Impact:** Critical for realistic simulations

#### 6.1 Field Grid Upload
**Files:**
- `src/core/gpu/FieldArrayGPU.h` (new)
- `src/core/gpu/FieldArrayGPU.cu` (new)

**Tasks:**
- [ ] Upload field grids (Ex, Ey, Ez, Bx, By, Bz) to GPU
- [ ] Support 1D/2D/3D field arrays
- [ ] Use CUDA texture memory for fast interpolation
- [ ] Handle field updates (dynamic fields, RF waveforms)

**Implementation:**
```cuda
struct FieldArrayGPU {
    cudaTextureObject_t Ex_tex, Ey_tex, Ez_tex;
    cudaTextureObject_t Bx_tex, By_tex, Bz_tex;
    Vec3 origin;
    Vec3 spacing;
    int nx, ny, nz;
};

__device__ Vec3 interpolate_E_field(
    const FieldArrayGPU& fields,
    const Vec3& pos
);
```

**Testing:**
- [ ] Uniform field → analytical comparison
- [ ] Quadrupole field → CPU parity
- [ ] IMS drift field → mobility validation

#### 6.2 Integrate Fields into RK4 Kernel
**Files:**
- `src/core/gpu/integrate_rk4_batch.cu` (modify)
- `src/core/gpu/GPUIntegrationHelper.cpp` (modify)

**Tasks:**
- [ ] Pass FieldArrayGPU to kernel
- [ ] Replace zero fields with interpolation
- [ ] Update GPUIntegrationHelper to extract fields from ForceRegistry
- [ ] Test with real instrument configs (IMS, Orbitrap, LQIT)

**Current TODO Location:**
```cpp
// src/core/gpu/GPUIntegrationHelper.cpp:80
// TODO Phase 3: Extract fields from ForceRegistry
// For now, use zero fields (testing phase)
Vec3 E_field = {0.0, 0.0, 0.0};
Vec3 B_field = {0.0, 0.0, 0.0};
```

---

### Phase 7: Additional Integrators 🟡 MEDIUM PRIORITY
**Effort:** ~1-2 days  
**Impact:** Required for magnetic instruments (Orbitrap, FTICR)

#### 7.1 Boris Pusher (Magnetic Fields)
**Files:**
- `src/core/gpu/integrate_boris_batch.cu` (new)
- `src/core/gpu/integrate_boris_batch.cuh` (new)

**Why Boris?**
- Better energy conservation for strong B-fields
- Standard in plasma physics
- Current RK4 has poor stability for Lorentz force

**Implementation:**
```cuda
__global__ void integrate_boris_batch_kernel(
    /* inputs */
    const double* x_in, const double* y_in, const double* z_in,
    const double* vx_in, const double* vy_in, const double* vz_in,
    /* fields */
    const FieldArrayGPU& E_field,
    const FieldArrayGPU& B_field,
    /* params */
    double dt, int N
);
```

**Algorithm:**
1. Half-step velocity with E-field: v⁻ = v + (q/m)E·dt/2
2. Rotate velocity with B-field: v⁺ = rotate(v⁻, B)
3. Half-step velocity with E-field: v' = v⁺ + (q/m)E·dt/2
4. Update position: x' = x + v'·dt

**Testing:**
- [ ] Circular motion in B-field (analytical)
- [ ] Orbitrap frequency validation
- [ ] FTICR cyclotron validation

#### 7.2 RK45 Adaptive Integrator
**Files:**
- `src/core/gpu/integrate_rk45_batch.cu` (new)
- `src/core/gpu/integrate_rk45_batch.cuh` (new)

**Challenges:**
- Variable timesteps → hard to parallelize
- Error estimation per ion
- Possible solution: Fixed timestep mode with error monitoring

**Priority:** LOW (most ICARION sims use fixed timestep)

---

### Phase 8: GPU Collision Handling 🟡 MEDIUM PRIORITY
**Effort:** ~3-4 days  
**Impact:** 20-30% of CPU time in collision-heavy sims

#### 8.1 EHSS Collision Kernel
**Files:**
- `src/core/gpu/collisions/ehss_collision_kernel.cu` (new)
- `src/core/gpu/collisions/GPUCollisionHandler.h` (new)
- `src/core/gpu/collisions/GPUCollisionHandler.cpp` (new)

**Tasks:**
- [ ] Port EHSS algorithm to GPU
- [ ] GPU random number generation (cuRAND)
- [ ] Velocity scattering (isotropic/forward)
- [ ] Energy loss calculation
- [ ] Statistics accumulation (collision count, energy loss)

**Challenges:**
- **RNG State:** Each ion needs persistent RNG state
  - Solution: Allocate RNG state array on GPU
  - Initialize once, update per collision
- **Atomic Operations:** Collision statistics
  - Solution: Per-block reduction, then atomic add
- **Divergence:** Active vs inactive ions
  - Solution: Warp-aware branching

**Implementation:**
```cuda
__global__ void ehss_collision_kernel(
    /* ion state */
    IonStateGPU ions,
    /* collision params */
    double dt,
    double pressure,
    double temperature,
    /* RNG state */
    curandState* rng_states,
    /* output stats */
    uint64_t* collision_counts
);
```

**Testing:**
- [ ] Thermalization to Maxwell-Boltzmann (300K)
- [ ] Mean free path validation
- [ ] Energy distribution after 1000 collisions

#### 8.2 HSS Collision Kernel
**Files:**
- `src/core/gpu/collisions/hss_collision_kernel.cu` (new)

**Similar to EHSS, but:**
- Hard-sphere model (no energy loss)
- Simpler scattering angle calculation
- Better for high-energy ions

---

### Phase 9: GPU Boundary Checking 🟢 LOW PRIORITY
**Effort:** ~1 day  
**Impact:** ~5% of CPU time

#### 9.1 Geometric Boundary Kernels
**Files:**
- `src/core/gpu/boundaries/boundary_kernels.cu` (new)
- `src/core/gpu/boundaries/GPUBoundaryChecker.h` (new)

**Tasks:**
- [ ] Cylinder boundary check
- [ ] Box boundary check
- [ ] Aperture transition check
- [ ] Domain assignment

**Implementation:**
```cuda
__device__ bool is_inside_cylinder(
    const Vec3& pos,
    double r_max,
    double z_min,
    double z_max
);

__global__ void check_boundaries_kernel(
    IonStateGPU ions,
    BoundaryParams params,
    bool* active_out,
    int* domain_out
);
```

**Testing:**
- [ ] Ion termination at walls
- [ ] Domain transitions through apertures
- [ ] Correct domain assignment

---

## 📊 Priority & Timeline

### Immediate (Next 1-2 weeks)
1. **Phase 6: Field Interpolation** 🔴
   - Critical for any realistic simulation
   - Blocks testing with real instruments
   - **Target:** Complete by Dec 6, 2025

### Short-term (Next 2-4 weeks)
2. **Phase 7.1: Boris Pusher** 🟡
   - Required for Orbitrap, FTICR
   - Better than RK4 for magnetic fields
   - **Target:** Complete by Dec 13, 2025

3. **Phase 8.1: EHSS Collisions** 🟡
   - Significant speedup for IMS, drift tube
   - Most common collision model
   - **Target:** Complete by Dec 20, 2025

### Medium-term (1-2 months)
4. **Phase 8.2: HSS Collisions** 🟢
   - Completes collision coverage
   - **Target:** Q1 2026

5. **Phase 9: Boundaries** 🟢
   - Minor speedup, but completes GPU coverage
   - **Target:** Q1 2026

### Long-term (Future)
6. **Phase 7.2: RK45** ⚪
   - Nice-to-have, not critical
   - Most sims use fixed timestep
   - **Target:** Q2 2026 or later

---

## 🎯 Success Criteria

### Phase 6 (Fields)
- ✅ IMS simulation matches CPU results (mobility, arrival time)
- ✅ Orbitrap frequency within 0.1% of CPU
- ✅ LQIT stability diagram matches CPU

### Phase 7 (Boris)
- ✅ Circular orbit in B-field (analytical comparison)
- ✅ Energy conservation < 1e-6 over 1000 orbits
- ✅ Orbitrap and FTICR benchmarks pass

### Phase 8 (Collisions)
- ✅ Thermalization to 300K within 100 µs
- ✅ Maxwell-Boltzmann distribution (KS test p > 0.05)
- ✅ Mean free path within 1% of theory

### Phase 9 (Boundaries)
- ✅ 100% domain transitions match CPU
- ✅ No false positives/negatives in boundary checks

---

## 📈 Expected Performance Gains

| Component | CPU Baseline | GPU Target | Speedup | Status |
|-----------|-------------|------------|---------|--------|
| **Integration (RK4)** | 100 ms | 5-10 ms | 10-20× | ✅ Done |
| **Integration (Boris)** | 120 ms | 5-10 ms | 12-24× | 🔴 TODO |
| **Collisions (EHSS)** | 50 ms | 5-10 ms | 5-10× | 🔴 TODO |
| **Field Evaluation** | 10 ms | 1-2 ms | 5-10× | 🔴 TODO |
| **Boundaries** | 5 ms | 0.5 ms | 10× | 🔴 TODO |
| **Total (N=100k)** | 285 ms | 20-30 ms | **10-15×** | 50% Done |

**Note:** Speedups are for N=100k ions. Threshold effects reduce gains for N<5000.

---

## 🔧 Infrastructure Needed

### For Phase 6 (Fields)
- [ ] CUDA texture memory support
- [ ] 3D interpolation routines
- [ ] Field update pipeline (CPU→GPU)
- [ ] Support for time-varying fields (RF waveforms)

### For Phase 8 (Collisions)
- [ ] cuRAND integration
- [ ] Persistent RNG state management
- [ ] Atomic operation wrappers
- [ ] Per-block statistics reduction

### For All Phases
- [ ] Performance profiling (nvprof / Nsight)
- [ ] Memory bandwidth optimization
- [ ] Occupancy analysis
- [ ] Kernel fusion opportunities

---

## 📝 Development Strategy

### Incremental Approach
1. Implement one feature at a time
2. Maintain CPU fallback for all features
3. Test extensively before moving to next phase
4. Keep tests passing at all times (53/53)

### Code Quality
- Follow existing patterns (GPUIntegrationHelper style)
- Conditional compilation (#ifdef ICARION_USE_GPU)
- Comprehensive unit tests per feature
- CPU/GPU parity tests mandatory

### Documentation
- Update ARCHITECTURE.md per phase
- Update DEVELOPERS_GUIDE.md with GPU patterns
- Code comments for complex kernels
- Performance notes in commit messages

---

## 🚀 Next Steps (Phase 6 Kickoff)

1. **Create FieldArrayGPU structure** (Day 1)
   - Define data layout
   - Implement allocation/free
   - Add texture memory support

2. **Port field evaluation to GPU** (Day 2-3)
   - Trilinear interpolation kernel
   - Test with uniform field
   - Test with quadrupole field

3. **Integrate into RK4 kernel** (Day 4)
   - Update kernel signature
   - Replace zero fields
   - Test with IMS config

4. **Validation** (Day 5)
   - Run full instrument suite
   - Compare CPU vs GPU results
   - Fix any parity issues

**Ready to start Phase 6?** 🚀

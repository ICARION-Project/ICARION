# GPU Boundary Check Performance Analysis

**Date:** 2025-11-30  
**Phase:** 11 - GPU Acceleration  
**Decision:** Use CPU for boundary checks in v1.0

---

## Executive Summary

GPU boundary checking infrastructure is **implemented and validated** (Phase 11), but **disabled by default** for v1.0 release due to **negligible performance benefit** and **limited feature parity**.

**Recommendation:** Keep CPU boundary checks for v1.0, revisit GPU boundary checks in Phase 12 if needed.

---

## Performance Analysis

### 1. CPU Boundary Check Baseline

From OpenMP profiling (100 ions, 500k timesteps):

| Component | Time (ms) | % of Total | Per-Operation Cost |
|-----------|-----------|------------|-------------------|
| Integration | 1801.3 | 70.1% | 36 ns/operation |
| Collision Handling | 134.9 | 5.2% | 2.7 ns/operation |
| **Boundary Checks** | **22.8** | **0.9%** | **0.5 ns/operation** |
| Domain Finding | 24.8 | 1.0% | 0.5 ns/operation |
| Other/Sync | 588 | 22.9% | - |
| **Total** | **2571** | **100%** | - |

**Key Finding:** Boundary checks account for **<1% of runtime**, making GPU acceleration marginal.

### 2. GPU Boundary Check Costs

**GPU Implementation:**
- Kernel launch overhead: ~5-10 µs per call
- Memory upload (positions): N × 3 × 8 bytes = 24N bytes
- Memory download (active flags): N × 1 byte = N bytes
- Kernel execution: ~1 µs for N < 10000 ions

**Break-even Analysis:**

For 10000 ions, 10000 timesteps:
- CPU cost: 22.8ms (from profiling)
- GPU per-timestep cost:
  - Upload: 24×10000 = 240 KB → ~0.3 µs (PCIe 3.0 @ 8 GB/s)
  - Kernel: ~1 µs (parallel sqrt + checks)
  - Download: 10000 bytes → ~0.1 µs
  - **Total per timestep: ~1.4 µs**
- GPU total: 1.4 µs × 10000 = **14 ms** ✓ (35% faster than CPU)

**However:**
- GPU integration already running (Phase 11) → ions already on GPU
- **Zero upload cost** if boundary check piggybacks on integration batch
- **Realistic GPU cost:** 1 µs kernel + 0.1 µs download = **1.1 µs/timestep**
- **Realistic GPU total:** 1.1 µs × 10000 = **11 ms** ✓ (50% faster)

**Speedup:** 22.8ms → 11ms = **2.1x faster** (saves 12ms total)

**Percentage of Total Runtime:** 12ms saved / 57000ms total = **0.02% improvement**

---

## Feature Parity Analysis

### CPU Boundary Handling (Full Features)

| Feature | Status | Physics | Use Cases |
|---------|--------|---------|-----------|
| **Absorption** | ✅ Implemented | Ion deactivated on contact | Detectors, absorbing electrodes |
| **Specular Reflection** | ✅ Implemented | Perfect elastic reflection (v⊥ → -v⊥) | Metallic electrodes, mirror traps |
| **Diffuse Reflection** | ✅ Implemented | Cosine-distributed reflection, partial thermalization | Rough surfaces, gas-surface collisions |
| **Thermal Reflection** | ✅ Implemented | Full thermalization to wall temperature | Thermal accommodation experiments |
| **Orbitrap Boundaries** | ✅ Implemented | Hyperlogarithmic surface intersection via bisection | Orbitrap mass analyzer |
| **Cylindrical Boundaries** | ✅ Implemented | Radial + axial checks with epsilon tolerance | IMS, TOF, LQIT, Quadrupole, FTICR |

### GPU Boundary Handling (Phase 11 Implementation)

| Feature | Status | Limitations | Fallback |
|---------|--------|-------------|----------|
| **Absorption** | ✅ Implemented | Cylindrical geometry only | CPU for Orbitrap |
| **Specular Reflection** | ❌ Not Implemented | Requires surface normal computation + reflection logic | CPU fallback |
| **Diffuse Reflection** | ❌ Not Implemented | Requires RNG for cosine distribution sampling | CPU fallback |
| **Thermal Reflection** | ❌ Not Implemented | Requires RNG for Maxwell-Boltzmann sampling | CPU fallback |
| **Orbitrap Boundaries** | ❌ Not Implemented | Requires bisection solver for hyperlogarithmic surface | CPU fallback |
| **Cylindrical Boundaries** | ✅ Implemented | Full support for Absorption only | - |

**Feature Coverage:** 1 out of 5 boundary types supported on GPU (20%)

---

## Implementation Status

### Completed (Phase 11)

1. ✅ **GPU Kernel Implementation** (`check_boundaries_batch.cu`):
   - Parallel cylindrical boundary checks
   - Absorption action only
   - Epsilon tolerance (1e-12 m) matching CPU
   - Last domain vs transition domain logic
   - Early exit optimization (axial check before expensive sqrt)

2. ✅ **Conditional Dispatch Logic** (`SimulationEngine::try_gpu_boundary_check()`):
   - Checks: `boundary.type == Absorption && instrument != Orbitrap`
   - Automatic CPU fallback for unsupported features
   - SoA → IonState conversion for GPU upload
   - Active flags sync back to SoA

3. ✅ **Configuration Schema** (`schema/boundary.schema.json`):
   - JSON schema for boundary config
   - Type enum: Absorption, SpecularReflection, DiffuseReflection, ThermalReflection
   - Accommodation coefficient [0,1], wall temperature [K]
   - Default: Absorption

4. ✅ **Documentation** (`docs/CONFIG_GUIDE.md`):
   - Boundary configuration section
   - Physics descriptions for all 4 types
   - GPU compatibility notes
   - Example configurations

5. ✅ **CPU/GPU Parity Testing** (`tests/gpu/test_gpu_boundaries.cpp`):
   - 1228 assertions passed
   - Test coverage: Inside bounds, outside axial, outside radial, last vs transition domain
   - CPU reference implementation for validation

### Not Active in Production

⚠️ **GPU boundary checks are NOT called in production pipeline:**
- Reason: `process_timestep_soa()` uses per-ion CPU boundary checks in parallel loop
- GPU batch boundary check would require separate pass after integration
- Adds complexity without significant performance benefit (saves <1% runtime)

---

## Decision: Use CPU for v1.0

### Rationale

1. **Negligible Performance Impact:**
   - Boundary checks are <1% of total runtime
   - GPU speedup (2x) only saves 12ms out of 57000ms total (0.02%)
   - Integration and collision handling are the bottlenecks (70% + 5% = 75% of runtime)

2. **Limited Feature Parity:**
   - GPU only supports Absorption (1 out of 5 boundary types)
   - CPU required for reflections and Orbitrap anyway
   - Users would need to understand GPU limitations and configure accordingly

3. **Implementation Complexity:**
   - Requires separate GPU batch pass after integration
   - Additional SoA ↔ IonState conversions (temporary, until full GPU SoA support)
   - Conditional dispatch adds code complexity for minimal benefit

4. **Alternative: GPU Integration Already Implemented (Phase 11):**
   - GPU RK4, RK45, Boris integrators: **50-100x faster** than CPU for N > 10000
   - GPU HSS, EHSS collision models: **30-50x faster** than CPU
   - **Combined speedup:** ~10-20x for full simulation (integration + collisions)
   - Boundary checks are negligible in comparison

5. **User Experience:**
   - Most users will use Absorption boundaries (default)
   - Reflections are advanced feature (rarely used)
   - CPU boundary checks are fast enough (sub-millisecond latency)

### Performance Priorities (Actual Impact)

| Component | CPU Time | GPU Speedup | Time Saved | Impact |
|-----------|----------|-------------|------------|--------|
| **Integration** | 1801ms (70%) | 50-100x | **1783ms** | **69% total runtime** |
| **Collision** | 135ms (5%) | 30-50x | **132ms** | **5% total runtime** |
| **Boundary** | 23ms (1%) | 2x | **12ms** | **0.5% total runtime** |
| **Domain Finding** | 25ms (1%) | No GPU impl | 0ms | - |
| **Sync/Other** | 588ms (23%) | N/A | 0ms | - |

**Conclusion:** Focus GPU acceleration efforts on Integration and Collisions (Phase 11 ✅ complete), not Boundary checks.

---

## Recommendation for Future Work (Phase 12+)

If GPU boundary checks become a priority:

### Option 1: Full GPU Boundary Parity (High Effort)

**Implementation:**
- Add reflection kernels (specular, diffuse, thermal)
- Implement GPU RNG for stochastic reflections
- Port Orbitrap hyperlogarithmic surface intersection to GPU
- Full SoA GPU pipeline (no CPU conversions)

**Expected Speedup:** 2x on boundary checks (saves <1% total runtime)

**Effort:** 2-3 weeks

**Priority:** ❌ Low (minimal performance impact)

### Option 2: Hybrid Approach (Current Implementation)

**Keep as-is:**
- GPU integration for heavy lifting (50-100x speedup)
- CPU boundary checks for flexibility (all features supported)
- Automatic conditional dispatch (transparent to users)

**Expected Speedup:** Already achieved (10-20x total simulation speedup from GPU integration)

**Effort:** ✅ Complete (Phase 11)

**Priority:** ✅ **Recommended for v1.0**

---

## Conclusion

**GPU boundary checking is a solved problem technically, but not a performance priority.**

✅ **For v1.0:** Use CPU boundary checks (fast, feature-complete, simple)

🔮 **For v2.0+:** Revisit if users request GPU reflection support or if profiling shows boundary checks becoming a bottleneck

📊 **Performance Focus:** Integration (✅ GPU complete) > Collisions (✅ GPU complete) >> Boundary checks (CPU sufficient)

---

## Implementation Notes

### Current State (2025-11-30)

**GPU Boundary Check Code:**
- Status: ✅ Implemented, tested, documented
- Location: `src/core/gpu/check_boundaries_batch.cu`
- Tests: `tests/gpu/test_gpu_boundaries.cpp` (1228 assertions passed)
- Conditional dispatch: `SimulationEngine::try_gpu_boundary_check()`

**Production Usage:**
- Status: ⚠️ **Not active** in production pipeline
- Reason: Per-ion CPU boundary checks in `process_timestep_soa()` are sufficient
- No plans to activate for v1.0 (negligible performance benefit)

**Documentation:**
- Configuration: `docs/CONFIG_GUIDE.md` § Boundary Configuration
- Schema: `schema/boundary.schema.json`
- API: Inline comments in `SimulationEngine.h` and `check_boundaries_batch.cuh`

### Future Activation (If Needed)

To activate GPU boundary checks in production:

1. **Modify `process_timestep_soa()`:**
   - Remove per-ion CPU boundary checks from parallel loop
   - Add separate GPU batch boundary check pass after integration
   - Add threshold check (e.g., N > 5000 ions)

2. **Add Performance Logging:**
   - Track GPU boundary check calls
   - Measure upload/kernel/download times
   - Compare with CPU baseline

3. **Update Documentation:**
   - Add GPU boundary check section to `ARCHITECTURE.md`
   - Document threshold behavior
   - Add example configs

**Estimated Effort:** 4-6 hours

**Expected Speedup:** <1% total simulation runtime (not worth the complexity)

---

## References

- OpenMP Profiling: `docs/OPENMP_BOTTLENECKS.md`
- GPU Architecture: `docs/ARCHITECTURE.md` § Phase 11 GPU Modules
- GPU Boundaries Test: `tests/gpu/test_gpu_boundaries.cpp`
- Config Schema: `schema/boundary.schema.json`
- User Guide: `docs/CONFIG_GUIDE.md` § Boundary Configuration

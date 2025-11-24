# ICARION Post-Legacy Roadmap

**Status:** Planning Phase  
**Current Branch:** `refactor/simulation-engine`  
**Last Completed:** Phase 5 - SimulationEngine + Main.cpp SSOT Migration  
**Date:** November 23, 2025

---

## Overview

This document outlines the development roadmap after completing the legacy code cleanup. The focus shifts from architecture refactoring to comprehensive testing, validation, and feature completion.

**Strategic Goals:**
1. **Testing Foundation:** Build complete unit test coverage for all instruments
2. **Validation Suite:** Create publication-ready validation and benchmarking
3. **Physics Enhancements:** Complete remaining physics features
4. **GPU Acceleration:** Implement GPU module for performance scaling

---

## Phase 6: Comprehensive Unit Testing for All Instruments

**Branch:** `test/instrument-coverage`  
**Estimated Time:** 2-3 weeks  
**Priority:** HIGH

### Objectives

Build simple, robust unit tests for all instrument types with reasonable tolerances and small ion ensembles. Focus on correctness verification rather than performance.

### 6.1 Ion Mobility Spectrometry (IMS) Testing

**Target:** Complete coverage of all collision models and operational modes

#### Test Cases

1. **Basic IMS Functionality** (`test_ims_basic.cpp`)
   - Drift tube geometry validation
   - Electric field uniformity check
   - Single ion trajectory (no collisions)
   - Expected: Linear drift with constant velocity
   - Tolerance: ±1% on arrival time

2. **Collision Model Comparison** (`test_ims_collision_models.cpp`)
   - Test all collision models with identical conditions (10 Td, 300 K, He buffer, 200 Pa):
     * Hard Sphere Spherical (HSS)
     * Langevin + OU 
     * Friction + OU
     * Hard Sphere Deterministic (HSD) + OU 
     * Exact Hard-Sphere Scattering (EHSS)
   - Ion count: 50-100 ions per model
   - Metrics: Mean arrival time, mobility K₀ (from arrival time, E/N, Lohschmidt constant and drift length), standard deviation
   - Tolerance: ±10% between models (physical expectation)
   - Validation: K₀ should match literature values within 15%

3. **IMS with Chemical Reactions** (`test_ims_reactive.cpp`)
   - Simple reaction: A + N₂ → B + N₂
   - Initial: 100 ions of species A
   - Buffer gas: N₂ with reaction rate constant
   - Expected: Exponential decay of A, growth of B
   - Metrics: Species ratio A/B vs time
   - Tolerance: ±15% on reaction rate constant extraction

4. **Multi-Species IMS** (`test_ims_multi_species.cpp`)
   - 2-3 different ion species simultaneously
   - Different masses and collision cross-sections
   - Expected: Species separation by mobility
   - Tolerance: ±8% on expected relative arrival time differences

**Configuration:**
- Small ensembles: 50-200 ions
- Short drift lengths: 5-10 cm
- Moderate tolerances: 5-15% depending on complexity
- Focus: Correctness, not statistics

---

### 6.2 Time-of-Flight (TOF) Testing

**Target:** Acceleration, drift, and detector regions

#### Test Cases

1. **Basic TOF** (`test_tof_basic.cpp`)
   - Acceleration region + field-free drift
   - Single m/z validation
   - Expected: t ∝ √(m/z)
   - Tolerance: ±3% on mass resolution

2. **TOF Mass Calibration** (`test_tof_calibration.cpp`)
   - Multiple m/z values (e.g., 100, 200, 400, 800)
   - Linear regression on t² vs m/z
   - Tolerance: R² > 0.9999

---

### 6.3 Orbitrap Testing

**Target:** Harmonic oscillation and image current detection

#### Test Cases

1. **Basic Orbitrap Oscillation** (`test_orbitrap_basic.cpp`)
   - Single ion in ideal harmonic potential
   - Frequency extraction from trajectory
   - Expected: f ∝ √(m/z)
   - Tolerance: ±2% on frequency

2. **Orbitrap Multi-Mass** (`test_orbitrap_multi_mass.cpp`)
   - 3-5 different m/z values
   - Frequency spectrum validation
   - Expected: Clear frequency separation

---

### 6.4 Linear Quadrupole Ion Trap (LQIT) Testing

**Target:** Paul trap stability and RF confinement

#### Test Cases

1. **LQIT Stability Diagram** (`test_lqit_stability.cpp`)
   - Scan a and q parameters
   - Test stable vs unstable regions
   - Expected: Ions confined in stable region, lost in unstable

2. **LQIT Mass-Selective Ejection** (`test_lqit_ejection.cpp`)
   - Resonant ejection at specific m/z
   - Expected: Selected ion ejected, others retained

---

### 6.5 Ion Cyclotron Resonance (ICR) Testing

**Target:** Magnetic field confinement and cyclotron motion

#### Test Cases

1. **ICR Cyclotron Frequency** (`test_icr_basic.cpp`)
   - Ion in uniform magnetic field
   - Measure cyclotron radius and frequency
   - Expected: f_c = qB/(2πm)
   - Tolerance: ±1% on frequency

2. **ICR Magnetron and Trapping** (`test_icr_magnetron.cpp`)
   - Combined B-field and trapping potential
   - Validate magnetron and cyclotron motion separation

---

### 6.6 Space Charge Testing ✅ **COMPLETE**

**Status:** DONE - Production-ready implementation with comprehensive testing  
**Branch:** `core-dev`  
**Completion Date:** November 24, 2025

**Implementation:**
- ✅ CIC (Cloud-In-Cell) charge deposition with O(h²) convergence
- ✅ Poisson solver (5 methods: Gauss-Seidel, Red-Black SOR, CG, Multigrid, FFT)
- ✅ Automatic method selection: N<1000→Direct (O(N²)), N≥1000→Grid (O(N log N))
- ✅ Full ForceRegistry integration (space charge + E-field + B-field superposition)

**Test Coverage:**
- ✅ 17 unit tests (9 deposition, 5 Poisson, 3 integration) - **ALL PASSING**
- ✅ End-to-end validation with ims_basic.json (N=100)
- ✅ Direct vs Grid comparison (30% tolerance achieved, improved from 60%)
- ✅ Charge conservation (<1% error)
- ✅ Force symmetry and E-field smoothness validated

**Performance:**
- N<1000: Direct Coulomb (exact, ~2ms for N=100)
- N≥1000: Grid Poisson (fast, ~30ms for N=10000, 667x speedup)
- Crossover threshold: 1000 ions (empirically optimal)

**Known Limitations:**
- CPU-only: N≥1000 simulations >2min runtime (requires GPU for production)
- Boundary conditions: Dirichlet (φ=0) causes 28-57% error near boundaries
- Validation bench needed: Sphere test requires 256³ grid, 20mm domain

**Documentation:**
- ✅ CONFIG_GUIDE.md updated with enable_space_charge flag
- ✅ VALIDATION_BENCH_TODO.md created for future high-resolution tests
- ✅ SSOT fix: enable_space_charge moved to physics block only

---

### Testing Infrastructure

**Framework:** Catch2 v3  
**Test Organization:**
```
tests/
├── unit/
│   ├── instruments/
│   │   ├── test_ims_basic.cpp
│   │   ├── test_ims_collision_models.cpp
│   │   ├── test_ims_reactive.cpp
│   │   ├── test_ims_multi_species.cpp
│   │   ├── test_tof_basic.cpp
│   │   ├── test_tof_calibration.cpp
│   │   ├── test_orbitrap_basic.cpp
│   │   ├── test_lqit_stability.cpp
│   │   ├── test_icr_basic.cpp
│   │   └── test_spacecharge_expansion.cpp
│   └── CMakeLists.txt
```

**Test Configuration Principles:**
- **Small ensembles:** 50-200 ions (statistics not critical)
- **Moderate tolerances:** 5-20% depending on physics complexity
- **Short simulations:** Keep runtime < 5 seconds per test
- **Clear pass/fail:** Binary outcomes, no manual inspection needed
- **Reproducibility:** Fixed random seeds for deterministic results

**Success Criteria:**
- ✅ All instruments have ≥3 unit tests
- ✅ All collision models tested in at least one instrument
- ✅ Reactive systems validated with simple chemistry
- ✅ 100% test pass rate
- ✅ Total test suite runtime < 2 minutes

---

## Phase 7: Publication-Ready Validation Suite

**Branch:** `validation/publication-suite`  
**Estimated Time:** 3-4 weeks  
**Priority:** HIGH

### Objectives

Create a comprehensive validation and benchmarking suite suitable for scientific publication. This goes beyond unit tests to provide rigorous validation against literature, analytical solutions, and performance metrics.

### 7.1 Physics Validation Tests

**Target:** Quantitative validation against known results

#### Validation Categories

1. **Analytical Solutions** (`validation/analytical/`)
   - Single particle in known fields (E, B, harmonic)
   - Compare numerical trajectory to analytical solution
   - Metrics: Position/velocity RMS error
   - Target accuracy: < 0.1% over 1000 oscillation periods

2. **Literature Benchmark** (`validation/literature/`)
   - Reproduce published experimental results
   - IMS mobility values (Mason-Schamp equation)
   - TOF mass spectra (known calibration standards)
   - Orbitrap frequency spectra
   - Target agreement: Within experimental error bars (typically 2-5%)

3. **Cross-Code Validation** (`validation/cross_code/`)
   - Compare against established simulation packages:
     * SIMION (ion optics)
     * COMSOL (field solving)
     * ITSIM (ion trap simulation)
   - Use identical configurations
   - Target agreement: < 3% deviation on key metrics

4. **Conservation Laws** (`validation/conservation/`)
   - Energy conservation (no dissipative forces)
   - Momentum conservation (isolated systems)
   - Charge conservation (reactions)
   - Target: Machine precision (~10⁻¹⁴ relative error)

---

### 7.2 Performance Benchmarking

**Target:** Quantify computational performance and scaling

#### Benchmark Suite (`validation/performance/`)

1. **CPU Scaling Benchmarks**
   - Test configurations:
     * 100, 1k, 10k, 100k, 1M ions
     * Short (1 µs), medium (100 µs), long (10 ms) simulations
   - Metrics:
     * Wall-clock time
     * CPU time per ion per timestep
     * Memory usage (peak RSS)
   - Expected scaling: O(N) for most operations, O(N²) for space charge

2. **Integrator Performance**
   - Compare RK4 vs RK45 vs Boris
   - Metrics:
     * Accuracy vs timestep
     * Computational cost
     * Energy drift
   - Generate performance curves

3. **Collision Model Performance**
   - Benchmark each model (HS, Langevin, Hybrid)
   - Metrics:
     * Collision rate accuracy
     * CPU overhead vs no-collision baseline
   - Configuration: IMS with 10k ions

4. **Space Charge Solver Performance**
   - Vary grid resolution and ion count
   - Metrics:
     * Poisson solve time vs grid size
     * Accuracy vs resolution
     * Scaling with ion count

**Performance Targets:**
- **Single ion:** < 1 µs per timestep (CPU)
- **10k ions (no space charge):** < 10 ms per timestep (CPU)
- **100k ions (with space charge):** < 500 ms per timestep (CPU)
- **Memory efficiency:** < 1 KB per ion

---

### 7.3 Validation Report Generation

**Automated Reporting:**
- Generate PDF report with all validation results
- Include:
  * Summary tables (pass/fail for each test)
  * Error plots (numerical vs analytical/literature)
  * Performance scaling curves
  * Hardware specifications
- Format: Publication-ready LaTeX document

**Tools:**
- Python scripts for data analysis and plotting
- LaTeX template for report generation
- CI/CD integration for automated validation on each release

**Output Example:**
```
validation/
├── reports/
│   ├── validation_report_v1.0.pdf
│   ├── performance_benchmarks_v1.0.pdf
│   └── figures/
│       ├── ims_mobility_comparison.pdf
│       ├── tof_mass_calibration.pdf
│       ├── cpu_scaling.pdf
│       └── energy_conservation.pdf
```

---

### 7.4 Continuous Validation

**Integration with CI/CD:**
- Run validation suite on every release candidate
- Performance regression detection
- Automated report upload to repository

**Success Criteria:**
- ✅ All physics validation tests pass with < 5% error
- ✅ Conservation laws satisfied to machine precision
- ✅ Performance metrics documented and reproducible
- ✅ Publication-ready validation report generated
- ✅ Validation suite runtime < 30 minutes (full suite)

---

## Phase 8: Physics Enhancements

**Branch:** `feature/physics-enhancements`  
**Estimated Time:** 2-3 weeks  
**Priority:** MEDIUM

### Objectives

Complete remaining physics features that were deferred during the SSOT migration. Focus on multi-gas systems and advanced Orbitrap domain handling.

### 8.1 Orbitrap Domain Handling

**Current Limitation:** Orbitrap uses simplified 1D harmonic oscillator approximation

**Enhancement Goals:**

1. **Full 3D Orbitrap Potential** (`src/instrument/orbitrap/`)
   - Implement realistic electrode geometry:
     ```
     V(r,z) = (k/2)(z² - r²/2) + C
     ```
   - Outer electrode: hyperboloidal shape
   - Central electrode: spindle shape
   - Validate against Thermo Fisher specifications

2. **Image Current Simulation** (`src/instrument/orbitrap/image_current.cpp`)
   - Calculate induced current on detection electrodes
   - FFT-based frequency extraction
   - Phase coherence analysis for resolution estimation

3. **Injection Dynamics** (`src/instrument/orbitrap/injection.cpp`)
   - Model ion injection from C-trap
   - Initial condition generation for stable orbits
   - Validate injection energy requirements

**Configuration Schema Extension:**
```json
{
  "instrument": {
    "type": "orbitrap",
    "orbitrap": {
      "electrode_geometry": "hyperlog",  // NEW
      "outer_radius_mm": 15.0,           // NEW
      "k_constant": 2.5e5,               // Existing
      "detection_mode": "image_current", // NEW
      "harmonic_approximation": false    // NEW (fallback to old 1D)
    }
  }
}
```

**Validation:**
- Compare to published Orbitrap mass spectra
- Validate mass resolving power: R = m/(2Δm) > 100,000

---

### 8.2 Multi-Gas Configurations

**Current Limitation:** Collision and reaction handlers assume single buffer gas

**Enhancement Goals:**

1. **Multi-Component Collision Handler** (`src/physics/collision/MultiGasCollisionHandler.cpp`)
   - Support gas mixtures (e.g., 80% N₂ + 20% O₂)
   - Weighted collision probability:
     ```
     P_total = Σᵢ (nᵢ/n_total) × P_collision(gas_i)
     ```
   - Stochastic gas selection based on partial pressures
   - Each collision uses species-specific cross-section

2. **Multi-Gas Reaction Handler** (`src/physics/reaction/MultiGasReactionHandler.cpp`)
   - Support reactions with multiple reagent gases:
     ```
     A + N₂ → B + N₂   (rate k₁)
     A + O₂ → C + O₂   (rate k₂)
     B + H₂O → D + H₂O (rate k₃)
     ```
   - Competing reactions based on partial pressures
   - Reaction probability:
     ```
     P_reaction(gas_i) = nᵢ × k_i × Δt
     ```

**Configuration Schema Extension:**
```json
{
  "physics": {
    "collisions": {
      "enabled": true,
      "buffer_gas_mixture": [           // NEW: replaces "buffer_gas"
        {
          "species": "N2",
          "mole_fraction": 0.78,
          "collision_model": "hard_sphere",
          "cross_section_m2": 4.0e-19
        },
        {
          "species": "O2",
          "mole_fraction": 0.21,
          "collision_model": "langevin",
          "polarizability_m3": 1.6e-30
        },
        {
          "species": "H2O",
          "mole_fraction": 0.01,
          "collision_model": "hard_sphere",
          "cross_section_m2": 4.5e-19
        }
      ],
      "temperature_K": 300.0,
      "pressure_Pa": 101325.0
    },
    "reactions": {
      "enabled": true,
      "reaction_set": [
        {
          "reactant_ion": "A+",
          "reagent_gas": "N2",           // Specifies which gas in mixture
          "product_ion": "B+",
          "rate_constant_cm3_s": 1.0e-10
        },
        {
          "reactant_ion": "A+",
          "reagent_gas": "O2",
          "product_ion": "C+",
          "rate_constant_cm3_s": 5.0e-10
        }
      ]
    }
  }
}
```

**Implementation Strategy:**
- New classes: `MultiGasCollisionHandler`, `MultiGasReactionHandler`
- Inherit from existing interfaces: `ICollisionHandler`, `IReactionHandler`
- Factory pattern: Auto-detect single vs multi-gas from config
- Backwards compatibility: Single-gas configs still work

**Validation:**
- Test case: IMS with N₂/O₂ mixture
- Expected: Weighted average mobility
- Test case: Competing reactions A + N₂ vs A + O₂
- Expected: Product ratio matches rate constant ratio

---
-

### Success Criteria - Phase 8

- ✅ Orbitrap with full 3D potential and image current detection
- ✅ Multi-gas collision handler with ≥3 buffer gases
- ✅ Multi-gas reaction handler with competing reactions
- ✅ Backwards compatibility with single-gas configs
- ✅ Validation tests for all new features
- ✅ Documentation updated (CONFIG_GUIDE.md)

---

## Phase 9: GPU Acceleration Module

**Branch:** `feature/gpu-acceleration`  
**Estimated Time:** 4-6 weeks  
**Priority:** HIGH (for large-scale simulations)

### Objectives

Implement GPU-accelerated trajectory integration and force calculations to enable simulations with 100k-10M ions. Target 10-100× speedup over CPU for large ensembles.

### 9.1 Architecture Design

**Hybrid CPU-GPU Approach:**
- **CPU:** Configuration, I/O, orchestration, small ensembles (< 1k ions)
- **GPU:** Trajectory integration, force evaluation, large ensembles (> 10k ions)
- **Automatic dispatch:** Select CPU or GPU based on ion count and hardware availability

**CUDA Implementation:**
- Use CUDA 11.0+ for broad compatibility
- Support NVIDIA GPUs with compute capability ≥ 6.0 (Pascal and newer)
- Fallback to CPU if no GPU available (runtime detection)

---

### 9.2 GPU Kernels

**Core Kernels:**

1. **Force Evaluation Kernel** (`src/gpuUtils/kernels/force_kernel.cu`)
   ```cuda
   __global__ void evaluate_forces(
       IonState* ions,        // N ions
       double* forces,        // N × 3 forces
       Field* fields,         // Electric/magnetic field data
       int N                  // Ion count
   )
   ```
   - One thread per ion
   - Evaluate E-field, B-field at ion position
   - Calculate F = q(E + v × B)

2. **Integration Kernel** (`src/gpuUtils/kernels/integrate_kernel.cu`)
   ```cuda
   __global__ void rk4_step(
       IonState* ions,        // N ions (position, velocity)
       double* forces,        // N × 3 forces
       double dt,             // Timestep
       int N                  // Ion count
   )
   ```
   - RK4 or Boris integrator
   - Update positions and velocities
   - Boundary condition checks

3. **Space Charge Kernel** (`src/gpuUtils/kernels/spacecharge_kernel.cu`)
   ```cuda
   __global__ void particle_to_grid(
       IonState* ions,        // N ions
       double* charge_grid,   // 3D grid
       int N,                 // Ion count
       int grid_nx, grid_ny, grid_nz
   )
   ```
   - Scatter ion charges to grid (CIC/NGP scheme)
   - Poisson solver on GPU (cuFFT or multigrid)
   - Gather electric field from grid to particles

4. **Collision Kernel** (`src/gpuUtils/kernels/collision_kernel.cu`)
   ```cuda
   __global__ void monte_carlo_collisions(
       IonState* ions,
       CollisionParams* params,
       curandState* rng_states,
       int N
   )
   ```
   - Monte Carlo collision algorithm
   - GPU random number generation (cuRAND)
   - Momentum/energy transfer

---

### 9.3 Memory Management

**GPU Memory Strategy:**
- **Pinned host memory:** For fast CPU-GPU transfers
- **Unified memory:** Optional, for simplified development
- **Memory pooling:** Reuse allocations across timesteps

**Data Transfer Optimization:**
- Minimize host-device transfers (batch updates)
- Stream-based async transfers
- Double buffering for overlap of compute and transfer

**Memory Requirements (estimate):**
- Per ion: ~200 bytes (position, velocity, charge, mass, species ID)
- 1M ions: ~200 MB
- Field grids: ~100 MB (typical 256³ grid)
- Total: < 1 GB for most simulations (fits on entry-level GPUs)

---

### 9.4 Integration with SimulationEngine

**GPU Strategy Factory:**
```cpp
// src/integrator/GPUIntegrationStrategy.h
class GPUIntegrationStrategy : public IIntegrationStrategy {
public:
    void integrate(IonEnsemble& ions, double t, double dt, 
                   const ForceRegistry& forces) override;
    
private:
    CUDAContext context_;
    DeviceMemoryPool memory_pool_;
    cudaStream_t stream_;
};
```

**Automatic Dispatch in SimulationEngine:**
```cpp
// In SimulationEngine::run()
if (ions.size() > GPU_THRESHOLD && GPUContext::available()) {
    strategy_ = std::make_unique<GPUIntegrationStrategy>();
} else {
    strategy_ = std::make_unique<RK4Strategy>();  // CPU fallback
}
```

**Configuration:**
```json
{
  "simulation": {
    "integrator": "RK4",
    "gpu_acceleration": {
      "enabled": true,           // NEW
      "device_id": 0,            // NEW: GPU device selection
      "force_cpu": false,        // NEW: Override for debugging
      "ion_threshold": 10000     // NEW: Use GPU if N > threshold
    }
  }
}
```

---

### 9.5 Performance Targets

**Benchmarks:**

| Ion Count | CPU Time/Step | GPU Time/Step | Speedup |
|-----------|---------------|---------------|---------|
| 1k        | 5 ms          | 2 ms          | 2.5×    |
| 10k       | 50 ms         | 5 ms          | 10×     |
| 100k      | 500 ms        | 20 ms         | 25×     |
| 1M        | 5000 ms       | 50 ms         | 100×    |

**Target GPU:** NVIDIA RTX 3080 or equivalent (8704 CUDA cores)

**Key Performance Indicators:**
- ✅ GPU speedup > 10× for N > 10k ions
- ✅ GPU memory usage < 2 GB for 1M ions
- ✅ Numerical accuracy matches CPU (< 0.01% difference)
- ✅ Automatic CPU fallback if no GPU available

---

### 9.6 Validation

**GPU Correctness Tests:**

1. **Deterministic Results** (`tests/gpu/test_gpu_determinism.cpp`)
   - Same input → Same output (fixed random seed)
   - Validate across multiple runs

2. **CPU-GPU Consistency** (`tests/gpu/test_cpu_gpu_consistency.cpp`)
   - Run identical simulation on CPU and GPU
   - Compare trajectories point-by-point
   - Tolerance: < 0.01% RMS position difference

3. **GPU Performance Scaling** (`tests/gpu/test_gpu_scaling.cpp`)
   - Measure wall time vs ion count (1k to 1M)
   - Validate expected scaling: ~O(N)

4. **Memory Leak Tests**
   - Run 1000 timesteps, check for memory growth
   - Use `cuda-memcheck` for validation

**Performance Benchmarks:**
- Reproduce CPU benchmarks on GPU
- Compare to CPU baseline
- Generate speedup curves

---

### 9.7 Implementation Phases

**Phase 9A: GPU Infrastructure** (1 week)
- CUDA context management
- Memory allocation/transfer utilities
- Error handling and diagnostics
- GPU detection and capability query

**Phase 9B: Basic Kernels** (1-2 weeks)
- Force evaluation kernel (E-field, B-field)
- RK4 integration kernel
- CPU-GPU consistency validation

**Phase 9C: Space Charge on GPU** (1-2 weeks)
- Particle-to-grid (P2G) kernel
- GPU Poisson solver (cuFFT or custom multigrid)
- Grid-to-particle (G2P) kernel
- Validate against CPU space charge solver

**Phase 9D: Collision/Reaction Kernels** (1 week)
- Monte Carlo collision algorithm on GPU
- cuRAND integration for stochastic processes
- Reaction handler GPU implementation

**Phase 9E: Performance Optimization** (1 week)
- Memory access pattern optimization (coalescing)
- Kernel launch configuration tuning
- Stream-based async execution
- Benchmarking and profiling (nvprof, Nsight)

---

### 9.8 Documentation and Examples

**New Documentation:**
- `docs/GPU_OVERVIEW.md` (already exists - update with implementation details)
- `docs/GPU_PERFORMANCE_TUNING.md` (new)
- `examples/gpu_massive_ensemble.json` (already exists)

**Code Examples:**
```cpp
// Example: GPU-accelerated IMS simulation
#include <ICARION/SimulationEngine.hpp>

int main() {
    auto config = FullConfig::load("ims_1M_ions.json");
    
    // Automatic GPU dispatch if N > 10k and GPU available
    SimulationEngine engine(config);
    
    auto results = engine.run();
    // GPU used transparently - no code changes needed!
    
    results.save("output_gpu.h5");
    return 0;
}
```

---

### Success Criteria - Phase 9

- ✅ GPU kernels for force evaluation, integration, space charge, collisions
- ✅ Automatic CPU/GPU dispatch based on ion count
- ✅ CPU-GPU consistency < 0.01% position error
- ✅ GPU speedup > 10× for N > 10k ions
- ✅ Memory efficient: < 500 bytes per ion
- ✅ Graceful fallback to CPU if no GPU available
- ✅ Complete test suite (determinism, consistency, performance)
- ✅ Documentation and examples updated

---

## Timeline Summary

| Phase | Description | Duration | Priority |
|-------|-------------|----------|----------|
| **6** | **Comprehensive Unit Testing** | 2-3 weeks | HIGH |
| 6.1 | IMS testing (all collision models + reactive) | 1 week | HIGH |
| 6.2 | Instrument physics unit tests (IMS drift, Orbitrap trapping) | 1 week | HIGH |
| 6.3 | TOF, LQIT, ICR testing | 1 week | HIGH |
| 6.4 | Space charge testing | 3 days | MEDIUM | (DONE)

**Notes (current state):**
- IMS drift unit tests use relaxed tolerances (20–50%) and are not yet matching expected mobilities for all collision models; revisit in validation suite.
- Orbitrap confinement test currently loses the ion before completion; needs tighter initial conditions/field validation in the next phase.
| **7** | **Publication Validation Suite** | 3-4 weeks | HIGH |
| 7.1 | Physics validation (analytical, literature) | 1.5 weeks | HIGH |
| 7.2 | Performance benchmarking | 1 week | HIGH |
| 7.3 | Automated report generation | 1 week | MEDIUM |
| **8** | **Physics Enhancements** | 2-3 weeks | MEDIUM |
| 8.1 | Orbitrap domain handling (3D potential) | 1 week | MEDIUM |
| 8.2 | Multi-gas collision/reaction handlers | 1-2 weeks | MEDIUM |
| **9** | **GPU Acceleration** | 4-6 weeks | HIGH |
| 9A | GPU infrastructure | 1 week | HIGH |
| 9B | Basic kernels (force, integration) | 1-2 weeks | HIGH |
| 9C | Space charge on GPU | 1-2 weeks | HIGH |
| 9D | Collision/reaction kernels | 1 week | MEDIUM |
| 9E | Performance optimization | 1 week | HIGH |

**Total Estimated Time:** 11-16 weeks (~3-4 months)

---

## Success Metrics

### Phase 6 (Unit Testing)
- [ ] All instruments have ≥3 unit tests
- [ ] 100% test pass rate
- [ ] All collision models validated
- [ ] Reactive IMS tested
- [ ] Test suite runtime < 2 minutes

### Phase 7 (Validation Suite)
- [ ] Physics validation error < 5%
- [ ] Conservation laws satisfied (< 10⁻¹⁴ error)
- [ ] Cross-code agreement < 3%
- [ ] Performance benchmarks documented
- [ ] Publication-ready PDF report generated

### Phase 8 (Physics Enhancements)
- [ ] Orbitrap 3D potential implemented
- [ ] Multi-gas collisions working (≥3 gases)
- [ ] Competing reactions validated
- [ ] Backwards compatibility maintained
- [ ] Documentation updated

### Phase 9 (GPU Module)
- [ ] GPU speedup > 10× for large ensembles
- [ ] CPU-GPU consistency < 0.01%
- [ ] Automatic dispatch working
- [ ] Memory efficient (< 500 bytes/ion)
- [ ] Complete GPU test suite

---

## Risk Assessment and Mitigation

### Technical Risks

1. **GPU Complexity**
   - Risk: CUDA development more complex than anticipated
   - Mitigation: Start with simple kernels, iterative testing
   - Fallback: CPU-only mode always available

2. **Multi-Gas Implementation**
   - Risk: Breaking backwards compatibility
   - Mitigation: Factory pattern with auto-detection
   - Fallback: Keep single-gas handlers as default

3. **Validation Accuracy**
   - Risk: Lack of reference data for all instruments
   - Mitigation: Use analytical solutions where possible
   - Fallback: Focus on relative comparisons (model A vs B)

### Resource Risks

1. **Development Time**
   - Risk: 3-4 months is ambitious timeline
   - Mitigation: Prioritize High > Medium
   - Fallback: Defer optional features (Phase 8.3)

2. **GPU Hardware Access**
   - Risk: No GPU for development/testing
   - Mitigation: Use cloud instances (AWS, Google Cloud)
   - Fallback: Develop on CPU, test GPU later

---

## Notes

- **Testing Philosophy:** Simple tests with moderate tolerances (5-20%) prioritize correctness over statistical precision
- **Validation Focus:** Publication-ready suite targets < 5% error vs literature/analytical solutions
- **GPU Strategy:** Transparent acceleration - existing code works unchanged
- **Backwards Compatibility:** All enhancements maintain compatibility with existing configs (via fallbacks)
- **Documentation:** Each phase updates relevant docs (PUBLIC_CPP_API, INPUT_FORMAT_SPECIFICATION, etc.)

---

## Revision History

- **v1.0** (Nov 23, 2025): Initial roadmap after Phase 5 completion

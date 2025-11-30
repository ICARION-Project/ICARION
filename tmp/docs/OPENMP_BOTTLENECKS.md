# OpenMP Performance Bottlenecks Analysis

## Executive Summary

OpenMP parallelization in ICARION currently provides **no speedup** and often makes simulations **slower**. This document analyzes the root causes and proposes solutions.

## Test Results

### Small Problem Size (100 ions, 500 µs)

| Cores | Total Time | Speedup | Overhead |
|-------|------------|---------|----------|
| 1     | 2.57s      | 1.00x   | baseline |
| 2     | 2.82s      | 0.91x   | +0.25s   |
| 4     | 2.60s      | 0.99x   | +0.03s   |
| 8     | 2.82s      | 0.91x   | +0.25s   |
| 16    | 2.90s      | 0.89x   | +0.33s   |

**Result:** No speedup. Synchronization overhead (500k timesteps × 0.5µs = 250ms) dominates.

### Medium Problem Size (1000 ions, 500 µs, no space charge)

| Cores | Total Time | Speedup | Overhead |
|-------|------------|---------|----------|
| 1     | 50.6s      | 1.00x   | baseline |
| 2     | 89.7s      | 0.56x   | +39s     |
| 4     | 80.7s      | 0.63x   | +30s     |
| 8     | 87.5s      | 0.58x   | +37s     |
| 16    | 108.6s     | 0.47x   | +58s     |

**Result:** Massive slowdown! False sharing and cache coherency issues dominate.

### Large Problem Size (10000 ions, 10 µs, no space charge)

| Cores | Total Time | Speedup | Notes |
|-------|------------|---------|-------|
| 1     | 56.4s      | 1.00x   | baseline |
| 4     | >90s       | <0.63x  | timeout |
| 8     | >120s      | <0.47x  | timeout |

**Result:** Catastrophic slowdown! Memory bandwidth saturation and cache thrashing.

## Detailed Profiling Analysis

### Component Breakdown (100 ions, 1 core vs 8 cores)

| Component | 1 Core (ms) | 8 Cores (ms) | Overhead |
|-----------|-------------|--------------|----------|
| Integration | 1801.3 | 1868.2 | +67ms (+3.7%) |
| Collision | 134.9 | 138.1 | +3ms (+2.2%) |
| Domain Finding | 24.8 | 25.9 | +1ms (+4.4%) |
| Boundary Checks | 22.8 | 23.5 | +1ms (+3.1%) |
| **Missing/Sync** | 393 | 588 | **+195ms (+50%)** |

**Key Finding:** OpenMP synchronization overhead grows from 393ms to 588ms (+195ms).

### Per-Operation Cost

**1000 ions × 500k timesteps = 500M operations:**
- 1 core: 14s collision → 28ns/call
- 8 cores: (timeout after equivalent work)

**10000 ions × 50k timesteps = 500M operations:**
- 1 core: 57s collision → 114ns/call

**Cache Effect:** 4x slowdown (28ns → 114ns) when ion array doesn't fit in L2/L3 cache!

## Root Causes

### 1. **False Sharing** (Primary Issue)
```cpp
#pragma omp parallel for
for (int i = 0; i < n_ions; ++i) {
    IonState& ion = ions[i];  // ← Adjacent threads write to adjacent memory
    ion.pos = ...;  // ← Cache line ping-pong between cores
    ion.vel = ...;
}
```

**Problem:** 
- `IonState` is ~200 bytes (cache line = 64 bytes)
- Each ion spans 3-4 cache lines
- Adjacent threads (e.g., ions 0-1249 on core 0, ions 1250-2499 on core 1) write to cache lines that overlap
- CPU must invalidate and transfer cache lines between cores on every write
- With 500k timesteps, this causes **millions of cache misses**

### 2. **Memory Bandwidth Saturation**
- 10000 ions × 200 bytes = 2 MB
- L2 cache per core: 256 KB → **cache miss on every access**
- L3 cache shared: 64 MB → sufficient, but **bandwidth-limited**
- Multiple cores reading/writing simultaneously → memory bus saturation
- AMD Ryzen 9 7950X: ~85 GB/s max bandwidth
- Estimated: 10000 ions × 200 bytes × 10000 timesteps × 2 (read+write) = 40 GB
- At 8 cores: 40 GB / 56s = 0.7 GB/s → **FAR below hardware limit!**
- **Conclusion:** Not bandwidth, but **cache coherency protocol overhead**

### 3. **Synchronization Overhead**
```cpp
// Current implementation: Implicit barrier at end of every timestep
for (step = 0; step < 500000; ++step) {
    #pragma omp parallel for
    for (ion in ions) { ... }
    // ← BARRIER: All threads wait here (500k times!)
}
```

**Cost:** ~0.5µs per barrier × 500k steps = **250ms pure synchronization overhead**

### 4. **Work Distribution Imbalance**
- Static scheduling: `schedule(static, 256)`
- With 100 ions, 8 cores: 12.5 ions per core
- With 1000 ions, 8 cores: 125 ions per core
- With 10000 ions, 8 cores: 1250 ions per core
- **Problem:** Ions don't all remain active (some exit simulation)
- Some threads finish early, wait at barrier

### 5. **Cache Line Bouncing (MESI Protocol)**
When ion data is modified:
1. Core 0 writes ion[0] → owns cache line (Modified state)
2. Core 1 tries to read ion[1] (same cache line) → cache miss
3. Core 0 must flush cache line to L3
4. Core 1 loads cache line (Shared state)
5. Core 1 writes ion[1] → invalidates Core 0's cache
6. Repeat for every ion, every timestep → **billions of cache invalidations**

## Why SoA (Structure of Arrays) Solves This

### Current AoS (Array of Structures)
```cpp
struct IonState {
    Vec3 pos;      // 24 bytes
    Vec3 vel;      // 24 bytes
    double mass;   // 8 bytes
    double charge; // 8 bytes
    // ... ~200 bytes total
};
std::vector<IonState> ions(10000);  // ← Bad for parallelism!

// Access pattern:
ions[0].pos.x = ...;  // Cache line A
ions[1].pos.x = ...;  // Cache line A (FALSE SHARING!)
```

### Proposed SoA (Structure of Arrays)
```cpp
struct IonStateSoA {
    std::vector<double> pos_x;     // 10000 doubles contiguous
    std::vector<double> pos_y;     // 10000 doubles contiguous
    std::vector<double> pos_z;     // 10000 doubles contiguous
    std::vector<double> vel_x;     // etc.
    std::vector<double> vel_y;
    std::vector<double> vel_z;
    std::vector<double> mass;
    std::vector<double> charge;
    // ...
};

// Access pattern:
ions.pos_x[0] = ...;  // Cache line A
ions.pos_x[1] = ...;  // Cache line A (SAME DATA, OK!)
ions.pos_x[8] = ...;  // Cache line B (CLEAN SEPARATION!)
```

**Benefits:**
1. **No False Sharing:** Each thread works on different cache lines
   - Thread 0: `pos_x[0-1249]` → cache lines 0-156
   - Thread 1: `pos_x[1250-2499]` → cache lines 157-312
   - **Zero overlap!**

2. **Better Vectorization:** SIMD can process 4-8 doubles per instruction
   - AVX2: 4 doubles/cycle
   - AVX-512: 8 doubles/cycle
   - Compiler can auto-vectorize SoA loops easily

3. **Better Cache Utilization:** 
   - Reading `pos_x[0-7]` loads one cache line (64 bytes / 8 bytes = 8 doubles)
   - AoS: Reading 8 ions loads 8 cache lines (1600 bytes) but uses only 64 bytes of pos_x data
   - **25x more efficient cache usage**

4. **Predictable Memory Access:** Sequential reads → hardware prefetcher works perfectly

## Expected Speedup with SoA

**Conservative estimate:**
- Eliminate false sharing: 2-3x speedup
- Better vectorization: 1.5-2x speedup
- Better cache utilization: 1.5x speedup
- **Combined: 4-9x speedup at 8 cores**

**Realistic target:**
- 10000 ions, 10µs: 56s → 10-15s with 8 cores (4-6x speedup)
- 100 ions, 500µs: 2.57s → 2.57s (no change, too small)
- 1000 ions, 500µs: 50s → 12-15s with 8 cores (3-4x speedup)

## Implementation Roadmap

### Phase 1: Hybrid Approach (Backward Compatible)
1. Create `IonStateSoA` class alongside `IonState`
2. Add conversion methods: `toSoA()` / `fromSoA()`
3. Keep existing AoS interface for I/O and user code
4. Convert to SoA only in hot loops:
   ```cpp
   IonStateSoA ions_soa = convert_to_soa(ions);
   #pragma omp parallel for
   for (int i = 0; i < n_ions; ++i) {
       ions_soa.pos_x[i] += ions_soa.vel_x[i] * dt;
       // ...
   }
   convert_from_soa(ions_soa, ions);
   ```

### Phase 2: Full SoA Migration
1. Replace `std::vector<IonState>` with `IonStateSoA` throughout codebase
2. Update `HDF5Writer` to write directly from SoA
3. Update force computations to use SoA
4. Benchmark and optimize

### Phase 3: GPU Preparation
- SoA is also required for efficient GPU parallelization (CUDA/HIP)
- GPU kernels require coalesced memory access → SoA is mandatory

## Alternative Optimizations (Without SoA)

If SoA refactoring is too invasive, consider:

1. **Reduce Synchronization Frequency:**
   ```cpp
   #define SYNC_INTERVAL 100
   #pragma omp parallel
   for (step = 0; step < n_steps; ++step) {
       #pragma omp for nowait  // ← Skip barrier most of the time
       for (ion in ions) { ... }
       
       if (step % SYNC_INTERVAL == 0) {
           #pragma omp barrier  // ← Only sync every 100 steps
       }
   }
   ```
   **Expected:** Save ~200ms (500k barriers → 5k barriers)

2. **Cache Line Padding:**
   ```cpp
   struct alignas(64) IonState {  // ← Force 64-byte alignment
       Vec3 pos;
       Vec3 vel;
       char padding[64 - sizeof(pos) - sizeof(vel)];  // ← Waste space to avoid sharing
   };
   ```
   **Tradeoff:** 3x memory usage, but eliminates false sharing

3. **Larger Chunk Size:**
   ```cpp
   #pragma omp parallel for schedule(static, 1024)  // ← Larger chunks
   ```
   **Expected:** Reduce scheduling overhead, but doesn't fix false sharing

4. **Disable OpenMP for N < 10000:**
   ```cpp
   #pragma omp parallel for if(n_ions >= 10000)
   ```
   **Best short-term fix:** Avoid slowdown on small problems

## Conclusion

**Current Status:** OpenMP provides **negative speedup** due to false sharing and cache coherency overhead.

**Root Cause:** AoS (Array of Structures) memory layout is fundamentally incompatible with efficient parallelization.

**Solution:** Migrate to SoA (Structure of Arrays) to achieve 4-9x speedup on multi-core systems.

**Timeline:**
- Short-term: Disable OpenMP for N < 10000 ions
- Medium-term: Implement hybrid AoS/SoA with conversion
- Long-term: Full SoA migration + GPU support

**Critical:** SoA is **mandatory** for future GPU acceleration anyway, so this investment pays double dividends.

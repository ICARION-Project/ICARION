# Poisson Solver Validation Bench

**Purpose:** High-accuracy validation for publication-quality results  
**Runtime:** 30-60 minutes (not for CI, only for major validation)  
**Target:** <10% error vs analytical solutions

---

## 🎯 Why This Bench Exists

Regular unit tests are optimized for speed (128³ grids, ~15s runtime).  
This causes 28-42% errors due to **Dirichlet boundary condition artifacts**.

For **publication-ready validation**, we need:
- Much larger domains (20mm+) to minimize BC effects
- Higher resolution (256³-512³) for accurate gradients
- Multiple solver comparison (all 5 methods)
- Rigorous convergence analysis

---

## 📋 Test Suite

### 1. High-Resolution Point Charge (10 min)
**Grid:** 256³, 50μm cells, 12.8mm domain  
**Test Points:** r = 0.5, 1, 2, 3, 4 mm  
**Expected Accuracy:** <5% error  
**Validates:** Basic Poisson solver correctness

```cpp
TEST_CASE("ValidationBench: Point charge 256³", "[validation][expensive]") {
    const int N = 256;
    const double cell_size = 5e-5;
    const double domain_size = N * cell_size;  // 12.8mm
    
    // Single electron at center
    // Test at multiple radii
    // Check φ(r) = Q/(4πε₀r)
}
```

---

### 2. Uniform Charged Sphere (15 min)
**Grid:** 256³, 80μm cells, 20.48mm domain  
**Sphere:** 1mm radius, uniform charge distribution  
**Test Points:**
- Inside: r = 0.3, 0.5, 0.7, 0.9 mm
- Outside: r = 1.5, 2, 3, 4, 5 mm

**Expected Accuracy:** <10% error  
**Validates:** Distributed charge handling, boundary conditions

**Analytical Solutions:**
- Inside: φ(r) = Q/(8πε₀R) × (3 - r²/R²)
- Outside: φ(r) = Q/(4πε₀r)

```cpp
TEST_CASE("ValidationBench: Uniform sphere 256³", "[validation][expensive]") {
    const int N = 256;
    const double cell_size = 8e-5;  // 80μm
    const double domain_size = N * cell_size;  // 20.48mm
    
    const double sphere_radius = 0.001;  // 1mm
    const double total_charge = 1.6e-17;  // 100 elementary charges
    
    // Fill sphere with uniform charge density
    // Test inside and outside points
    // Validate against analytical solutions
}
```

---

### 3. Grid Convergence Study (20 min)
**Grids:** 64³, 128³, 256³ (same physical domain)  
**Expected:** Error reduces by ~4x with 2x refinement (O(h²))  

```cpp
TEST_CASE("ValidationBench: Convergence analysis", "[validation][expensive]") {
    std::vector<int> grid_sizes = {64, 128, 256};
    
    for (int N : grid_sizes) {
        // Keep domain size constant (16mm)
        // Measure error vs analytical
        // Verify O(h²) scaling
    }
    
    // Plot: log(error) vs log(h)
    // Slope should be ~2
}
```

---

### 4. Solver Method Comparison (15 min)
**Goal:** Verify all 5 solvers give same result  
**Solvers:** Gauss-Seidel, Red-Black SOR, CG, Multigrid, FFT  
**Grid:** 128³ (large enough to show differences)  

```cpp
TEST_CASE("ValidationBench: Solver comparison 128³", "[validation][expensive]") {
    // Same problem (point charge)
    // Solve with all 5 methods
    // Compare solutions (should agree to <0.1%)
    // Measure iteration counts
    // Measure wall time
}
```

**Expected Results:**
| Solver | Iterations | Time | Accuracy |
|--------|-----------|------|----------|
| Gauss-Seidel | ~3000 | 8s | Reference |
| Red-Black SOR | ~1500 | 6s | <0.1% |
| Conjugate Gradient | ~200 | 3s | <0.1% |
| Multigrid | ~50 | 2s | <0.1% |
| FFT | 1 | 1s | <0.1% |

---

## 🚀 Running the Validation Bench

```bash
# Build with optimization
cd build
cmake -DCMAKE_BUILD_TYPE=Release ..
make test_poisson_validation_bench

# Run full validation (30-60 minutes)
./tests/physics/spacecharge/test_poisson_validation_bench

# Run specific test
./test_poisson_validation_bench "[validation][sphere]"

# Generate plots
python3 ../validation/analyze_poisson_results.py
```

---

## 📊 Success Criteria

For **v1.0.0 publication readiness:**
- ✅ Point charge: <5% error at all test radii
- ✅ Sphere inside: <10% error
- ✅ Sphere outside: <10% error
- ✅ Convergence: O(h²) verified (slope 1.9-2.1)
- ✅ Solver agreement: <0.1% difference

For **future improvements (v1.1+):**
- Implement Neumann boundary conditions
- Implement analytical/Robin boundary conditions
- Adaptive mesh refinement near charges
- Higher-order discretization (4th order)

---

## 🔧 Implementation Priority

**Phase 1 (v1.0.0):** Basic validation bench
- [ ] Create `test_poisson_validation_bench.cpp`
- [ ] Implement 256³ point charge test
- [ ] Implement 256³ sphere test
- [ ] Run and document results

**Phase 2 (v1.1):** Advanced validation
- [ ] Convergence study with plots
- [ ] Solver comparison with performance metrics
- [ ] Neumann BC implementation
- [ ] Publish validation paper

---

## 📝 Notes

- **Memory:** 256³ double grid = 256³ × 8 bytes = 128 MB (manageable)
- **Parallelization:** All solvers support OpenMP
- **Reproducibility:** Use fixed RNG seed for charge placement
- **Documentation:** Save results to validation/poisson_results.json

---

**Created:** 2025-11-24  
**Status:** TODO - implement for v1.0.0 validation  
**Owner:** Space Charge Team

# ICARION Missing Features Summary

**Date:** November 29, 2025  
**Branch:** `feature/gpu-acceleration`  
**Status:** Field arrays complete ✅, 4 major features remaining

---

## 🎯 Overview

After completing Phase 7 (Field Arrays), ICARION still needs 4 major features before GPU acceleration can be fully leveraged:

1. **RK45 Adaptive Integrator** - Better accuracy for long simulations
2. **Boris Pusher** - Correct handling of magnetic fields
3. **Domain Boundaries** - Reflection, absorption, re-injection
4. **GPU Batch Processing** - Collision and integration acceleration

---

## 📊 Feature Priority Matrix

| Feature | Effort | Impact | GPU Benefit | Priority | Status |
|---------|--------|--------|-------------|----------|--------|
| **Domain Boundaries** | 3-4 days | 🔴 CRITICAL | Medium | 🥇 HIGHEST | Partial |
| **RK45 Adaptive** | 2-3 days | 🟡 HIGH | High | 🥈 HIGH | Not started |
| **Boris Pusher** | 2-3 days | 🟡 HIGH | High | 🥉 MEDIUM | Not started |
| **GPU Collisions** | 4-5 days | 🟢 MEDIUM | Very High | 🏅 HIGH | CPU complete |

---

## 1. Domain Boundaries (CRITICAL)

### Current State
- ✅ Boundary **detection** implemented
- ❌ No boundary **actions** (reflection, absorption, re-injection)
- ❌ Ions just pass through geometry walls

### What's Missing
```cpp
// Current: Detection only
bool is_outside = !domain.geometry.contains(ion.pos);
if (is_outside) {
    ion.active = false;  // Just deactivate
}

// Needed: Actions
if (boundary_checker->check_boundary(ion, domain)) {
    BoundaryAction* action = boundary_config.get_action(boundary_type);
    action->apply(ion, normal_vector);
    // action could be: Reflect, Absorb, ReInject, etc.
}
```

### Implementation Plan

#### Step 1: Define BoundaryAction Interface
**File:** `src/core/integrator/boundaries/BoundaryAction.h`
```cpp
class BoundaryAction {
public:
    virtual ~BoundaryAction() = default;
    virtual void apply(IonState& ion, const Vec3& normal) = 0;
    virtual std::string name() const = 0;
};
```

#### Step 2: Implement Reflection
**File:** `src/core/integrator/boundaries/ReflectionAction.h`
```cpp
class ReflectionAction : public BoundaryAction {
public:
    enum Type { SPECULAR, DIFFUSE, THERMAL };
    
    ReflectionAction(Type type, double accommodation_coeff = 1.0);
    void apply(IonState& ion, const Vec3& normal) override;
    
private:
    Type type_;
    double accommodation_coeff_;
};
```

**Physics:**
- **Specular:** `v_reflected = v - 2(v·n)n` (mirror reflection)
- **Diffuse:** Randomize velocity direction (cosine distribution)
- **Thermal:** Re-sample velocity from Maxwell-Boltzmann at wall temperature

#### Step 3: Implement Absorption
**File:** `src/core/integrator/boundaries/AbsorptionAction.cpp`
```cpp
class AbsorptionAction : public BoundaryAction {
public:
    void apply(IonState& ion, const Vec3& normal) override {
        ion.active = false;  // Remove ion from simulation
    }
    std::string name() const override { return "Absorption"; }
};
```

#### Step 4: Implement Re-Injection
**File:** `src/core/integrator/boundaries/ReInjectionAction.cpp`
```cpp
class ReInjectionAction : public BoundaryAction {
public:
    ReInjectionAction(Vec3 injection_point, double temperature_K);
    void apply(IonState& ion, const Vec3& normal) override {
        // Move to injection point
        ion.pos = injection_point_;
        // Re-sample thermal velocity
        ion.vel = sample_maxwell_boltzmann(temperature_K_);
        ion.active = true;
    }
    
private:
    Vec3 injection_point_;
    double temperature_K_;
};
```

#### Step 5: Integrate into SimulationEngine
**File:** `src/core/integrator/SimulationEngine.cpp`
```cpp
void SimulationEngine::check_domain_boundaries(IonState& ion, int current_domain) {
    const auto& domain = config_.domains[current_domain];
    
    // Check if ion left domain
    if (!domain.geometry.contains(ion.pos)) {
        Vec3 normal = domain.geometry.get_normal(ion.pos);
        
        // Get boundary action from config
        auto* action = domain.boundary_action.get();
        action->apply(ion, normal);
        
        // Log event
        logger_->debug("Boundary event: {} at pos ({},{},{})", 
                      action->name(), ion.pos.x, ion.pos.y, ion.pos.z);
    }
}
```

### Configuration Format
```json
{
  "domains": [{
    "name": "drift_region",
    "geometry": { ... },
    "boundary": {
      "type": "reflection",
      "reflection_type": "thermal",
      "temperature_K": 300.0,
      "accommodation_coeff": 0.8
    }
  }]
}
```

### Testing Requirements
1. **Specular Reflection:** 45° incident → 45° reflected
2. **Diffuse Reflection:** Uniform angular distribution (cos θ weighted)
3. **Energy Conservation:** `|v_before| ≈ |v_after|` for elastic
4. **Thermal Accommodation:** Velocity distribution → Maxwell-Boltzmann after many bounces

### Estimated Effort
- **Boundary actions:** 1 day
- **Integration:** 1 day
- **Testing & validation:** 1-2 days
- **Total:** 3-4 days

---

## 2. RK45 Adaptive Integrator (HIGH PRIORITY)

### Why Needed
- RK4 uses fixed timesteps → inefficient for varying dynamics
- Adaptive RK45 adjusts `dt` based on error estimation
- Better accuracy with fewer steps (10-100× speedup for stiff systems)

### Current State
- ✅ RK4 implemented and working
- ✅ Integration strategy interface defined
- ❌ No adaptive timestep control

### Algorithm: Runge-Kutta-Fehlberg (RK45)
```
k1 = f(t, y)
k2 = f(t + c2*dt, y + dt*(a21*k1))
k3 = f(t + c3*dt, y + dt*(a31*k1 + a32*k2))
k4 = f(t + c4*dt, y + dt*(a41*k1 + a42*k2 + a43*k3))
k5 = f(t + c5*dt, y + dt*(a51*k1 + a52*k2 + a53*k3 + a54*k4))
k6 = f(t + c6*dt, y + dt*(a61*k1 + a62*k2 + a63*k3 + a64*k4 + a65*k5))

y_4th = y + dt*(b1*k1 + b2*k2 + b3*k3 + b4*k4 + b5*k5)         # 4th order
y_5th = y + dt*(b1_*k1 + b2_*k2 + b3_*k3 + b4_*k4 + b5_*k5 + b6_*k6)  # 5th order

error = |y_5th - y_4th|
```

### Implementation Plan

#### Step 1: RK45Strategy Class
**File:** `src/core/integrator/strategies/RK45Strategy.h`
```cpp
class RK45Strategy : public IIntegrationStrategy {
public:
    RK45Strategy(double tolerance = 1e-6, 
                 double dt_min = 1e-12, 
                 double dt_max = 1e-6);
    
    void step(IonState& ion, double t, double dt, 
              const ForceRegistry& forces,
              const std::vector<IonState>& all_ions) override;
    
    double get_adaptive_dt() const { return dt_adaptive_; }
    size_t get_rejected_steps() const { return rejected_steps_; }
    
private:
    double tolerance_;
    double dt_min_, dt_max_;
    double dt_adaptive_;
    size_t rejected_steps_;
    
    std::array<Vec3, 6> k_pos_;  // k1-k6 for position
    std::array<Vec3, 6> k_vel_;  // k1-k6 for velocity
};
```

#### Step 2: Butcher Tableau (RK45 Coefficients)
```cpp
// Dormand-Prince coefficients (most common RK45 variant)
static constexpr double c2 = 1.0/5.0;
static constexpr double c3 = 3.0/10.0;
static constexpr double c4 = 4.0/5.0;
static constexpr double c5 = 8.0/9.0;
static constexpr double c6 = 1.0;

static constexpr double a21 = 1.0/5.0;
static constexpr double a31 = 3.0/40.0, a32 = 9.0/40.0;
// ... (full tableau ~20 constants)

static constexpr double b1 = 35.0/384.0;  // 5th order
static constexpr double b1_star = 5179.0/57600.0;  // 4th order
// ... (error estimation coefficients)
```

#### Step 3: Error Estimation & Step Control
```cpp
void RK45Strategy::step(IonState& ion, double t, double dt, 
                       const ForceRegistry& forces, 
                       const std::vector<IonState>& all_ions) {
    double dt_try = dt_adaptive_;
    bool accepted = false;
    
    while (!accepted) {
        // Compute k1-k6
        compute_stages(ion, t, dt_try, forces, all_ions);
        
        // 4th and 5th order solutions
        Vec3 pos_4 = compute_solution_4th();
        Vec3 pos_5 = compute_solution_5th();
        Vec3 vel_4 = compute_velocity_4th();
        Vec3 vel_5 = compute_velocity_5th();
        
        // Error estimate
        double error_pos = (pos_5 - pos_4).magnitude();
        double error_vel = (vel_5 - vel_4).magnitude();
        double error = std::max(error_pos, error_vel);
        
        // Accept or reject
        if (error < tolerance_) {
            ion.pos = pos_5;
            ion.vel = vel_5;
            accepted = true;
        } else {
            rejected_steps_++;
        }
        
        // Adjust timestep
        double safety = 0.9;
        double dt_new = safety * dt_try * std::pow(tolerance_ / error, 0.2);
        dt_adaptive_ = std::clamp(dt_new, dt_min_, dt_max_);
        dt_try = dt_adaptive_;
    }
}
```

### Testing Requirements
1. **Harmonic Oscillator:** `x(t) = A*cos(ωt)` → error < tolerance
2. **Kepler Orbit:** Energy conservation over 1000 orbits
3. **Stiff System:** Van der Pol oscillator (adaptive wins vs RK4)
4. **Step Rejection:** Verify rejected steps for large errors

### Configuration
```json
{
  "physics": {
    "integrator": "RK45",
    "rk45_tolerance": 1e-6,
    "rk45_dt_min": 1e-12,
    "rk45_dt_max": 1e-6
  }
}
```

### Estimated Effort
- **Implementation:** 1 day
- **Testing & validation:** 1-2 days
- **Total:** 2-3 days

---

## 3. Boris Pusher (MEDIUM PRIORITY)

### Why Needed
- Standard integrators (RK4, RK45) have energy drift for `E + B` fields
- Boris pusher is **symplectic** → conserves phase space volume
- Standard in PIC codes for magnetized plasmas

### Current State
- ✅ Electric field force working
- ✅ Magnetic field force exists but not tested
- ❌ No specialized B-field integrator

### Algorithm
```
1. v^- = v^n + (q/m) * E * (dt/2)          # Half electric impulse
2. t = (q/m) * B * (dt/2)                  # Rotation vector
3. s = 2t / (1 + |t|^2)                    # Scaled rotation
4. v^+ = v^- + (v^- + v^- × t) × s         # Rotation
5. v^{n+1} = v^+ + (q/m) * E * (dt/2)      # Half electric impulse
6. x^{n+1} = x^n + v^{n+1} * dt            # Position update
```

### Implementation Plan

#### Step 1: BorisPusherStrategy Class
**File:** `src/core/integrator/strategies/BorisPusherStrategy.h`
```cpp
class BorisPusherStrategy : public IIntegrationStrategy {
public:
    void step(IonState& ion, double t, double dt,
              const ForceRegistry& forces,
              const std::vector<IonState>& all_ions) override;
    
private:
    Vec3 rotate_velocity(const Vec3& v_minus, const Vec3& B, double q_over_m, double dt);
};
```

#### Step 2: Rotation Implementation
```cpp
Vec3 BorisPusherStrategy::rotate_velocity(const Vec3& v_minus, const Vec3& B, 
                                         double q_over_m, double dt) {
    Vec3 t = B * (q_over_m * dt / 2.0);
    double t_mag_sq = t.dot(t);
    Vec3 s = t * (2.0 / (1.0 + t_mag_sq));
    
    Vec3 v_prime = v_minus + v_minus.cross(t);
    Vec3 v_plus = v_minus + v_prime.cross(s);
    
    return v_plus;
}
```

### Testing Requirements
1. **Uniform B-field:** Circular motion, period = `2πm/(qB)`
2. **Energy Conservation:** `|v|` constant over 1000 gyrations
3. **E + B Drift:** `v_drift = E × B / B²`
4. **Compare vs RK4:** Boris should have better energy conservation

### Configuration
```json
{
  "physics": {
    "integrator": "Boris",
    "magnetic_field": {
      "type": "uniform",
      "B_T": [0.0, 0.0, 1.0]
    }
  }
}
```

### Estimated Effort
- **Implementation:** 1 day
- **Testing & validation:** 1-2 days
- **Total:** 2-3 days

---

## 4. GPU Collision Handler (HIGH PERFORMANCE IMPACT)

### Why High Priority
- **15-25% of CPU time** in typical simulations
- **5-20× speedup potential** with GPU
- Already validated on CPU (EHSS, HSS, OU models)
- Good candidate for first GPU acceleration

### Current State
- ✅ EHSS collision model implemented and validated
- ✅ HSS collision model implemented and validated
- ✅ Thermalization tests passing
- ❌ No GPU implementation
- ❌ No batch processing

### GPU Challenges
1. **Random Number Generation:** Need GPU RNG (cuRAND)
2. **Cross-Section Lookup:** Large tables → texture memory
3. **Divergent Branching:** Collision vs no-collision
4. **Memory Bandwidth:** Ion state read/write

### Implementation Plan

#### Step 1: GPU RNG Integration
**File:** `src/core/gpu/GPURandom.h`
```cpp
class GPURandom {
public:
    GPURandom(size_t n_ions, unsigned long seed);
    ~GPURandom();
    
    // Generate random numbers for all ions
    void generate_uniform(float* d_output, size_t count);
    void generate_gaussian(float* d_output, size_t count);
    
private:
    curandGenerator_t generator_;
    size_t n_ions_;
};
```

#### Step 2: Collision Data Structure (GPU-friendly)
```cpp
struct CollisionDataGPU {
    float* d_cross_sections;      // [n_species x n_energy_bins]
    float* d_energy_bins;          // [n_energy_bins]
    int n_energy_bins;
    float gas_density_m3;          // Number density [m^-3]
    float gas_temperature_K;
    float gas_mass_kg;
};
```

#### Step 3: EHSS Collision Kernel
**File:** `src/core/gpu/kernels/collision_kernels.cu`
```cuda
__global__ void process_collisions_ehss_kernel(
    IonStateGPU* d_ions,
    int n_ions,
    const CollisionDataGPU* d_collision_data,
    float* d_random_uniform,
    float* d_random_gaussian,
    double dt
) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_ions) return;
    
    IonStateGPU& ion = d_ions[idx];
    if (!ion.active) return;
    
    // Compute collision probability
    float v_rel = sqrtf(ion.vel.x*ion.vel.x + 
                       ion.vel.y*ion.vel.y + 
                       ion.vel.z*ion.vel.z);
    float cross_section = lookup_cross_section(d_collision_data, v_rel);
    float P_coll = 1.0f - expf(-d_collision_data->gas_density_m3 * 
                               cross_section * v_rel * dt);
    
    // Collision decision
    if (d_random_uniform[idx] < P_coll) {
        // Apply EHSS collision
        apply_ehss_collision(ion, d_collision_data, 
                           &d_random_gaussian[idx*6]);
    }
}
```

#### Step 4: GPUCollisionHandler Wrapper
**File:** `src/core/physics/collisions/GPUCollisionHandler.h`
```cpp
class GPUCollisionHandler : public ICollisionHandler {
public:
    GPUCollisionHandler(
        std::shared_ptr<ICollisionHandler> cpu_handler,
        size_t gpu_threshold = 10000
    );
    
    void process_collisions(
        std::vector<IonState>& ions,
        double dt,
        const DomainConfig& domain
    ) override;
    
private:
    std::shared_ptr<ICollisionHandler> cpu_handler_;
    std::unique_ptr<GPURandom> rng_;
    std::unique_ptr<GPUMemoryPool> memory_pool_;
    size_t gpu_threshold_;
    
    void process_collisions_gpu(
        std::vector<IonState>& ions,
        double dt,
        const DomainConfig& domain
    );
};
```

### Testing Requirements
1. **Bit-Exact Validation:** GPU matches CPU for same RNG seed
2. **Thermalization:** GPU reaches Maxwell-Boltzmann (same as CPU tests)
3. **Performance:** GPU > 10× faster for N > 10k ions
4. **Scaling:** Test 1k → 1M ions

### Estimated Effort
- **GPU RNG setup:** 1 day
- **Collision kernels:** 2 days
- **Wrapper & integration:** 1 day
- **Testing & validation:** 1-2 days
- **Total:** 4-5 days

---

## 🎯 Recommended Implementation Order

### Phase A: Critical CPU Features (Week 1-2)
1. **Domain Boundaries** (3-4 days) - CRITICAL for realistic physics
2. **RK45 Adaptive** (2-3 days) - Big accuracy win
3. **Boris Pusher** (2-3 days) - Enables magnetic simulations

### Phase B: GPU Acceleration (Week 3-5)
4. **GPU Collisions** (4-5 days) - Biggest performance win
5. **GPU Integration** (3-4 days) - RK4/RK45 batch processing
6. **GPU Boundaries** (2-3 days) - Complete GPU path

### Total Effort: ~22-31 days (4-6 weeks)

---

## 🔧 Getting Started

### For Domain Boundaries (NEXT STEP)
```bash
cd /home/chsch95/ICARION
git checkout feature/gpu-acceleration  # Or create new branch

# Create files
mkdir -p src/core/integrator/boundaries
touch src/core/integrator/boundaries/BoundaryAction.h
touch src/core/integrator/boundaries/ReflectionAction.h
touch src/core/integrator/boundaries/ReflectionAction.cpp
touch src/core/integrator/boundaries/AbsorptionAction.cpp
touch tests/integrator/test_boundary_actions.cpp

# Update CMakeLists.txt
# Implement classes
# Write tests
# Integrate into SimulationEngine
```

### For RK45 Integrator (ALTERNATIVE START)
```bash
cd /home/chsch95/ICARION
git checkout feature/gpu-acceleration

# Create files
touch src/core/integrator/strategies/RK45Strategy.h
touch src/core/integrator/strategies/RK45Strategy.cpp
touch tests/integrator/test_rk45_adaptive.cpp

# Update CMakeLists.txt
# Implement RK45 algorithm with Butcher tableau
# Write validation tests (harmonic oscillator, Kepler)
# Add configuration support
```

---

## 📚 References

### RK45
- Dormand, J. R., & Prince, P. J. (1980). "A family of embedded Runge-Kutta formulae"
- Press et al., *Numerical Recipes* Chapter 17.2

### Boris Pusher
- Boris, J. P. (1970). "Relativistic plasma simulation-optimization of a hybrid code"
- Birdsall & Langdon, *Plasma Physics via Computer Simulation*

### Boundaries
- Bird, G. A. (1994). *Molecular Gas Dynamics and the Direct Simulation of Gas Flows*
- Cercignani, C. (1988). *The Boltzmann Equation and Its Applications*

---

**Last Updated:** November 29, 2025  
**Next Milestone:** Complete Domain Boundaries (3-4 days)

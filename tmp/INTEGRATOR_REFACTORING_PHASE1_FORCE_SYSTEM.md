# 🔧 PHASE 1: FORCE SYSTEM REFACTORING

**Branch:** `refactor/force-system`  
**Duration:** 1 week (5 working days)  
**Status:** 🔵 Ready to Start  
**Dependencies:** None (self-contained)

---

## 🎯 OBJECTIVES

### Primary Goals
1. ✅ **Modular force computation** with clean interfaces
2. ✅ **Eliminate instrument-specific if/else chains**
3. ✅ **Enable force composition** (add/remove forces dynamically)
4. ✅ **Full unit test coverage** (>90%)
5. ✅ **Zero performance regression** (<5% overhead)

### Secondary Goals
- Clean separation between physics and field sampling
- Preparation for GPU acceleration (CUDA kernels)
- Foundation for Phase 4 (integrators need derivatives)

---

## 📐 ARCHITECTURE DESIGN

### Core Components

#### 1. **IForce Interface** (Abstract Base Class)

```cpp
// src/core/physics/forces/IForce.h

namespace ICARION::physics {

/**
 * @brief Force contribution interface
 * 
 * All forces (electric, magnetic, damping, space charge) implement this.
 * Enables modular force composition via ForceRegistry.
 */
class IForce {
public:
    virtual ~IForce() = default;
    
    /**
     * @brief Compute force contribution for single ion
     * @param ion Current ion state
     * @param t Current time [s]
     * @param context Optional context (e.g., all ions for space charge)
     * @return Force vector [N]
     */
    virtual Vec3 compute(
        const IonState& ion,
        double t,
        const ForceContext& context
    ) const = 0;
    
    /**
     * @brief Check if this force applies to given ion
     * @return true if force should be computed
     */
    virtual bool applies_to(const IonState& ion) const { 
        return true; 
    }
    
    /**
     * @brief Get force name (for logging/debugging)
     */
    virtual std::string name() const = 0;
};

} // namespace ICARION::physics
```

**Design Rationale:**
- Simple interface (single responsibility)
- Const-correct (forces don't modify state)
- Virtual destructor (proper cleanup)
- Optional `applies_to()` for conditional forces

---

#### 2. **ForceContext** (Context Struct)

```cpp
// src/core/physics/forces/ForceContext.h

namespace ICARION::physics {

/**
 * @brief Context for force computation
 * 
 * Contains shared data needed by multiple forces (e.g., all ions for space charge).
 */
struct ForceContext {
    // === Field Provider (optional) ===
    const IFieldProvider* field_provider = nullptr;
    
    // === Domain Configuration ===
    const DomainConfig* domain = nullptr;
    
    // === Ion Ensemble (for space charge) ===
    const std::vector<IonState>* all_ions = nullptr;
    
    // === Environment (for damping forces) ===
    double temperature_K = 300.0;
    double pressure_Pa = 101325.0;
    Vec3 gas_velocity_ms{0, 0, 0};
};

} // namespace ICARION::physics
```

**Design Rationale:**
- Single struct avoids parameter explosion
- Pointers allow optional data (nullptr = not needed)
- No ownership (context is transient)

---

#### 3. **ElectricFieldForce** (Primary Force)

```cpp
// src/core/physics/forces/ElectricFieldForce.h

namespace ICARION::physics {

/**
 * @brief Electric field force: F = q * E
 * 
 * Supports multiple instruments (IMS, FTICR, LQIT, TOF, Orbitrap).
 * Uses field provider if available, otherwise analytical fields.
 */
class ElectricFieldForce : public IForce {
public:
    /**
     * @brief Construct electric field force
     * @param domain Domain configuration (instrument type, geometry, fields)
     * @param field_provider Optional field provider (nullptr = analytical)
     */
    ElectricFieldForce(
        const DomainConfig& domain,
        std::shared_ptr<IFieldProvider> field_provider = nullptr
    );
    
    Vec3 compute(
        const IonState& ion, 
        double t, 
        const ForceContext& context
    ) const override;
    
    std::string name() const override { return "ElectricField"; }
    
private:
    const DomainConfig& domain_;
    std::shared_ptr<IFieldProvider> field_provider_;
    
    // === Instrument-specific field calculations ===
    Vec3 compute_ims_field(const IonState& ion, double t) const;
    Vec3 compute_fticr_field(const IonState& ion, double t) const;
    Vec3 compute_lqit_field(const IonState& ion, double t) const;
    Vec3 compute_tof_field(const IonState& ion, double t) const;
    Vec3 compute_orbitrap_field(const IonState& ion, double t) const;
};

} // namespace ICARION::physics
```

**Implementation Sketch:**

```cpp
// src/core/physics/forces/ElectricFieldForce.cpp

Vec3 ElectricFieldForce::compute(
    const IonState& ion, 
    double t, 
    const ForceContext& context
) const {
    Vec3 E_field{0, 0, 0};
    
    // Priority 1: Use field provider (grid-based or analytical)
    if (field_provider_) {
        E_field = field_provider_->sample_field(ion.pos, t);
    } 
    // Priority 2: Use context field provider
    else if (context.field_provider) {
        E_field = context.field_provider->sample_field(ion.pos, t);
    }
    // Priority 3: Instrument-specific analytical fields
    else {
        switch (domain_.instrument) {
            case Instrument::IMS:
                E_field = compute_ims_field(ion, t);
                break;
            case Instrument::FTICR:
                E_field = compute_fticr_field(ion, t);
                break;
            case Instrument::LQIT:
                E_field = compute_lqit_field(ion, t);
                break;
            case Instrument::TOF:
                E_field = compute_tof_field(ion, t);
                break;
            case Instrument::Orbitrap:
                E_field = compute_orbitrap_field(ion, t);
                break;
            default:
                throw std::runtime_error("Unknown instrument type");
        }
    }
    
    // F = q * E
    return E_field * ion.ion_charge_C;
}

Vec3 ElectricFieldForce::compute_ims_field(const IonState& ion, double t) const {
    // Axial DC field
    double E_axial = domain_.fields.dc.axial_V / domain_.geometry.length_m;
    
    // Optional radial quadrupole field
    Vec3 E{0, 0, E_axial};
    
    if (std::abs(domain_.fields.dc.quad_V) > 1e-9) {
        E.x = domain_.fields.dc.quad_V * ion.pos.x / std::pow(domain_.geometry.radius_m, 2);
        E.y = -domain_.fields.dc.quad_V * ion.pos.y / std::pow(domain_.geometry.radius_m, 2);
    }
    
    return E;
}

// ... other instrument implementations
```

---

#### 4. **MagneticFieldForce** (Lorentz Force)

```cpp
// src/core/physics/forces/MagneticFieldForce.h

namespace ICARION::physics {

/**
 * @brief Magnetic field force: F = q * (v × B)
 * 
 * Used in FTICR, Orbitrap, and other magnetic confinement instruments.
 */
class MagneticFieldForce : public IForce {
public:
    MagneticFieldForce(const MagneticFieldConfig& config);
    
    Vec3 compute(const IonState& ion, double t, const ForceContext& context) const override;
    
    bool applies_to(const IonState& ion) const override {
        return config_.enabled;  // Only if B-field enabled
    }
    
    std::string name() const override { return "MagneticField"; }
    
private:
    const MagneticFieldConfig& config_;
};

} // namespace ICARION::physics
```

**Implementation:**

```cpp
Vec3 MagneticFieldForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    // Lorentz force: F = q * (v × B)
    Vec3 B{0, 0, config_.field_strength_T};  // Typically axial
    
    // Optional gradient (for FTICR shimming)
    if (std::abs(config_.field_gradient_T_m) > 1e-9) {
        B.z += config_.field_gradient_T_m * ion.pos.z;
    }
    
    return cross(ion.vel, B) * ion.ion_charge_C;
}
```

---

#### 5. **DampingForce** (Friction/Langevin)

```cpp
// src/core/physics/forces/DampingForce.h

namespace ICARION::physics {

/**
 * @brief Damping force for collision models
 * 
 * Supports:
 * - Friction model: F = -γ * v
 * - Langevin model: F = -γ * v + sqrt(2*k_B*T*γ) * ξ(t)
 */
class DampingForce : public IForce {
public:
    enum class Model { Friction, Langevin };
    
    DampingForce(
        Model model,
        double damping_coefficient_gamma,
        double temperature_K = 300.0
    );
    
    Vec3 compute(const IonState& ion, double t, const ForceContext& context) const override;
    
    std::string name() const override { 
        return model_ == Model::Friction ? "Friction" : "Langevin"; 
    }
    
private:
    Model model_;
    double gamma_;
    double temperature_K_;
    mutable EhssRng rng_;  // For Langevin stochastic force
};

} // namespace ICARION::physics
```

---

#### 6. **ForceRegistry** (Manager)

```cpp
// src/core/physics/forces/ForceRegistry.h

namespace ICARION::physics {

/**
 * @brief Manages all active forces for simulation
 * 
 * Allows modular force composition without changing integrator.
 */
class ForceRegistry {
public:
    /**
     * @brief Add force to registry
     * @param force Force implementation (ownership transferred)
     */
    void add_force(std::unique_ptr<IForce> force);
    
    /**
     * @brief Compute total force on ion
     * @param ion Current ion state
     * @param t Current time [s]
     * @param context Optional context (e.g., all ions)
     * @return Total force vector [N]
     */
    Vec3 compute_total_force(
        const IonState& ion,
        double t,
        const ForceContext& context = {}
    ) const;
    
    /**
     * @brief Get all registered forces
     */
    const std::vector<std::unique_ptr<IForce>>& forces() const { 
        return forces_; 
    }
    
    /**
     * @brief Clear all forces
     */
    void clear() { forces_.clear(); }
    
private:
    std::vector<std::unique_ptr<IForce>> forces_;
};

} // namespace ICARION::physics
```

**Implementation:**

```cpp
void ForceRegistry::add_force(std::unique_ptr<IForce> force) {
    forces_.push_back(std::move(force));
}

Vec3 ForceRegistry::compute_total_force(
    const IonState& ion,
    double t,
    const ForceContext& context
) const {
    Vec3 total_force{0, 0, 0};
    
    for (const auto& force : forces_) {
        if (force->applies_to(ion)) {
            total_force += force->compute(ion, t, context);
        }
    }
    
    return total_force;
}
```

---

## 📁 FILE STRUCTURE

```
src/core/physics/forces/
├── IForce.h                    # Force interface (abstract base class)
├── ForceContext.h              # Context struct for force computation
├── ElectricFieldForce.h        # Electric field force header
├── ElectricFieldForce.cpp      # Electric field force implementation
├── MagneticFieldForce.h        # Magnetic field force header
├── MagneticFieldForce.cpp      # Magnetic field force implementation
├── DampingForce.h              # Damping/friction force header
├── DampingForce.cpp            # Damping/friction force implementation
├── SpaceChargeForce.h          # Space charge force header
├── SpaceChargeForce.cpp        # Space charge force implementation
├── ForceRegistry.h             # Force manager header
└── ForceRegistry.cpp           # Force manager implementation

tests/physics/forces/
├── CMakeLists.txt              # Test build configuration
├── test_electric_field_force.cpp    # Unit tests for electric force
├── test_magnetic_field_force.cpp    # Unit tests for magnetic force
├── test_damping_force.cpp           # Unit tests for damping force
├── test_space_charge_force.cpp      # Unit tests for space charge
└── test_force_registry.cpp          # Unit tests for force registry
```

---

## ✅ IMPLEMENTATION CHECKLIST

### Day 1: Interface & Context
- [ ] Create `src/core/physics/forces/` directory
- [ ] Implement `IForce.h` interface
- [ ] Implement `ForceContext.h` struct
- [ ] Write basic documentation

### Day 2: Electric Field Force
- [ ] Implement `ElectricFieldForce.h`
- [ ] Implement `ElectricFieldForce.cpp`
- [ ] Add IMS field calculation
- [ ] Add FTICR field calculation
- [ ] Add LQIT field calculation
- [ ] Add TOF field calculation
- [ ] Add Orbitrap field calculation
- [ ] Write unit tests (`test_electric_field_force.cpp`)

### Day 3: Magnetic & Damping Forces
- [ ] Implement `MagneticFieldForce.h` + `.cpp`
- [ ] Implement `DampingForce.h` + `.cpp`
- [ ] Write unit tests for both

### Day 4: Space Charge & Registry
- [ ] Implement `SpaceChargeForce.h` + `.cpp`
- [ ] Implement `ForceRegistry.h` + `.cpp`
- [ ] Write unit tests

### Day 5: Integration & Benchmarking
- [ ] Integration tests (multiple forces together)
- [ ] Performance benchmarks
- [ ] Code review
- [ ] Merge to `integrator-refactor`

---

## 🧪 TESTING STRATEGY

### Unit Tests (Per Force)

```cpp
// tests/physics/forces/test_electric_field_force.cpp

TEST_CASE("ElectricFieldForce - IMS uniform field", "[forces][electric]") {
    // Setup
    DomainConfig domain;
    domain.instrument = Instrument::IMS;
    domain.geometry.length_m = 0.1;
    domain.fields.dc.axial_V = 250.0;
    
    ElectricFieldForce force(domain);
    
    // Ion at origin
    IonState ion;
    ion.pos = {0, 0, 0.05};  // Midpoint
    ion.ion_charge_C = 1.602e-19;  // +1 charge
    
    // Compute force
    ForceContext ctx;
    Vec3 F = force.compute(ion, 0.0, ctx);
    
    // Verify
    double expected_E = 250.0 / 0.1;  // V/m
    double expected_F = ion.ion_charge_C * expected_E;
    
    REQUIRE(F.x == Approx(0.0));
    REQUIRE(F.y == Approx(0.0));
    REQUIRE(F.z == Approx(expected_F));
}
```

### Integration Tests

```cpp
TEST_CASE("ForceRegistry - multiple forces", "[forces][registry]") {
    // Setup domain
    DomainConfig domain;
    // ... configure domain
    
    // Create registry
    ForceRegistry registry;
    
    // Add electric force
    registry.add_force(std::make_unique<ElectricFieldForce>(domain));
    
    // Add magnetic force
    MagneticFieldConfig mag_config;
    mag_config.enabled = true;
    mag_config.field_strength_T = 3.0;
    registry.add_force(std::make_unique<MagneticFieldForce>(mag_config));
    
    // Ion with velocity
    IonState ion;
    ion.pos = {0, 0, 0};
    ion.vel = {100, 0, 0};  // m/s
    ion.ion_charge_C = 1.602e-19;
    
    // Compute total force
    ForceContext ctx;
    Vec3 F_total = registry.compute_total_force(ion, 0.0, ctx);
    
    // Verify (should be sum of electric + Lorentz)
    REQUIRE(F_total.magnitude() > 0.0);
}
```

### Benchmarks

```cpp
BENCHMARK("ForceRegistry - 10 forces, 1000 ions") {
    // Setup registry with 10 forces
    ForceRegistry registry;
    // ... add forces
    
    // 1000 ions
    std::vector<IonState> ions(1000);
    
    // Compute forces
    ForceContext ctx;
    for (auto& ion : ions) {
        registry.compute_total_force(ion, 0.0, ctx);
    }
};
```

---

## 🎯 SUCCESS CRITERIA

### Functional Requirements
- ✅ All forces implement `IForce` interface
- ✅ `ForceRegistry::compute_total_force()` works correctly
- ✅ Instrument-specific fields match legacy behavior
- ✅ Magnetic force computes Lorentz force correctly
- ✅ Damping force supports friction and Langevin models

### Non-Functional Requirements
- ✅ Unit test coverage >90%
- ✅ Performance regression <5% (compared to legacy)
- ✅ Zero memory leaks (valgrind clean)
- ✅ Clean code (no warnings with `-Wall -Wextra`)

### Integration Requirements
- ✅ Can be used by integrators (Phase 4)
- ✅ Works with existing field providers
- ✅ Compatible with GPU acceleration (future)

---

## 🔄 MIGRATION NOTES

### Legacy Code to Replace
- `src/core/physics/computeAccelerations.cpp` (partially)
- Instrument-specific field calculations (scattered)
- Force computation in `integrate_one_step()` (partially)

### Legacy Code to Keep (for now)
- Space charge solver (Poisson equation)
- Field providers (grid-based, analytical)
- Collision/reaction handling (Phase 2-3)

---

## 📚 REFERENCES

### Legacy Code
- `src/core/physics/computeAccelerations.cpp`
- `src/core/physics/fields/defineFields.cpp`
- `src/core/integrator/integrator_helpers.cpp`

### Related Documentation
- Phase 4: Integration Strategies (needs ForceRegistry)
- Phase 5: SimulationEngine (uses ForceRegistry)

---

**Last Updated:** 2025-11-21  
**Status:** 🔵 Ready to Start  
**Next Step:** Create branch `refactor/force-system`

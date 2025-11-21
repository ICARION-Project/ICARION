# COLLISION SYSTEM REFACTORING PLAN

**Branch:** `refactor/collision-system`  
**Date:** November 21, 2025  
**Goal:** Extract stochastic collision logic into modular `ICollisionHandler` architecture (SSOT-compliant)

---

## 🎯 PROBLEM STATEMENT

### Current Issues:

1. **Inline collision logic** in `integrate_one_step()` (~150 lines of mixed EHSS/HSS code)
2. **Wrong naming:** `HSMC` should be `HSS` (Hard-Sphere Stochastic)
3. **Duplicate damping definitions:** `defineCollisionForces.cpp` duplicates `DampingForce`
4. **No separation:** Deterministic (continuous) vs Stochastic (discrete) collisions mixed

---

## 📋 COLLISION MODEL CLASSIFICATION

### **Deterministic Models** (Continuous Damping Force)

| Model | Implementation | Thermal Noise |
|-------|---------------|---------------|
| **Friction** | `DampingForce` (mobility-based) | Optional: + OU |
| **Langevin** | `DampingForce` (velocity-dependent) | Optional: + OU |
| **HardSphere/HSD** | `DampingForce` (collision frequency) | Optional: + OU |

**→ Already implemented in `src/core/physics/forces/DampingForce.h`!** ✅

**Note:** These models apply **continuous friction forces** F = -γ·m·v

---

### **Stochastic Models** (Discrete Collision Events)

| Model | Collision Type | Geometry | Implementation |
|-------|---------------|----------|----------------|
| **EHSS** | Structure-resolved | Atom-centered spheres | `ICollisionHandler` |
| **HSS** | Isotropic scattering | Single effective sphere | `ICollisionHandler` |
| **OU** | Thermal velocity kicks | N/A (add-on) | `ICollisionHandler` |

**→ NOT yet modularized (inline in `integrator_helpers.cpp`)!** ❌

**Note:** These models apply **discrete velocity changes** during collision events

---

## 🏗️ TARGET ARCHITECTURE

### **Component Hierarchy:**

```
CollisionModel Enum:
  ├─ NoCollisions    (no collisions - vacuum)
  ├─ Friction        (deterministic - uses DampingForce)
  ├─ Langevin        (deterministic - uses DampingForce)
  ├─ HSD             (deterministic - uses DampingForce)
  ├─ EHSS            (stochastic - uses ICollisionHandler)
  ├─ HSS             (stochastic - uses ICollisionHandler)  ✅ NOT HSMC!
  └─ OU              (stochastic kicks - uses ICollisionHandler)

Force System:
  └─ DampingForce (already implemented!)
       ├─ Friction model
       ├─ Langevin model
       └─ HardSphere model

Collision System (NEW):
  └─ ICollisionHandler (interface)
       ├─ EHSSCollisionHandler  (structure-resolved scattering)
       ├─ HSSCollisionHandler   (isotropic scattering)
       └─ OUCollisionHandler    (Ornstein-Uhlenbeck thermal kicks)

Factory:
  └─ CollisionHandlerFactory
       └─ create(config.physics.collision_model)
            ├─ Returns nullptr for deterministic models (use DampingForce)
            └─ Returns handler for stochastic models (EHSS, HSS, OU)
```

---

## 🎯 SSOT DESIGN PRINCIPLES

### ✅ **What We DO:**

1. **Direct config references:** `ICollisionHandler::handle_collision(const EnvironmentConfig& env)`
2. **No parameter duplication:** Read from `env.temperature_K`, `env.particle_density_m_3`, etc.
3. **Factory uses config directly:** `factory.create(config.physics.collision_model)`
4. **Clean integration:** `collision_handler->handle_collision(ion, dt, rng, domain.environment)`

### ❌ **What We DON'T DO:**

1. ❌ No `CollisionContext` struct (parameter duplication!)
2. ❌ No `CollisionConfig` struct (use `PhysicsConfig` directly)
3. ❌ No struct-to-struct conversions
4. ❌ No parameter copies in constructors

---

## 📦 FILE STRUCTURE

```
src/core/physics/collisions/
├── ICollisionHandler.h              ✅ Interface (SSOT-compliant)
├── CollisionStats.h                 ✅ Statistics struct
├── EHSSCollisionHandler.h           ✅ EHSS implementation
├── EHSSCollisionHandler.cpp
├── HSSCollisionHandler.h            ✅ HSS implementation (NOT HSMC!)
├── HSSCollisionHandler.cpp
├── OUCollisionHandler.h             ✅ OU thermal kicks
├── OUCollisionHandler.cpp
├── CollisionHandlerFactory.h        ✅ Factory
└── CollisionHandlerFactory.cpp

tests/physics/collisions/
├── CMakeLists.txt
├── test_ehss_collision_handler.cpp
├── test_hss_collision_handler.cpp
├── test_ou_collision_handler.cpp
└── test_collision_factory.cpp
```

---

## 🚀 IMPLEMENTATION PHASES

### **Phase 2A: Rename HSMC → HSS (1 day)**

**Goal:** Fix wrong naming throughout codebase

**Tasks:**
1. ✅ Update `src/core/config/types/PhysicsEnums.h`: `HSMC` → `HSS`
2. ✅ Add backward-compat alias in `EnumMapper.h`: `"hsmc"` → `CollisionModel::HSS`
3. ✅ Update all references in collision code
4. ✅ Update tests
5. ✅ Update documentation

**Files to Change:**
- `src/core/config/types/PhysicsEnums.h`
- `src/core/config/conversion/EnumMapper.h`
- `src/core/param/paramUtils.h` (legacy - will be removed later)
- `src/core/integrator/integrator_helpers.cpp`
- `tests/**/*.cpp` (grep for HSMC)
- `docs/ARCHITECTURE.md`

**Commit:** `refactor: Rename HSMC → HSS (correct collision model naming)`

---

### **Phase 2B: Delete Duplicate Damping Code (1 day)**

**Goal:** Remove duplicate force definitions (already in `DampingForce`)

**Tasks:**
1. ❌ **DELETE** `src/core/physics/collisions/defineCollisionForces.cpp`
2. ❌ **DELETE** `src/core/physics/collisions/defineCollisionForces.h`
3. ✅ Verify `DampingForce` covers all cases (Friction, Langevin, HardSphere)
4. ✅ Update CMakeLists.txt to remove deleted files
5. ✅ Verify no other files reference deleted headers

**Rationale:** These functions duplicate `DampingForce::compute()` logic (SSOT violation!)

**Commit:** `refactor: Remove duplicate damping force definitions (use DampingForce)`

---

### **Phase 2C: Implement Stochastic Collision Handlers (3 days)**

**Goal:** Create modular `ICollisionHandler` architecture

#### **Day 1: Interface + EHSS**

**Files to Create:**
```cpp
// src/core/physics/collisions/ICollisionHandler.h
class ICollisionHandler {
public:
    virtual bool handle_collision(
        IonState& ion,
        double dt,
        EhssRng& rng,
        const config::EnvironmentConfig& env  // ✅ SSOT!
    ) = 0;
    virtual std::string name() const = 0;
};

// src/core/physics/collisions/EHSSCollisionHandler.h
class EHSSCollisionHandler : public ICollisionHandler {
public:
    explicit EHSSCollisionHandler(
        const GeometryMap& geometry_map,  // Reference, not copy!
        bool enable_logging = false
    );
    
    bool handle_collision(...) override;
    
private:
    const GeometryMap& geometry_map_;  // ✅ SSOT!
    bool enable_logging_;
    CollisionStats stats_;
    
    double compute_effective_ccs(const IonState& ion, double neutral_radius) const;
    Vec3 sample_post_collision_velocity(...) const;
};
```

**Tests:**
- `tests/physics/collisions/test_ehss_collision_handler.cpp`
- Verify SSOT: Pass `EnvironmentConfig` directly, check no parameter copies

**Commit:** `feat: Add ICollisionHandler interface and EHSS implementation`

---

#### **Day 2: HSS + OU Handlers**

**Files to Create:**
```cpp
// src/core/physics/collisions/HSSCollisionHandler.h
class HSSCollisionHandler : public ICollisionHandler {
public:
    explicit HSSCollisionHandler(bool enable_logging = false);
    bool handle_collision(...) override;
    
private:
    bool enable_logging_;
    CollisionStats stats_;
    Vec3 sample_post_collision_velocity(...) const;
};

// src/core/physics/collisions/OUCollisionHandler.h
class OUCollisionHandler : public ICollisionHandler {
public:
    explicit OUCollisionHandler(double gamma_coefficient);
    bool handle_collision(...) override;
    
private:
    double gamma_;  // Must match DampingForce gamma!
};
```

**Key Points:**
- **HSS:** Uses `ion.CCS_m2` (effective cross-section), isotropic scattering
- **OU:** Calls existing `apply_ou_velocity_kick()`, reads `env.temperature_K`

**Tests:**
- `tests/physics/collisions/test_hss_collision_handler.cpp`
- `tests/physics/collisions/test_ou_collision_handler.cpp`

**Commit:** `feat: Add HSS and OU collision handlers`

---

#### **Day 3: Factory**

**Files to Create:**
```cpp
// src/core/physics/collisions/CollisionHandlerFactory.h
class CollisionHandlerFactory {
public:
    static std::unique_ptr<ICollisionHandler> create(
        const config::PhysicsConfig& config,  // ✅ SSOT!
        const GeometryMap* geometry_map = nullptr,
        double gamma_for_ou = 0.0
    );
};

// Implementation:
std::unique_ptr<ICollisionHandler> CollisionHandlerFactory::create(...) {
    switch (config.collision_model) {
        case CollisionModel::EHSS:
            return std::make_unique<EHSSCollisionHandler>(*geometry_map, ...);
        case CollisionModel::HSS:
            return std::make_unique<HSSCollisionHandler>(...);
        case CollisionModel::OU:
            return std::make_unique<OUCollisionHandler>(gamma_for_ou);
        case CollisionModel::Friction:
        case CollisionModel::Langevin:
        case CollisionModel::HSD:
        case CollisionModel::NoCollisions:
            return nullptr;  // ✅ Use DampingForce!
    }
}
```

**Tests:**
- `tests/physics/collisions/test_collision_factory.cpp`
- Verify factory returns correct handlers
- Verify factory returns `nullptr` for deterministic models

**Commit:** `feat: Add CollisionHandlerFactory with SSOT design`

---

### **Phase 2D: Integrate into Integrator (2 days)**

**Goal:** Replace inline collision logic with `ICollisionHandler` calls

#### **Day 1: Refactor `integrate_one_step()`**

**Current Code (inline logic):**
```cpp
// src/core/integrator/integrator_helpers.cpp::integrate_one_step()
if (gParams.collisionModel == CollisionModel::EHSS) {
    // 50 lines of EHSS logic
    auto it = geometry_map.find(y.species_id);
    y.vel = collide_ehss_cpu_geometry_given_neutral(...);
} else if (gParams.collisionModel == CollisionModel::HSMC) {
    // 20 lines of HSS logic
    y.vel = collide_hs_cpu(...);
}
```

**Refactored Code:**
```cpp
// New signature:
void integrate_one_step(
    IonState& ion,
    const config::FullConfig& config,           // ✅ SSOT!
    const config::DomainConfig& domain,         // ✅ SSOT!
    physics::ICollisionHandler* collision_handler,  // ✅ Injected!
    physics::ForceRegistry* force_registry,
    EhssRng& local_rng
) {
    const double dt = config.simulation.dt_s;
    
    // === 1. Stochastic Collision ===
    if (collision_handler) {
        collision_handler->handle_collision(
            ion, dt, local_rng, 
            domain.environment  // ✅ SSOT!
        );
    }
    
    // === 2. Force Computation (includes DampingForce) ===
    physics::ForceContext force_ctx{
        .domain = &domain,
        .temperature_K = domain.environment.temperature_K,
        .pressure_Pa = domain.environment.pressure_Pa,
        // ...
    };
    Vec3 total_force = force_registry->compute_total_force(ion, t, force_ctx);
    
    // === 3. Integration ===
    // ... RK4/RK45 step ...
}
```

**Changes:**
1. ✅ Remove all inline collision logic
2. ✅ Add `ICollisionHandler*` parameter (nullable!)
3. ✅ Pass `domain.environment` directly (SSOT!)
4. ✅ Update call sites to create handler via factory

**Commit:** `refactor: Replace inline collision logic with ICollisionHandler`

---

#### **Day 2: Integration Tests**

**Tests to Add:**
- `tests/integration/test_collision_integration.cpp`
- Verify EHSS + ElectricFieldForce works
- Verify HSS + MagneticFieldForce works
- Verify Langevin (DampingForce) + OU (handler) works
- Verify deterministic models (no handler) work

**Commit:** `test: Add integration tests for collision system`

---

## ✅ VALIDATION CHECKLIST

### **SSOT Compliance:**
- [ ] All handlers read from `const EnvironmentConfig&` (no copies)
- [ ] Factory uses `config.physics.collision_model` (direct access)
- [ ] No `CollisionContext` or `CollisionConfig` structs
- [ ] No parameter duplication in constructors

### **Functionality:**
- [ ] EHSS: Structure-resolved scattering works
- [ ] HSS: Isotropic scattering works
- [ ] OU: Thermal kicks work (with DampingForce)
- [ ] Deterministic models: Use DampingForce (no handler)
- [ ] Factory: Returns correct handlers for each model
- [ ] Integration: `integrate_one_step()` works with all models

### **Tests:**
- [ ] Unit tests for each handler (EHSS, HSS, OU)
- [ ] Unit tests for factory
- [ ] Integration tests for collision + forces
- [ ] All tests pass (100% coverage)

### **Documentation:**
- [ ] Update `docs/ARCHITECTURE.md` (collision system section)
- [ ] Update `docs/DEVELOPERS_GUIDE.md` (adding collision handlers)
- [ ] Add SSOT warnings where applicable

---

## 📊 METRICS

**Before:**
- Inline collision logic: ~150 lines in `integrator_helpers.cpp`
- Duplicate damping code: ~200 lines in `defineCollisionForces.cpp`
- SSOT violations: 2 (parameter structs, duplicate forces)

**After:**
- Modular handlers: ~500 lines (across 6 files)
- Inline logic: 0 lines (removed!)
- Duplicate code: 0 lines (removed!)
- SSOT violations: 0 ✅

**Test Coverage:**
- Unit tests: 15+ test cases (3 handlers + factory)
- Integration tests: 5+ test cases
- Total assertions: ~80+

---

## 🎯 SUCCESS CRITERIA

1. ✅ All collision logic extracted into `ICollisionHandler` implementations
2. ✅ SSOT: All handlers read from `EnvironmentConfig` (no copies)
3. ✅ Factory returns correct handlers (or nullptr for deterministic models)
4. ✅ `integrate_one_step()` uses injected handler (no inline logic)
5. ✅ All tests pass (100% green)
6. ✅ Documentation updated

---

## 🚀 NEXT PHASE

**Phase 3:** Reaction System Refactoring
- Extract reaction logic into `IReactionHandler`
- Similar SSOT pattern (read from `PhysicsConfig`)
- Modular reaction types (ion-neutral, charge transfer, etc.)

---

**Start Date:** November 21, 2025  
**Estimated Completion:** November 28, 2025 (7 days)  
**Branch:** `refactor/collision-system`  
**Merge Target:** `integrator-refactor`

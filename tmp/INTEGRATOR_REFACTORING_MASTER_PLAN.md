# 🔧 INTEGRATOR REFACTORING - MASTER PLAN

**Project:** ICARION Core Integrator Modernization  
**Branch:** `integrator-refactor` (base)  
**Start Date:** 2025-11-21  
**Estimated Duration:** 5 weeks  
**Status:** 🟢 Planning Complete - Ready to Start

---

## 📋 EXECUTIVE SUMMARY

### Current Problems
- ❌ **Monolithic functions** (260+ lines)
- ❌ **12+ parameter functions**
- ❌ **Mixed responsibilities** (integration + I/O + physics)
- ❌ **Global state** (static counters)
- ❌ **Poor testability** (can't isolate components)
- ❌ **Instrument-specific if/else chains**
- ❌ **No abstractions** (collision/reaction handlers inline)

### Solution: Modular Architecture
- ✅ **Single Responsibility Principle**
- ✅ **Dependency Injection**
- ✅ **Interface-Based Design**
- ✅ **Unit Testable Components**
- ✅ **Modern C++17** features
- ✅ **Zero-cost abstractions**

---

## 🏗️ TARGET ARCHITECTURE

```
SimulationEngine (orchestrator)
    │
    ├── IntegrationStrategy (solver selection)
    │   ├── RK4Integrator
    │   ├── RK45Integrator
    │   └── BorisIntegrator
    │
    ├── ForceRegistry (computes total force)
    │   ├── ElectricFieldForce
    │   ├── MagneticFieldForce
    │   ├── SpaceChargeForce
    │   └── DampingForce
    │
    ├── CollisionHandler (interface)
    │   ├── EHSSCollisionHandler
    │   ├── HSMCCollisionHandler
    │   ├── LangevinCollisionHandler
    │   └── FrictionCollisionHandler
    │
    ├── ReactionEngine
    │   ├── ReactionRateCalculator
    │   └── ReactionEventQueue
    │
    ├── DomainManager
    │   ├── BoundaryConditions
    │   └── DomainTransition
    │
    └── OutputManager
        ├── TrajectoryWriter (HDF5Writer v2)
        └── StatisticsCollector
```

---

## 📅 PHASE BREAKDOWN

### **PHASE 1: Force System Refactoring** 
**Branch:** `refactor/force-system`  
**Duration:** 1 week  
**Status:** 🔵 Ready to Start

#### Goals
- Modular force computation with clean interfaces
- Eliminate instrument-specific if/else chains
- Enable force composition (add/remove forces dynamically)
- Full unit test coverage

#### Deliverables
1. `IForce` interface (`src/core/physics/forces/IForce.h`)
2. `ElectricFieldForce` implementation
3. `MagneticFieldForce` implementation
4. `DampingForce` implementation (friction/Langevin)
5. `SpaceChargeForce` implementation
6. `ForceRegistry` manager
7. Unit tests (>90% coverage)

#### Files to Create
```
src/core/physics/forces/
├── IForce.h                    (interface)
├── ForceContext.h              (context struct)
├── ElectricFieldForce.h        (header)
├── ElectricFieldForce.cpp      (impl)
├── MagneticFieldForce.h
├── MagneticFieldForce.cpp
├── DampingForce.h
├── DampingForce.cpp
├── SpaceChargeForce.h
├── SpaceChargeForce.cpp
├── ForceRegistry.h             (manager)
└── ForceRegistry.cpp

tests/physics/forces/
├── CMakeLists.txt
├── test_electric_field_force.cpp
├── test_magnetic_field_force.cpp
├── test_damping_force.cpp
└── test_force_registry.cpp
```

#### Success Criteria
- ✅ All forces implement `IForce` interface
- ✅ `ForceRegistry::compute_total_force()` works correctly
- ✅ Unit tests pass (90%+ coverage)
- ✅ Benchmarks show <5% performance regression
- ✅ Instrument-specific code isolated in force implementations

---

### **PHASE 2: Collision System Refactoring**
**Branch:** `refactor/collision-system`  
**Duration:** 1 week  
**Status:** ⏳ Waiting for Phase 1

#### Goals
- Interface-based collision handlers
- Remove collision logic from integrator
- Enable collision model swapping at runtime
- Separate physics from debug logging

#### Deliverables
1. `ICollisionHandler` interface
2. `EHSSCollisionHandler` implementation
3. `HSMCCollisionHandler` implementation
4. `LangevinCollisionHandler` implementation
5. `FrictionCollisionHandler` implementation
6. `NoCollisionHandler` (null object pattern)
7. `CollisionHandlerFactory`
8. Unit tests

#### Files to Create
```
src/core/physics/collisions/
├── ICollisionHandler.h         (interface)
├── CollisionContext.h          (context struct)
├── CollisionStats.h            (statistics)
├── EHSSCollisionHandler.h
├── EHSSCollisionHandler.cpp
├── HSMCCollisionHandler.h
├── HSMCCollisionHandler.cpp
├── LangevinCollisionHandler.h
├── LangevinCollisionHandler.cpp
├── FrictionCollisionHandler.h
├── FrictionCollisionHandler.cpp
├── NoCollisionHandler.h
├── CollisionHandlerFactory.h
└── CollisionHandlerFactory.cpp

tests/physics/collisions/
├── CMakeLists.txt
├── test_ehss_collision_handler.cpp
├── test_hsmc_collision_handler.cpp
├── test_langevin_collision_handler.cpp
└── test_collision_factory.cpp
```

#### Success Criteria
- ✅ All collision models implement `ICollisionHandler`
- ✅ Factory creates correct handler based on config
- ✅ Statistics tracking works
- ✅ Unit tests pass (90%+ coverage)
- ✅ No I/O in physics code (logging separated)

---

### **PHASE 3: Reaction System Refactoring**
**Branch:** `refactor/reaction-system`  
**Duration:** 3 days  
**Status:** ⏳ Waiting for Phase 1

#### Goals
- Efficient reaction engine with fast species lookup
- Event tracking for reaction history
- Separate rate calculation from event handling
- Optimized data structures (hash maps)

#### Deliverables
1. `ReactionEngine` class
2. `ReactionRateCalculator`
3. `ReactionEvent` struct (tracking)
4. Optimized species lookup (std::unordered_map)
5. Unit tests

#### Files to Create
```
src/core/physics/reactions/
├── ReactionEngine.h
├── ReactionEngine.cpp
├── ReactionRateCalculator.h
├── ReactionRateCalculator.cpp
├── ReactionEvent.h
├── ReactionContext.h
└── ReactionStats.h

tests/physics/reactions/
├── CMakeLists.txt
├── test_reaction_engine.cpp
└── test_reaction_rate_calculator.cpp
```

#### Success Criteria
- ✅ Reaction lookup O(1) (hash map)
- ✅ Event tracking works
- ✅ Unit tests pass (90%+ coverage)
- ✅ Performance better than legacy system

---

### **PHASE 4: Integration Strategy Refactoring**
**Branch:** `refactor/integration-strategies`  
**Duration:** 1 week  
**Status:** ⏳ Waiting for Phase 1 (needs ForceRegistry)

#### Goals
- Strategy pattern for different solvers
- Clean RK4, RK45, Boris implementations
- Separate adaptive logic from integration
- Enable solver swapping at runtime

#### Deliverables
1. `IIntegrationStrategy` interface
2. `RK4Integrator` implementation
3. `RK45Integrator` implementation (adaptive)
4. `BorisIntegrator` implementation (magnetic fields)
5. `IntegrationStrategyFactory`
6. Unit tests

#### Files to Create
```
src/core/integrator/strategies/
├── IIntegrationStrategy.h      (interface)
├── DerivativeFunction.h        (type alias)
├── RK4Integrator.h
├── RK4Integrator.cpp
├── RK45Integrator.h
├── RK45Integrator.cpp
├── BorisIntegrator.h
├── BorisIntegrator.cpp
├── IntegrationStrategyFactory.h
└── IntegrationStrategyFactory.cpp

tests/integrator/strategies/
├── CMakeLists.txt
├── test_rk4_integrator.cpp
├── test_rk45_integrator.cpp
├── test_boris_integrator.cpp
└── test_strategy_factory.cpp
```

#### Success Criteria
- ✅ All solvers implement `IIntegrationStrategy`
- ✅ RK45 adaptive stepping works correctly
- ✅ Boris integrator preserves energy
- ✅ Unit tests pass (95%+ coverage)
- ✅ Performance matches legacy implementations

---

### **PHASE 5: SimulationEngine**
**Branch:** `refactor/simulation-engine`  
**Duration:** 1 week  
**Status:** ⏳ Waiting for Phase 1-4

#### Goals
- Clean orchestrator with dependency injection
- Replace monolithic `integrate_trajectory()`
- Modular design (swap components easily)
- Separate concerns (physics / I/O / domain management)

#### Deliverables
1. `SimulationEngine` class
2. `OutputManager` (wraps HDF5Writer v2)
3. `DomainManager` (boundary conditions, transitions)
4. Integration tests
5. Performance benchmarks

#### Files to Create
```
src/core/integrator/
├── SimulationEngine.h
├── SimulationEngine.cpp
├── OutputManager.h
├── OutputManager.cpp
├── DomainManager.h
├── DomainManager.cpp
└── SimulationStats.h

tests/integrator/
├── CMakeLists.txt
├── test_simulation_engine.cpp
├── test_output_manager.cpp
└── test_domain_manager.cpp
```

#### Success Criteria
- ✅ `SimulationEngine::run()` works end-to-end
- ✅ All physics modules injected (no globals)
- ✅ HDF5Writer v2 integrated
- ✅ Integration tests pass
- ✅ Performance within 5% of legacy

---

### **PHASE 6: Legacy Migration & Cleanup**
**Branch:** `refactor/legacy-cleanup`  
**Duration:** 3 days  
**Status:** ⏳ Waiting for Phase 5

#### Goals
- Update `src/main/main.cpp` to use new engine
- Delete old integrator code
- Regression tests (ensure no breakage)
- Performance validation

#### Tasks
1. Update `main.cpp` to instantiate `SimulationEngine`
2. Delete old `integrate_trajectory()` function
3. Delete old `integrate_one_step()` function
4. Delete old collision/reaction handling code
5. Run all validation configs
6. Performance benchmarks (CPU/GPU parity)

#### Success Criteria
- ✅ All example configs run successfully
- ✅ Validation tests pass (CPU/GPU parity)
- ✅ Performance within 5% of baseline
- ✅ No memory leaks (valgrind clean)
- ✅ Code coverage >85%

---

## 🎯 IMPLEMENTATION STRATEGY

### Branch Structure
```
integrator-refactor (base)
    ├── refactor/force-system           (Phase 1)
    ├── refactor/collision-system       (Phase 2)
    ├── refactor/reaction-system        (Phase 3)
    ├── refactor/integration-strategies (Phase 4)
    ├── refactor/simulation-engine      (Phase 5)
    └── refactor/legacy-cleanup         (Phase 6)
```

### Workflow
1. **Branch out** from `integrator-refactor` for each phase
2. **Implement** features with TDD (test-first)
3. **Review** code and tests
4. **Merge** back to `integrator-refactor` when complete
5. **Repeat** for next phase

### Testing Strategy
- **Unit tests** for each component (90%+ coverage)
- **Integration tests** for SimulationEngine
- **Regression tests** against legacy integrator
- **Performance benchmarks** (ensure <5% overhead)

---

## 📊 PROGRESS TRACKING

### Phase Status
| Phase | Component | Status | Branch | Duration | Completion |
|-------|-----------|--------|--------|----------|-----------|
| 1 | Force System | 🔵 Ready | `refactor/force-system` | 1 week | 0% |
| 2 | Collision System | ⏳ Waiting | `refactor/collision-system` | 1 week | 0% |
| 3 | Reaction System | ⏳ Waiting | `refactor/reaction-system` | 3 days | 0% |
| 4 | Integration Strategies | ⏳ Waiting | `refactor/integration-strategies` | 1 week | 0% |
| 5 | SimulationEngine | ⏳ Waiting | `refactor/simulation-engine` | 1 week | 0% |
| 6 | Legacy Cleanup | ⏳ Waiting | `refactor/legacy-cleanup` | 3 days | 0% |

### Overall Progress
- **Total Duration:** ~5 weeks
- **Completion:** 0% (planning complete)
- **Next Step:** Start Phase 1 (Force System)

---

## 🚀 NEXT STEPS

1. ✅ **Create branch** `refactor/force-system`
2. ✅ **Implement** `IForce` interface
3. ✅ **Implement** `ElectricFieldForce`
4. ✅ **Implement** `ForceRegistry`
5. ✅ **Write tests** (TDD approach)
6. 🔄 **Merge** when Phase 1 complete

---

## 📚 REFERENCES

- **Legacy Code:**
  - `src/core/integrator/integrator.cpp` (260 lines)
  - `src/core/integrator/integrator_helpers.cpp` (230 lines)
  - `src/core/physics/computeAccelerations.cpp`
  - `src/core/physics/collisions/`

- **Related Work:**
  - HDF5Writer v2 refactoring (completed)
  - Config system refactoring (completed)
  - SHA256 hashing (completed)

- **Documentation:**
  - See `tmp/INTEGRATOR_REFACTORING_PHASE1_FORCE_SYSTEM.md`
  - See `tmp/INTEGRATOR_REFACTORING_PHASE2_COLLISION_SYSTEM.md`
  - See `tmp/INTEGRATOR_REFACTORING_PHASE3_REACTION_SYSTEM.md`
  - See `tmp/INTEGRATOR_REFACTORING_PHASE4_INTEGRATION_STRATEGIES.md`
  - See `tmp/INTEGRATOR_REFACTORING_PHASE5_SIMULATION_ENGINE.md`

---

**Last Updated:** 2025-11-21  
**Maintained By:** ICARION Development Team  
**Status:** 🟢 Active Development

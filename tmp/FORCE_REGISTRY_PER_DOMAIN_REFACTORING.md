# ForceRegistry per Domain Refactoring Plan

## Motivation

ChatGPT Feedback (Phase 12):
> "ForceRegistry sollte DomainConfig nicht in step() bekommen, sondern im Konstruktor.
> Das würde RK4/RK45 weiter vereinfachen."

**Problem:** Multi-Domain Simulation
- Aktuell: 1 ForceRegistry für alle Domains
- Domain wechselt pro Ion → domain muss Parameter bleiben
- Lösung: **1 ForceRegistry pro Domain**

## Benefits

### 1. Simplified Integrator API
```cpp
// Before
integrator->step(ion, t, dt, force_registry, domain, all_ions);

// After  
integrator->step(ion, t, dt, force_registry, all_ions);
```

### 2. Better Encapsulation
- ForceRegistry owns its domain configuration
- No need to pass domain through multiple layers
- Integrator doesn't need domain knowledge

### 3. Domain-Specific Forces
Each domain can have different forces:
- Domain 0 (LQIT): Electric field + RF field + damping
- Domain 1 (IMS): Electric field + drag force
- Domain 2 (TOF): Electric field only

### 4. Thread-Safe by Design
Each domain has its own registry → no shared state in parallel loops

## Implementation Steps

### Step 1: Update ForceRegistry Constructor

**File:** `src/core/physics/forces/ForceRegistry.h`

```cpp
class ForceRegistry {
public:
    // Remove default constructor
    // explicit ForceRegistry(const config::DomainConfig& domain);
    
    // Store domain as member
    const config::DomainConfig& domain() const { return domain_; }
    
private:
    const config::DomainConfig& domain_;  // Store reference
    std::vector<std::unique_ptr<IForce>> forces_;
};
```

**Changes:**
- ✅ Add domain_ member (const reference)
- ✅ Constructor requires domain parameter
- ✅ Add getter domain()
- ❌ Remove default constructor

### Step 2: Update IIntegrationStrategy Interface

**File:** `src/core/integrator/strategies/IIntegrationStrategy.h`

```cpp
class IIntegrationStrategy {
public:
    virtual void step(
        IonState& ion,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry,
        // REMOVED: const config::DomainConfig& domain,
        const std::vector<IonState>& all_ions
    ) = 0;
};
```

**Changes:**
- ❌ Remove domain parameter from step()
- ✅ Access via force_registry.domain() if needed
- ✅ Update all implementations (RK4, RK45, Boris)

### Step 3: Update Integration Strategies

**Files to modify:**
- `src/core/integrator/strategies/RK4Strategy.h/cpp`
- `src/core/integrator/strategies/RK45Strategy.h/cpp`
- `src/core/integrator/strategies/BorisStrategy.h/cpp`

**Example (RK4Strategy.cpp):**

```cpp
void RK4Strategy::step(
    IonState& ion,
    double t,
    double dt,
    const physics::ForceRegistry& force_registry,
    // REMOVED: const config::DomainConfig& domain,
    const std::vector<IonState>& all_ions
) {
    // Get domain from registry
    const auto& domain = force_registry.domain();
    
    // Create context
    physics::ForceContext ctx;
    ctx.domain = &domain;  // Get from registry
    ctx.all_ions = &all_ions;
    
    // Rest unchanged...
}
```

**Changes:**
- ❌ Remove domain parameter from step()
- ✅ Get domain via force_registry.domain()
- ✅ Update ForceContext initialization

### Step 4: Update SimulationEngine

**File:** `src/core/integrator/SimulationEngine.h`

```cpp
class SimulationEngine {
public:
    SimulationEngine(
        const config::FullConfig& config,
        // REMOVED: std::shared_ptr<physics::ForceRegistry> force_registry,
        std::shared_ptr<IIntegrationStrategy> integrator,
        // ...
    );
    
private:
    // NEW: One registry per domain
    std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries_;
    // OLD: Single registry
    // std::shared_ptr<physics::ForceRegistry> force_registry_;
};
```

**Changes:**
- ✅ Replace single force_registry_ with vector force_registries_
- ✅ Create one registry per domain in constructor
- ❌ Remove force_registry parameter from constructor

**File:** `src/core/integrator/SimulationEngine.cpp`

#### Constructor Changes:

```cpp
SimulationEngine::SimulationEngine(
    const config::FullConfig& config,
    // REMOVED: std::shared_ptr<physics::ForceRegistry> force_registry,
    std::shared_ptr<IIntegrationStrategy> integrator,
    std::shared_ptr<physics::ICollisionHandler> collision_handler,
    std::shared_ptr<physics::IReactionHandler> reaction_handler
) : config_(config),
    integrator_(integrator),
    collision_handler_(collision_handler),
    reaction_handler_(reaction_handler)
{
    // Create force registries (one per domain)
    force_registries_.reserve(config_.domains.size());
    for (const auto& domain : config_.domains) {
        auto registry = std::make_shared<physics::ForceRegistry>(domain);
        
        // TODO: Populate forces based on domain configuration
        // This is where we'd add domain-specific forces
        // Example:
        // if (domain.fields.has_electric_field) {
        //     registry->add_force(std::make_unique<ElectricFieldForce>(domain));
        // }
        
        force_registries_.push_back(registry);
    }
    
    // Validation
    if (!integrator_) {
        throw std::invalid_argument("SimulationEngine: IntegrationStrategy cannot be null");
    }
    if (force_registries_.empty()) {
        throw std::invalid_argument("SimulationEngine: No domains configured");
    }
    
    // Rest of initialization...
}
```

#### process_timestep() Changes:

```cpp
void SimulationEngine::process_timestep(std::vector<IonState>& ions, double dt) {
    const int n_ions = static_cast<int>(ions.size());
    
    // Ion-based RNG (already implemented)
    std::vector<EhssRng> rng_by_ion;
    rng_by_ion.reserve(n_ions);
    for (int i = 0; i < n_ions; ++i) {
        uint64_t ion_seed = config_.simulation.rng_seed + static_cast<uint64_t>(i);
        rng_by_ion.emplace_back(ion_seed);
    }
    
    #pragma omp parallel if(config_.simulation.enable_openmp)
    {
        #pragma omp for schedule(dynamic)
        for (int i = 0; i < n_ions; ++i) {
            IonState& ion = ions[i];
            EhssRng& ion_rng = rng_by_ion[i];
            
            if (!ion.active) continue;
            
            // 1. Find current domain
            int domain_idx = domain_manager_->find_domain_index(ion.pos);
            if (domain_idx < 0) {
                ion.active = false;
                continue;
            }
            
            // 2. Get domain-specific force registry
            auto& force_registry = *force_registries_[domain_idx];
            const auto& domain_config = config_.domains[domain_idx];
            
            // 3. Update domain properties
            if (ion.current_domain_index != domain_idx) {
                domain_manager_->update_domain_properties(ion, domain_idx);
            }
            
            // 4. Transform to local coordinates
            Vec3 pos_local = domain_manager_->global_to_local_pos(ion.pos, domain_idx);
            Vec3 vel_local = domain_manager_->global_to_local_vel(ion.vel, domain_idx);
            Vec3 pos_before = pos_local;
            
            // 5. Handle collisions
            if (collision_handler_) {
                IonState ion_local = ion;
                ion_local.pos = pos_local;
                ion_local.vel = vel_local;
                collision_handler_->handle_collision(ion_local, dt, ion_rng, domain_config.environment);
                pos_local = ion_local.pos;
                vel_local = ion_local.vel;
            }
            
            // 6. Handle reactions
            if (reaction_handler_ && !config_.reaction_db.reactions.empty()) {
                IonState ion_local = ion;
                ion_local.pos = pos_local;
                ion_local.vel = vel_local;
                bool reaction_occurred = reaction_handler_->handle_reaction(
                    ion_local, dt, ion_rng,
                    config_.reaction_db,
                    config_.species_db,
                    domain_config.environment
                );
                if (reaction_occurred) {
                    ion.species_id = ion_local.species_id;
                    ion.mass_kg = ion_local.mass_kg;
                    ion.ion_charge_C = ion_local.ion_charge_C;
                }
                pos_local = ion_local.pos;
                vel_local = ion_local.vel;
            }
            
            // 7. Integrate trajectory
            IonState ion_local = ion;
            ion_local.pos = pos_local;
            ion_local.vel = vel_local;
            
            // NEW: No domain parameter!
            integrator_->step(ion_local, current_time_, dt, force_registry, ions);
            // OLD: integrator_->step(ion_local, current_time_, dt, force_registry, domain_config, ions);
            
            pos_local = ion_local.pos;
            vel_local = ion_local.vel;
            
            // 8. Boundary checks (unchanged)
            // ...
        }
    }
}
```

**Key Changes:**
- ✅ Get domain-specific registry: `force_registries_[domain_idx]`
- ✅ Pass registry to integrator (no domain parameter)
- ✅ ForceRegistry already knows its domain

### Step 5: Update main.cpp

**File:** `src/main/main.cpp`

**Before:**
```cpp
// Create ForceRegistry (empty, will be populated by SimulationEngine)
auto force_registry = std::make_shared<physics::ForceRegistry>();

// Create SimulationEngine
auto engine = std::make_unique<integrator::SimulationEngine>(
    full_config,
    force_registry,  // ← Remove this
    integrator_strategy,
    collision_handler,
    reaction_handler
);
```

**After:**
```cpp
// ForceRegistry creation moved into SimulationEngine constructor
// No need to create it here anymore

// Create SimulationEngine
auto engine = std::make_unique<integrator::SimulationEngine>(
    full_config,
    // force_registry removed - created internally per domain
    integrator_strategy,
    collision_handler,
    reaction_handler
);
```

**Changes:**
- ❌ Remove ForceRegistry creation
- ❌ Remove force_registry parameter from SimulationEngine constructor
- ✅ SimulationEngine creates registries internally

### Step 6: Update Tests

**Files to update:**
- `tests/integrator/test_rk4_strategy.cpp`
- `tests/integrator/test_rk45_strategy.cpp`
- `tests/integrator/test_boris_strategy.cpp`
- `tests/integrator/test_domain_manager.cpp`
- `tests/integrator/test_simulation_engine.cpp`

**Example (test_rk4_strategy.cpp):**

**Before:**
```cpp
TEST_CASE("RK4 integration", "[integrator]") {
    config::DomainConfig domain = create_test_domain();
    physics::ForceRegistry registry;
    registry.add_force(std::make_unique<TestForce>());
    
    RK4Strategy rk4;
    IonState ion = create_test_ion();
    
    rk4.step(ion, 0.0, 1e-6, registry, domain, {});
}
```

**After:**
```cpp
TEST_CASE("RK4 integration", "[integrator]") {
    config::DomainConfig domain = create_test_domain();
    physics::ForceRegistry registry(domain);  // ← Pass domain in constructor
    registry.add_force(std::make_unique<TestForce>());
    
    RK4Strategy rk4;
    IonState ion = create_test_ion();
    
    rk4.step(ion, 0.0, 1e-6, registry, {});  // ← No domain parameter
}
```

**Changes:**
- ✅ Pass domain to ForceRegistry constructor
- ❌ Remove domain parameter from step()
- ✅ Update all test cases

### Step 7: Force Factory Pattern (Optional Enhancement)

**New File:** `src/core/physics/forces/ForceFactory.h`

```cpp
namespace ICARION::physics {

class ForceFactory {
public:
    /**
     * @brief Create forces for a domain based on configuration
     * @param domain Domain configuration
     * @return ForceRegistry populated with appropriate forces
     */
    static std::shared_ptr<ForceRegistry> create_for_domain(
        const config::DomainConfig& domain
    ) {
        auto registry = std::make_shared<ForceRegistry>(domain);
        
        // Add forces based on domain configuration
        if (domain.fields.electric_field.enabled) {
            registry->add_force(std::make_unique<ElectricFieldForce>(domain));
        }
        
        if (domain.fields.magnetic_field.enabled) {
            registry->add_force(std::make_unique<MagneticFieldForce>(domain.fields.magnetic_field));
        }
        
        if (domain.environment.enable_damping) {
            double gamma = calculate_damping_coefficient(domain.environment);
            registry->add_force(std::make_unique<DampingForce>(gamma));
        }
        
        if (domain.space_charge.enabled) {
            registry->add_force(std::make_unique<SpaceChargeForce>(
                domain.space_charge.grid_nx,
                domain.space_charge.grid_ny,
                domain.space_charge.grid_nz
            ));
        }
        
        return registry;
    }
};

} // namespace ICARION::physics
```

**Usage in SimulationEngine:**
```cpp
// Constructor
for (const auto& domain : config_.domains) {
    auto registry = ForceFactory::create_for_domain(domain);
    force_registries_.push_back(registry);
}
```

**Benefits:**
- Centralized force creation logic
- Easy to add new force types
- Domain-specific force configurations

## Testing Strategy

### Unit Tests

1. **ForceRegistry Constructor Test**
   - Verify domain is stored correctly
   - Verify domain() getter returns correct reference
   
2. **Integration Strategy Tests**
   - Update all strategy tests to pass domain in ForceRegistry constructor
   - Verify step() works without domain parameter
   - Test with different domain configurations

3. **SimulationEngine Tests**
   - Test multi-domain scenario
   - Verify correct registry is used per domain
   - Test domain transitions (ion crossing boundaries)

### Integration Tests

1. **Multi-Domain Simulation**
   - 2-3 domains with different configurations
   - Verify ions transition correctly
   - Check forces are applied from correct domain registry

2. **Performance Test**
   - Compare performance before/after refactoring
   - Should be same or slightly better (less parameter passing)

3. **Reproducibility Test**
   - Run same simulation with different OpenMP thread counts
   - Verify identical results (ion-based RNG ensures this)

## Migration Path

### Phase 1: Preparation (Safe)
- ✅ Implement ion-based RNG (already done!)
- ✅ Add ForceRegistry constructor overload (backward compatible)
- ✅ Add domain() getter to ForceRegistry

### Phase 2: Parallel Implementation (Safe)
- ✅ Create new IIntegrationStrategy2 interface (without domain parameter)
- ✅ Implement RK4Strategy2, RK45Strategy2, BorisStrategy2
- ✅ Test new strategies with ForceRegistry-per-domain
- ✅ Keep old interface working

### Phase 3: Migration (Breaking)
- ❌ Replace IIntegrationStrategy with new interface
- ❌ Update SimulationEngine to use force_registries_ vector
- ❌ Update main.cpp and tests
- ❌ Remove old constructors/interfaces

### Phase 4: Cleanup
- ❌ Remove deprecated code
- ❌ Update documentation
- ❌ Add ForceFactory for convenience

## Estimated Effort

- **Step 1-2 (Interface):** 30 min
- **Step 3 (Strategies):** 45 min
- **Step 4 (SimulationEngine):** 1-2 hours
- **Step 5 (main.cpp):** 15 min
- **Step 6 (Tests):** 1-2 hours
- **Step 7 (ForceFactory):** 30 min (optional)

**Total:** ~4-6 hours

## Risks & Mitigations

### Risk 1: Breaking Changes
- **Impact:** All existing code using SimulationEngine breaks
- **Mitigation:** Use parallel implementation (Phase 2) for gradual migration
- **Alternative:** Implement in major version update (v2.0)

### Risk 2: Performance Regression
- **Impact:** Multiple ForceRegistry instances might be slower
- **Mitigation:** Benchmark before/after
- **Expected:** Should be same or faster (less parameter passing, better cache locality)

### Risk 3: Test Coverage
- **Impact:** Missing edge cases in multi-domain scenarios
- **Mitigation:** Comprehensive integration tests with 2-3 domains
- **Validation:** Run existing example configs and verify output unchanged

## Decision

**Recommended:** ✅ Implement in Phase 13 (next major release)

**Reasoning:**
- Clean design with clear benefits
- Reasonable effort (~4-6 hours)
- Natural fit for multi-domain architecture
- Simplifies integrator interface
- Enables domain-specific force configurations

**Not urgent because:**
- Current design works fine
- No performance issues
- Ion-based RNG (Phase 12) already provides main reproducibility benefit

## Implementation Checklist

When implementing:

- [ ] Step 1: Update ForceRegistry with domain constructor
- [ ] Step 2: Update IIntegrationStrategy interface
- [ ] Step 3: Update RK4Strategy
- [ ] Step 4: Update RK45Strategy
- [ ] Step 5: Update BorisStrategy
- [ ] Step 6: Update SimulationEngine.h (add force_registries_)
- [ ] Step 7: Update SimulationEngine.cpp constructor
- [ ] Step 8: Update SimulationEngine::process_timestep()
- [ ] Step 9: Update main.cpp
- [ ] Step 10: Update test_rk4_strategy.cpp
- [ ] Step 11: Update test_rk45_strategy.cpp
- [ ] Step 12: Update test_boris_strategy.cpp
- [ ] Step 13: Update test_simulation_engine.cpp
- [ ] Step 14: Update test_domain_manager.cpp
- [ ] Step 15: Run full test suite (30/30 passing)
- [ ] Step 16: Run performance benchmark
- [ ] Step 17: Update documentation
- [ ] Step 18: (Optional) Implement ForceFactory
- [ ] Step 19: Commit and merge to core-dev

## References

- ChatGPT Feedback: Phase 12 Review
- Current Implementation: `src/core/integrator/SimulationEngine.cpp`
- Force System: `src/core/physics/forces/`
- Integration Strategies: `src/core/integrator/strategies/`


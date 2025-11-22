# PHASE 3: REACTION SYSTEM REFACTORING - SSOT-Compliant Design

**Status:** ⏳ In Progress  
**Branch:** `feature/reaction-system`  
**Start Date:** 2025-11-22  
**Estimated Duration:** 3 days

---

## 🎯 OBJECTIVES

### **Primary Goals:**
1. ✅ **SSOT Compliance:** Eliminate parameter duplication (no `ReactionContext` struct with copied T/n)
2. ✅ **Direct Database Access:** Use `config::ReactionDatabase` & `config::SpeciesDatabase` directly
3. ✅ **Modern Types:** Replace legacy `ReactionEntry` with `config::Reaction`
4. ✅ **Modular Design:** Extract reaction logic into `IReactionHandler` interface
5. ✅ **Multi-Order Support:** Handle 1st/2nd/3rd-order reactions with proper concentration handling

### **Success Criteria:**
- ✅ All reaction logic extracted into `IReactionHandler` implementations
- ✅ SSOT: Handler reads from `config::ReactionDatabase` / `config::SpeciesDatabase` / `config::EnvironmentConfig` (no copies)
- ✅ No legacy `ReactionEntry` structs (use `config::Reaction` directly)
- ✅ Factory pattern for handler creation (from `config.physics.enable_reactions`)
- ✅ `integrate_one_step()` uses injected handler (no inline logic)
- ✅ All tests pass (100% green)
- ✅ Legacy code deleted (`handle_reaction()`, `load_reactions()`, conversion in `main.cpp`)

---

## 📊 CURRENT STATE ANALYSIS

### **Problem: Scattered Reaction Logic**

#### **Location 1: `src/core/integrator/integrator_helpers.cpp:304`**
```cpp
void handle_reaction(IonState& y, EhssRng& rng, double dt, const GlobalParams& gParams,
                     const std::unordered_map<std::string, Species>& speciesDB,
                     const std::vector<ReactionEntry>& reaction_list) {
    // ❌ 50+ lines of inline reaction logic
    // ❌ Uses legacy GlobalParams
    // ❌ Uses legacy ReactionEntry struct
    // ❌ Mixed physics + conversion logic
}
```

#### **Location 2: `src/core/physics/reactions/reactionUtils.cpp:113`**
```cpp
std::vector<ReactionEntry> load_reactions(const GlobalParams& gParams) {
    // ❌ Loads from GlobalParams (legacy)
    // ❌ Returns vector<ReactionEntry> (not modern ReactionDatabase)
    // ❌ Duplicate of config::ReactionLoader
}
```

#### **Location 3: `src/main/main.cpp:359`**
```cpp
// ❌ Manual conversion: ReactionDatabase → vector<ReactionEntry>
for (const auto& rxn : full_config.reaction_db.reactions) {
    ReactionEntry entry;
    entry.reactant = rxn.reactant;
    entry.product = rxn.product;
    // ... 20 lines of conversion ...
    reaction_list.push_back(entry);
}
```

### **SSOT VIOLATIONS:**
1. 🔴 **Two reaction loaders** (`reactionUtils.cpp` vs `config::ReactionLoader`)
2. 🔴 **Manual struct conversion** (modern → legacy)
3. 🔴 **Legacy types propagated** (`ReactionEntry` instead of `config::Reaction`)
4. 🔴 **Parameter duplication** (temperature/density copied instead of referenced)

---

## ✅ SOLUTION: SSOT-Compliant Reaction System

### **Design Principles:**

| **Principle** | **Implementation** |
|--------------|-------------------|
| **SSOT Config** | Read directly from `config::ReactionDatabase` & `config::SpeciesDatabase` |
| **No Conversion** | Use modern types (`config::Reaction`, not `ReactionEntry`) |
| **No Parameter Duplication** | Pass `EnvironmentConfig&` directly (contains T, n) |
| **Factory Pattern** | Create handler from `config.physics.enable_reactions` |
| **Interface-Based** | `IReactionHandler` with modular implementations |

---

## 🏗️ ARCHITECTURE

### **Component Hierarchy:**

```
IReactionHandler (interface)
    ├── StochasticReactionHandler    (discrete reaction events)
    └── NoReactionHandler             (null object pattern)

ReactionHandlerFactory
    └── create(config.physics)  // Direct access to enable_reactions

Integration Point:
    integrate_one_step(ion, config, domain, collision_handler, reaction_handler, ...)
        └── reaction_handler->handle_reaction(ion, dt, rng, reaction_db, species_db, env)
```

### **Key Differences from Legacy:**
- ✅ **No `ReactionEntry`** (use `config::Reaction` directly)
- ✅ **No `load_reactions()`** (use `config.reaction_db` directly)
- ✅ **No struct conversion** (modern types throughout)
- ✅ **No parameter duplication** (pass `EnvironmentConfig&` not copied T/n values)

---

## 🔬 REACTION ORDER HANDLING

### **Strategy: Explicit Concentration Handling**

**Principle:** All concentrations in **SI units [m⁻³]** (SSOT: No unit ambiguity!)

### **Rate Constant Units by Reaction Order:**

| **Total Order** | **Rate Constant Unit** | **k_eff Unit** | **Example** |
|-----------------|------------------------|----------------|-------------|
| **1 (pseudo-first-order)** | [s⁻¹] | [s⁻¹] | Unimolecular decay |
| **2 (bimolecular)** | [m³/s] | [s⁻¹] | Ion + neutral → products |
| **3 (termolecular)** | [m⁶/s] | [s⁻¹] | Ion + neutral + M → cluster |

### **Effective Rate Calculation:**

```cpp
double k_eff = reaction.rate_constant_m3s;

for (const auto& term : reaction.order_terms) {
    double conc_m3 = (term.concentration_m3 > 0.0)
        ? term.concentration_m3              // Explicit concentration
        : env.particle_density_m_3;          // Fallback to buffer gas
    
    k_eff *= std::pow(conc_m3, term.exponent);
}

return k_eff;  // [s⁻¹]
```

### **Behavior:**

| **JSON Configuration** | **Interpretation** | **k_eff Formula** |
|------------------------|-------------------|-------------------|
| No `order` terms | Pseudo-first-order | `k_eff = k₀` [s⁻¹] |
| `concentration_m3 = 2.5e25` | Fixed neutral density | `k_eff = k₀ × (2.5e25)^n` [s⁻¹] |
| `concentration_m3 = -1.0` (or missing) | Use buffer gas density | `k_eff = k₀ × (n_gas)^n` [s⁻¹] |
| `concentration_m3 = 0` | Zero concentration (no neutral) | `k_eff = 0` [s⁻¹] |
| `concentration_m3 < -1.0` or `(-1.0, 0)` | **ERROR** (invalid) | Validation error |

---

## 📋 IMPLEMENTATION PLAN

### **Phase 3A: Design & Interface (Day 1)**

#### **File Structure:**

```
src/core/physics/reactions/
├── IReactionHandler.h                 ✅ Interface (SSOT-compliant)
├── StochasticReactionHandler.h        ✅ Implementation
├── StochasticReactionHandler.cpp
├── NoReactionHandler.h                ✅ Null object pattern
├── ReactionHandlerFactory.h           ✅ Factory
└── ReactionHandlerFactory.cpp

tests/physics/reactions/
├── CMakeLists.txt
├── test_stochastic_reaction_handler.cpp
├── test_reaction_factory.cpp
└── test_reaction_order.cpp
```

#### **Tasks:**

1. ✅ Create `IReactionHandler.h` (interface with 5 parameters: ion, dt, rng, reaction_db, species_db, env)
2. ✅ Create `StochasticReactionHandler.h/.cpp` (reads from `env.temperature_K`, `env.particle_density_m_3`)
3. ✅ Implement `compute_effective_rate()` (multi-order support)
4. ✅ Implement `find_applicable_reactions()` (filter by ion species)
5. ✅ Implement `update_ion_species()` (lookup in species_db)
6. ✅ Write unit tests (2nd-order, 3rd-order, buffer gas fallback)

**Estimated Time:** 2-3 hours

---

### **Phase 3B: Factory & Integration (Day 2)**

#### **Tasks:**

7. ✅ Create `ReactionHandlerFactory.h/.cpp`
8. ✅ Update `integrate_one_step()` (add `reaction_handler` parameter)
9. ✅ Pass 3 parameters directly: `config.reaction_db`, `config.species_db`, `domain.environment`
10. ✅ Update `integrate_trajectory()` (create handler from factory)
11. ✅ Integration tests (end-to-end with JSON config)

**Estimated Time:** 2 hours

---

### **Phase 3C: Factory Integration (Day 2-3)** ✅ COMPLETED

#### **Status:** ✅ Handler factory created and instantiated in integrate_trajectory

#### **Completed Tasks:**

12. ✅ Added `ReactionHandlerFactory` includes to `integrator.cpp`
13. ✅ Moved `PhysicsConfig` creation outside collision handler scope
14. ✅ Created reaction handler from factory in `integrate_trajectory()`
15. ✅ Added logger output for handler creation
16. ✅ All 23 tests passing

#### **Commits:**
- `9791baa`: Factory integration in integrate_trajectory
- `d9a0871`: integrate_one_step signature with reaction_handler parameter
- `3a35748`: Competing channels fix (weighted selection)
- `8a35840`: Handler system (IReactionHandler, StochasticReactionHandler, Factory)

**Time Spent:** 3 hours

---

### **Phase 3D: Database Migration (Day 4)** ⏳ TODO

#### **BLOCKER IDENTIFIED:**

**Type Mismatch Prevents Full Integration:**
- `integrator.cpp` uses `ICARION::io::Species` (from speciesLoader.h)
- `integrator_helpers.cpp` uses `Species` (from reactionUtils.h)
- `reaction_handler` expects `config::SpeciesDatabase` (from SpeciesConfig.h)

**Current Data Flow (Suboptimal):**
```
JSON → config::SpeciesDatabase (SSOT)
     → Adapter in main.cpp (converts to legacy)
     → ICARION::io::SpeciesDatabase (integrator.cpp)
     → std::unordered_map<string, Species> (integrator_helpers.cpp)
     → ❌ Cannot convert back to config::SpeciesDatabase for reaction_handler
```

**Target Data Flow (After Phase 3D):**
```
JSON → config::SpeciesDatabase (SSOT)
     → integrate_trajectory() directly
     → reaction_handler (no conversion!)
```

#### **Migration Plan:**

1. ❌ **Unify Species Types:**
   - Remove `ICARION::io::Species` (speciesLoader.h)
   - Remove `Species` (reactionUtils.h)
   - Keep only `config::SpeciesProperties` (SpeciesConfig.h)

2. ❌ **Update integrate_trajectory Signature:**
   ```cpp
   // OLD:
   integrate_trajectory(..., const ICARION::io::SpeciesDatabase& speciesDB, ...)
   
   // NEW:
   integrate_trajectory(..., const config::SpeciesDatabase& species_db, ...)
   ```

3. ❌ **Wire reaction_handler to integrate_one_step:**
   ```cpp
   // Currently disabled due to type mismatch:
   if (reaction_handler != nullptr) {
       reaction_handler->handle_reaction(ion, dt, rng, reaction_db, species_db, env);
   }
   ```

4. ❌ **Delete Legacy Code:**
   - `src/core/physics/reactions/reactionUtils.cpp:113-142` (load_reactions)
   - `src/core/integrator/integrator_helpers.cpp:316` (handle_reaction)
   - `src/main/main.cpp:359-400` (ReactionEntry conversion)

5. ❌ **Update main.cpp:**
   - Remove legacy adapter (SSOT config → Legacy types)
   - Pass SSOT databases directly to integrator

#### **Dependencies:**

- ⚠️ **Force System** must also migrate to SSOT (currently uses GlobalParams)
- ⚠️ **All systems** using legacy GlobalParams must be modernized first

#### **Estimated Time:** 4-6 hours (complex cross-cutting change)

---

### **Phase 3E: Full SSOT Integration (Day 5)** ⏳ FUTURE

**Remove Legacy Adapter Entirely:**
- Delete main.cpp conversion code (SSOT → Legacy)
- Integrate force system with SSOT configs
- Remove GlobalParams (replace with modern config types)

---

## 🧪 TEST STRATEGY

### **Unit Tests:**

#### **Test 1: Second-Order Reaction**
```cpp
TEST_CASE("StochasticReactionHandler: Second-order reaction") {
    // Reaction: H3O+ + NH3 → NH4+
    // k₀ = 2.0e-9 [m³/s]
    // [NH3] = 1e20 [m⁻³]
    // k_eff = 2e-9 × 1e20 = 2e11 [s⁻¹]
    // P ≈ 1.0 (certain)
}
```

#### **Test 2: Third-Order Reaction**
```cpp
TEST_CASE("StochasticReactionHandler: Third-order reaction") {
    // Reaction: H3O+ + H2O + He → H5O2+ + He
    // k₀ = 1.2e-28 [m⁶/s]
    // [H2O] = 2.5e25 [m⁻³], [He] = 2.5e25 [m⁻³]
    // k_eff = 1.2e-28 × 2.5e25 × 2.5e25 = 7.5e22 [s⁻¹]
}
```

#### **Test 3: Buffer Gas Fallback**
```cpp
TEST_CASE("StochasticReactionHandler: Buffer gas fallback") {
    // Order term with concentration_m3 = 0
    // → Should use env.particle_density_m_3
}
```

#### **Test 4: SSOT Compliance**
```cpp
TEST_CASE("StochasticReactionHandler: No parameter duplication") {
    // Verify temperature/density read from EnvironmentConfig
    // NOT from copied struct fields
}
```

### **Integration Tests:**

#### **Test 5: End-to-End with JSON**
```cpp
TEST_CASE("Reaction system integration") {
    // Load full config from JSON
    // Create handler from factory
    // Run 1000 timesteps
    // Verify species conversion statistics
}
```

---

## 📊 SSOT COMPLIANCE CHECKLIST

### **Phase 3: Reaction System**

- [ ] **IReactionHandler** takes direct references (`ReactionDatabase&`, `SpeciesDatabase&`, `EnvironmentConfig&`)
- [ ] **No `ReactionContext` struct** (SSOT violation avoided!)
- [ ] **No parameter duplication** (reads from `env.temperature_K`, `env.particle_density_m_3`)
- [ ] **No `ReactionEntry` struct** (uses `config::Reaction` directly)
- [ ] **No legacy `load_reactions()`** (uses `config.reaction_db`)
- [ ] **ReactionHandlerFactory** uses `config.physics.enable_reactions` (direct)
- [ ] **`integrate_one_step()`** passes databases directly (no struct creation)
- [ ] **All tests passing** (100% green)

---

## 📈 PROGRESS TRACKING

### **Day 1: Interface & Handler** ⏳
- [ ] Create `IReactionHandler.h`
- [ ] Create `StochasticReactionHandler.h/.cpp`
- [ ] Implement `compute_effective_rate()`
- [ ] Implement `find_applicable_reactions()`
- [ ] Implement `update_ion_species()`
- [ ] Write unit tests (3 test cases)

### **Day 2: Factory & Integration** ⏳
- [ ] Create `ReactionHandlerFactory.h/.cpp`
- [ ] Update `integrate_one_step()`
- [ ] Update `integrate_trajectory()`
- [ ] Integration tests

### **Day 3: Cleanup** ⏳
- [ ] Delete legacy `handle_reaction()`
- [ ] Delete legacy `load_reactions()`
- [ ] Delete conversion code in `main.cpp`
- [ ] Build & test
- [ ] Commit & merge

---

## 🎯 VALIDATION CRITERIA

### **Before Merge:**

1. ✅ **Code Quality:**
   - All handlers implement `IReactionHandler` interface
   - SSOT: No parameter duplication (checked manually)
   - Factory returns correct handler type
   - No legacy `ReactionEntry` usage (grep verified)

2. ✅ **Testing:**
   - All unit tests passing (100%)
   - Integration tests passing
   - CTest suite passing (21+ tests)
   - No memory leaks (valgrind clean)

3. ✅ **Documentation:**
   - Handler interfaces documented (Doxygen)
   - Reaction order handling explained
   - Migration guide for users (if needed)

4. ✅ **Cleanup:**
   - Legacy code deleted (3 locations)
   - No dead code remaining
   - Build clean (no warnings)

---

## 📝 NOTES

### **Key Design Decisions:**

1. **No `ReactionContext` struct** → Avoids parameter duplication (T, n already in `EnvironmentConfig`)
2. **Multi-order support** → Handles 1st/2nd/3rd-order reactions with explicit concentration handling
3. **Buffer gas fallback** → `concentration_m3 = 0` → use `env.particle_density_m_3`
4. **Direct database access** → No intermediate copies (SSOT!)
5. **Same pattern as Collision System** → Consistent architecture across codebase

### **Lessons from Phase 2 (Collision System):**

- ✅ Interface-based design enables testing and modularity
- ✅ Factory pattern simplifies handler creation
- ✅ Direct config references avoid SSOT violations
- ✅ Extensive unit tests catch integration bugs early
- ✅ Delete legacy code immediately (don't leave it commented out)

---

## 🚀 NEXT STEPS AFTER PHASE 3

**Phase 4: Force Registry & Integration Strategy** (Future work)
- Extract force computation logic
- Create `ForceRegistry` with `IForceProvider` interface
- Implement `IntegrationStrategy` (RK4, Verlet, etc.)
- Complete integrator refactoring (40% → 70%)

---

## 🚀 FUTURE ENHANCEMENTS (Phase 3 - Later)

### **Phase 3B+: Competing Reaction Channels** 🔴 **CRITICAL**

**Problem:** Current implementation tests reactions sequentially → First reaction biased!

**Physics:** When ion has multiple reactions (k₁, k₂, ..., kₙ):
- **Total rate:** k_total = Σ kᵢ
- **Total probability:** P_total = 1 - exp(-k_total × dt)
- **Channel selection:** P(channel i) = kᵢ / k_total

**Algorithm:**
```cpp
// Step 1: Compute k_total
double k_total = 0.0;
for (auto& rxn : reactions) {
    k_effs.push_back(compute_k_eff(rxn));
    k_total += k_effs.back();
}

// Step 2: Total reaction probability
double P = 1 - exp(-k_total * dt);
if (rng() >= P) return false;

// Step 3: Weighted channel selection
double r = rng() * k_total;
double cumulative = 0.0;
for (size_t i = 0; i < reactions.size(); ++i) {
    cumulative += k_effs[i];
    if (r < cumulative) {
        // React via channel i
        return true;
    }
}
```

**Implementation:** Phase 3B (before integrator integration)

---

### **Phase 3D: Temperature-Dependent Rate Constants**

**Motivation:** Ion-molecule reactions often have T-dependence:
- **Arrhenius:** k(T) = A × exp(-Eₐ / (kB T)) (activated reactions)
- **Capture:** k(T) = C × (T/300)ⁿ (ion-dipole capture)

**Config Extension:**
```cpp
struct Reaction {
    enum class RateType { Constant, Arrhenius, Capture } rate_type;
    double rate_constant_m3s;  // k₀ (or A for Arrhenius)
    
    // Arrhenius parameters
    double activation_energy_J;  // Eₐ
    
    // Capture parameters
    double temperature_exponent;  // n
};
```

**JSON Example:**
```json
{
  "id": "rxn_arrhenius",
  "reactant": "H3O+",
  "product": "NH4+",
  "rate_type": "arrhenius",
  "rate_constant_m3s": 2.0e-9,
  "activation_energy_J": 1.5e-20
}
```

**Implementation:** Phase 3D (after basic integration works)

---

### **Phase 3E: Multi-Gas Mixtures**

**Motivation:** Realistic environments have gas mixtures (e.g., He + H₂O, N₂ + VOCs)

**Config Extension:**
```cpp
struct EnvironmentConfig {
    struct GasComponent {
        std::string species_id;
        double mole_fraction;
        double density_m_3;
    };
    
    std::vector<GasComponent> gas_mixture;
    
    double get_density(const std::string& species_id) const;
};
```

**JSON Example:**
```json
{
  "environment": {
    "temperature_K": 300,
    "pressure_Pa": 101325,
    "gas_mixture": [
      {"species": "He", "mole_fraction": 0.95},
      {"species": "H2O", "mole_fraction": 0.05}
    ]
  },
  "reactions": [
    {
      "reactant": "H3O+",
      "product": "H5O2+",
      "rate_constant_m3s": 1.2e-28,
      "order": [
        {"species": "H2O", "exponent": 1, "concentration_m3": -1.0},
        {"species": "He", "exponent": 1, "concentration_m3": -1.0}
      ]
    }
  ]
}
```

**Handler Adaptation:**
```cpp
double conc_m3 = (term.concentration_m3 < 0.0)
    ? env.get_density(term.species)  // ✅ Lookup in gas_mixture!
    : term.concentration_m3;
```

**Implementation:** Phase 3E (requires EnvironmentConfig redesign)

---

## 📊 PRIORITY RANKING

| **Priority** | **Feature** | **Effort** | **Impact** | **Phase** |
|-------------|-----------|----------|----------|---------|
| 🔴 **P0 CRITICAL** | Competing channels | 1 hour | **HIGH** (correctness!) | 3B |
| 🟡 **P1 High** | Temperature dependence | 2 hours | Medium (common reactions) | 3D |
| 🟢 **P2 Medium** | Multi-gas mixtures | 3 hours | Medium (realistic sims) | 3E |

---

**Document Created:** 2025-11-22  
**Last Updated:** 2025-11-22  
**Author:** GitHub Copilot + User (chsch95)

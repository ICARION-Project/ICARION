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

### **Phase 3C: Legacy Cleanup (Day 3)**

#### **Files to DELETE:**

1. ❌ `src/core/physics/reactions/reactionUtils.cpp:113-142` (`load_reactions()` function)
2. ❌ `src/core/integrator/integrator_helpers.cpp:304` (`handle_reaction()` function)
3. ❌ `src/main/main.cpp:359-400` (manual `ReactionEntry` conversion)

#### **Files to UPDATE:**

```diff
// src/core/physics/reactions/reactionUtils.h

-/**
- * @brief Load ion-molecule reactions from JSON file
- * 
- * @deprecated Use config::ReactionLoader instead
- */
-std::vector<ReactionEntry> load_reactions(const GlobalParams& gParams);

+// DEPRECATED: Use config::ReactionLoader and ReactionHandlerFactory instead
+// This function removed in Phase 3 refactor (2025-11-22)
```

#### **Tasks:**

12. ✅ Delete legacy `handle_reaction()`
13. ✅ Delete legacy `load_reactions()`
14. ✅ Delete conversion code in `main.cpp`
15. ✅ Update deprecation comments
16. ✅ Build & test (all 21+ tests passing)

**Estimated Time:** 1 hour

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

**Document Created:** 2025-11-22  
**Last Updated:** 2025-11-22  
**Author:** GitHub Copilot + User (chsch95)

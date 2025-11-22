# Phase 3 Completion Plan - Option C (Hybrid Approach)

**Strategy:** Complete reaction system in current branch, defer database migration to separate branches

**Date:** 2025-11-22  
**Branch Strategy:** 3 sequential feature branches  
**Estimated Total Time:** 12-16 hours (split over 3 branches)

---

## 🎯 BRANCH 1: feature/reaction-system (CURRENT) ✅

**Goal:** Complete reaction system refactoring as far as possible without touching Force System or Database types

**Status:** ⏳ In Progress (90% complete)  
**Time Remaining:** 1-2 hours  
**Merge Target:** `main` or `integrator-refactor`

### Tasks Remaining:

#### 1. Mark Legacy Code as Deprecated (30 min)
- [ ] `src/core/integrator/integrator_helpers.cpp:316` - Add `@deprecated` comment to `handle_reaction()`
- [ ] `src/core/physics/reactions/reactionUtils.cpp:113` - Add `@deprecated` comment to `load_reactions()`
- [ ] `src/core/physics/reactions/reactionUtils.cpp:56` - Add `@deprecated` comment to `load_speciesDB()`
- [ ] `src/core/physics/reactions/reactionUtils.h` - Add deprecation notices in header

#### 2. Comment Adapter Code in main.cpp (10 min)
- [ ] Add `// TODO Phase 3D: Remove after database migration` to ReactionEntry conversion
- [ ] Add `// TODO Phase 3D: Remove after database migration` to Species conversion

#### 3. Add Edge-Case Tests (30 min)
- [ ] Test: Zero reactions available (handler should return false)
- [ ] Test: Very large k_eff (numerical stability check)
- [ ] Test: Very small k_eff (rare event probability)

#### 4. Update Documentation (20 min)
- [ ] `docs/ARCHITECTURE.md` - Add Reaction System section
- [ ] `README.md` - Update Phase 3 status
- [ ] `REACTION_SYSTEM_REFACTORING_PLAN.md` - Mark Phase 3C complete

#### 5. Final Commit & Merge (10 min)
- [ ] Commit: "feat(reactions): Complete Phase 3C - Modern handler system ready"
- [ ] Push to remote
- [ ] Create PR: "Phase 3: Modern Reaction Handler System"
- [ ] Merge to main (or integrator-refactor)
- [ ] Tag: `v1.1-reactions-modern-handler`

**Deliverables:**
- ✅ Modern reaction handler system (IReactionHandler, StochasticReactionHandler, Factory)
- ✅ Competing channels algorithm (weighted selection)
- ✅ 10+ unit tests (all passing)
- ✅ Legacy code marked as deprecated (but still functional)
- ✅ Documentation complete
- ✅ All 23 tests passing

---

## 🎯 BRANCH 2: feature/database-unification (FUTURE)

**Goal:** Unify Species types and wire reaction_handler to integrate_one_step

**Status:** 📋 Planned (not started)  
**Estimated Time:** 4-6 hours  
**Dependencies:** Branch 1 must be merged first

### Scope:

#### 1. Unify Species Types (2h)
- Remove `ICARION::io::Species` (src/core/io/speciesLoader.h)
- Remove `Species` (src/core/physics/reactions/reactionUtils.h)
- Keep only `config::SpeciesProperties` (src/core/config/types/SpeciesConfig.h)
- Update all references

#### 2. Update integrate_trajectory Signature (1h)
```cpp
// OLD:
integrate_trajectory(..., const ICARION::io::SpeciesDatabase& speciesDB, ...)

// NEW:
integrate_trajectory(..., const config::SpeciesDatabase& species_db, ...)
```

#### 3. Wire reaction_handler to integrate_one_step (1h)
```cpp
// Enable this code in integrator_helpers.cpp:
if (reaction_handler != nullptr) {
    reaction_handler->handle_reaction(ion, dt, rng, reaction_db, species_db, env);
} else {
    handle_reaction(ion, local_rng, dt_local, gParams, speciesDB, reaction_list);
}
```

#### 4. Update main.cpp (1h)
- Remove ReactionEntry conversion (lines ~359-400)
- Pass SSOT config::SpeciesDatabase directly
- Remove legacy adapter calls

#### 5. Delete Legacy Code (30min)
- Delete `handle_reaction()` from integrator_helpers.cpp
- Delete `load_reactions()` from reactionUtils.cpp
- Delete `load_speciesDB()` from reactionUtils.cpp

**Deliverables:**
- ✅ Single Species type (config::SpeciesProperties)
- ✅ reaction_handler fully wired and active
- ✅ Legacy reaction code deleted
- ✅ All tests passing

**Blockers:**
- None (independent of Force System)

---

## 🎯 BRANCH 3: feature/force-system-ssot (FUTURE)

**Goal:** Migrate Force System to SSOT config types, eliminate GlobalParams

**Status:** 📋 Planned (not started)  
**Estimated Time:** 6-8 hours  
**Dependencies:** Branch 2 must be merged first

### Scope:

#### 1. Create ForceConfig Types (2h)
- `config::ForceConfig` (electric field, magnetic field, damping)
- `config::ElectricFieldConfig` (DC, RF, custom fields)
- `config::MagneticFieldConfig` (B field parameters)

#### 2. Update compute_accelerations() (2h)
```cpp
// OLD:
IonState compute_accelerations(double t, const IonState& y, 
    const GlobalParams& gParams, const InstrumentDomain& dom, ...);

// NEW:
IonState compute_accelerations(double t, const IonState& y,
    const config::ForceConfig& force_config, const InstrumentDomain& dom, ...);
```

#### 3. Eliminate GlobalParams (3h)
- Replace all GlobalParams references with config types
- Update integrate_trajectory() signature
- Update main.cpp to use modern configs

#### 4. Delete LegacyAdapter (1h)
- Remove `src/core/config/adapter/LegacyAdapter.cpp`
- Remove conversion code from main.cpp
- Full SSOT compliance achieved

**Deliverables:**
- ✅ GlobalParams completely eliminated
- ✅ All systems use SSOT config types
- ✅ LegacyAdapter deleted
- ✅ Full SSOT compliance

---

## 📊 Progress Tracking

| Branch | Status | Commits | Tests | Time |
|--------|--------|---------|-------|------|
| **feature/reaction-system** | ⏳ 90% | 5/6 | 23/23 ✅ | 1-2h |
| **feature/database-unification** | 📋 Planned | 0/? | - | 4-6h |
| **feature/force-system-ssot** | 📋 Planned | 0/? | - | 6-8h |

**Total Progress:** 30% complete (Branch 1 nearly done)

---

## ✅ Success Criteria

### Branch 1 (Reaction System):
- [x] Modern handler system exists
- [x] Factory creates handlers
- [x] Competing channels algorithm correct
- [x] 10+ unit tests passing
- [ ] Legacy code marked deprecated
- [ ] Documentation complete

### Branch 2 (Database Unification):
- [ ] Single Species type
- [ ] reaction_handler wired and active
- [ ] Legacy reaction code deleted
- [ ] All tests passing

### Branch 3 (Force System SSOT):
- [ ] GlobalParams eliminated
- [ ] LegacyAdapter deleted
- [ ] Full SSOT compliance
- [ ] All tests passing

---

## 🚀 Immediate Next Steps (Branch 1)

1. **Mark legacy functions as deprecated** (30 min)
2. **Add TODO comments to main.cpp** (10 min)
3. **Add 3 edge-case tests** (30 min)
4. **Update ARCHITECTURE.md** (20 min)
5. **Commit and merge** (10 min)

**Total:** ~1.5 hours to complete Branch 1

---

## 📝 Notes

**Why 3 branches?**
- Clear separation of concerns
- Smaller, reviewable PRs
- Tests stay green between merges
- Easier to revert if issues found
- Each branch has clear deliverables

**Why not do it all now?**
- Force System is complex (6-8h)
- Avoid scope creep
- Keep feature/reaction-system focused
- Better git history

**Current Blocker:**
Type mismatch prevents wiring reaction_handler:
- `ICARION::io::Species` (integrator.cpp)
- `Species` (reactionUtils.h)
- `config::SpeciesProperties` (SpeciesConfig.h)

**Resolution:** Branch 2 will unify these types

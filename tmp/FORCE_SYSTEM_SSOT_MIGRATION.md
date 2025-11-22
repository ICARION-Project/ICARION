# FORCE SYSTEM SSOT MIGRATION PLAN

**Branch:** `refactor/force-system-ssot`  
**Base:** `integrator-refactor` (includes reaction system)  
**Goal:** Remove parameter structs, use FullConfig/DomainConfig directly  

---

## 🎯 OBJECTIVE

**Eliminate SSOT violations in Force System:**

### **Current State (SSOT Violations):**
- ❌ `MagneticFieldParams` duplicates `DomainConfig.fields.magnetic`
- ❌ `AnalyticalFieldParams` duplicates `DomainConfig.fields.*`
- ❌ `DampingParams` duplicates `DomainConfig.environment.*`
- ❌ `compute_accelerations()` uses `GlobalParams` + `InstrumentDomain` (legacy)

### **Target State (SSOT Compliant):**
- ✅ Forces take `const DomainConfig&` reference
- ✅ Read directly from `domain.fields.dc.axial_V`, etc.
- ✅ No parameter structs (deleted)
- ✅ `compute_accelerations()` uses `DomainConfig` (modern)

---

## 📋 IMPLEMENTATION STEPS

### **Step 1: Update MagneticFieldForce (30 min)**

**Files:**
- `src/core/physics/forces/MagneticFieldForce.h`
- `src/core/physics/forces/MagneticFieldForce.cpp`

**Changes:**
1. Delete `MagneticFieldParams` struct
2. Add constructor: `MagneticFieldForce(const config::MagneticFieldConfig& magnetic)`
3. Store reference: `const config::MagneticFieldConfig& magnetic_;`
4. Update `compute()` to read from `magnetic_.field_strength_T`

---

### **Step 2: Update ElectricFieldForce (45 min)**

**Files:**
- `src/core/physics/forces/ElectricFieldForce.h`
- `src/core/physics/forces/ElectricFieldForce.cpp`

**Changes:**
1. Delete `AnalyticalFieldParams` struct
2. Add constructor: `ElectricFieldForce(config::Instrument, const config::GeometryConfig&, const config::FieldsConfig&)`
3. Store references to `instrument_`, `geometry_`, `fields_`
4. Update all instrument-specific methods to read from `fields_.dc.*`, `fields_.rf.*`

---

### **Step 3: Update DampingForce (30 min)**

**Files:**
- `src/core/physics/forces/DampingForce.h`
- `src/core/physics/forces/DampingForce.cpp`

**Changes:**
1. Delete `DampingParams` struct
2. Add constructor: `DampingForce(const config::EnvironmentConfig&, DampingModel)`
3. Store reference: `const config::EnvironmentConfig& env_;`
4. Update `compute()` to read from `env_.pressure_Pa`, `env_.temperature_K`

---

### **Step 4: Update compute_accelerations() signature (45 min)**

**Files:**
- `src/core/physics/computeAccelerations.h`
- `src/core/physics/computeAccelerations.cpp`

**Changes:**
1. Change signature:
   ```cpp
   // OLD:
   IonState compute_accelerations(..., const GlobalParams&, const InstrumentDomain&, ...);
   
   // NEW:
   IonState compute_accelerations(..., const config::DomainConfig&, ...);
   ```

2. Update all field/force computations to use `domain.fields.*`, `domain.environment.*`

---

### **Step 5: Update integrator_helpers.cpp (60 min)**

**Files:**
- `src/core/integrator/integrator_helpers.h`
- `src/core/integrator/integrator_helpers.cpp`

**Changes:**
1. Update `integrate_one_step()`:
   - Change signature: `const config::DomainConfig& domain` (not `InstrumentDomain`)
   - Update `compute_accelerations()` call

2. Update `integrate_trajectory()`:
   - Change signature: `const config::FullConfig& config` (not `GlobalParams + vector<InstrumentDomain>`)
   - Update all calls

---

### **Step 6: Update integrator.cpp (30 min)**

**Files:**
- `src/core/integrator/integrator.cpp`

**Changes:**
1. Update all `compute_accelerations()` calls (RK4, DOPRI5)
2. Use `DomainConfig` instead of `InstrumentDomain`

---

### **Step 7: Update main.cpp (30 min)**

**Files:**
- `src/main/main.cpp`

**Changes:**
1. Delete `LegacyAdapter::to_global_params()` call
2. Delete `LegacyAdapter::to_instrument_domains()` call
3. Pass `full_config` directly to `integrate_trajectory()`

---

### **Step 8: Update all tests (60 min)**

**Files:**
- `tests/physics/forces/test_magnetic_field_force.cpp`
- `tests/physics/test_compute_accelerations.cpp`
- `tests/integrator/test_integrator_helpers.cpp`

**Changes:**
1. Replace parameter structs with config references
2. Update all constructor calls
3. Verify tests still pass

---

### **Step 9: Delete legacy code (15 min)**

**Files to DELETE:**
- `src/core/config/adapter/LegacyAdapter.h`
- `src/core/config/adapter/LegacyAdapter.cpp`
- `src/core/param/paramUtils.h` (if GlobalParams/InstrumentDomain deleted)
- `src/core/param/paramUtils.cpp`

**Verify:**
- No references to `GlobalParams` remain
- No references to `InstrumentDomain` remain
- No references to `MagneticFieldParams`, `AnalyticalFieldParams`, `DampingParams`

---

### **Step 10: Update documentation (30 min)**

**Files:**
- `docs/ARCHITECTURE.md`
- `docs/DEVELOPERS_GUIDE.md`

**Changes:**
1. Update "Force System" section (SSOT-compliant architecture)
2. Add migration notes
3. Update code examples

---

## ⏱️ TIME ESTIMATE

| Step | Task | Time |
|------|------|------|
| 1 | MagneticFieldForce | 30 min |
| 2 | ElectricFieldForce | 45 min |
| 3 | DampingForce | 30 min |
| 4 | compute_accelerations() | 45 min |
| 5 | integrator_helpers | 60 min |
| 6 | integrator.cpp | 30 min |
| 7 | main.cpp | 30 min |
| 8 | Tests | 60 min |
| 9 | Delete legacy | 15 min |
| 10 | Documentation | 30 min |
| **TOTAL** | | **6h 15min** |

---

## ✅ SUCCESS CRITERIA

- [ ] No `MagneticFieldParams`, `AnalyticalFieldParams`, `DampingParams` structs
- [ ] No `GlobalParams` in force/integrator code
- [ ] No `InstrumentDomain` in force/integrator code
- [ ] Forces take `const DomainConfig&` or sub-config references
- [ ] All parameters read from config (no copies)
- [ ] All tests passing (100%)
- [ ] No compiler warnings
- [ ] Documentation updated

---

## 🚀 EXECUTION

**Start:** 2025-11-22  
**Branch:** `refactor/force-system-ssot`  
**Commit Strategy:** One commit per step (10 commits total)


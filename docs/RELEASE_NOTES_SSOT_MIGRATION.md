# ICARION v1.1 - Force System SSOT Migration

**Release Date:** TBD  
**Breaking Changes:** YES  
**Migration Required:** YES (if using force constructors directly)

---

## 🎯 Overview

ICARION v1.1 migrates the Force System to follow the **Single Source of Truth (SSOT)** principle. Forces now read directly from config references instead of copying data into parameter structs.

**Key Changes:**
- ✅ Parameter structs deleted (~200 lines of duplication removed)
- ✅ Forces use `const config&` references
- ✅ No data duplication
- ✅ Cleaner interfaces

---

## 🚨 Breaking Changes

### 1. Force Constructors Changed

#### MagneticFieldForce

```cpp
// ❌ OLD (v1.0 - REMOVED):
#include "core/physics/forces/MagneticFieldParams.h"

MagneticFieldParams params;
params.uniform_field_T = Vec3{0, 0, 1.5};
params.gradient_T_m = Vec3{0, 0, 0.1};
params.enabled = true;

MagneticFieldForce force(params);

// ✅ NEW (v1.1):
#include "core/config/types/DomainConfig.h"

const config::MagneticFieldConfig& magnetic = domain.fields.magnetic;
MagneticFieldForce force(magnetic);
```

#### ElectricFieldForce

```cpp
// ❌ OLD (v1.0 - REMOVED):
#include "core/physics/forces/AnalyticalFieldParams.h"

AnalyticalFieldParams params;
params.instrument_type = InstrumentType::IMS;
params.dc_axial_voltage_V = 1000.0;
params.length_m = 0.1;

ElectricFieldForce force(params);

// ✅ NEW (v1.1):
#include "core/config/types/DomainConfig.h"

const config::DomainConfig& domain = ...; // From config loader
ElectricFieldForce force(domain);
```

#### DampingForce

```cpp
// ❌ OLD (v1.0 - REMOVED):
#include "core/physics/forces/DampingParams.h"

DampingParams params;
params.pressure_Pa = 101325.0;
params.temperature_K = 300.0;
params.gas_mass_kg = 28.0 * AMU_TO_KG;

DampingForce force(params, DampingModel::HardSphere);

// ✅ NEW (v1.1):
#include "core/config/types/EnvironmentConfig.h"

const config::EnvironmentConfig& env = domain.environment;
DampingForce force(env, DampingModel::HardSphere);
```

### 2. Deleted Files

The following parameter struct files have been **permanently deleted**:

```
src/core/physics/forces/MagneticFieldParams.h     (DELETED)
src/core/physics/forces/AnalyticalFieldParams.h   (DELETED)
src/core/physics/forces/DampingParams.h           (DELETED)
```

If your code includes these headers, **compilation will fail**.

---

## 📦 Migration Guide

### Step 1: Update Includes

Replace parameter struct includes with config includes:

```cpp
// ❌ Remove these:
#include "core/physics/forces/MagneticFieldParams.h"
#include "core/physics/forces/AnalyticalFieldParams.h"
#include "core/physics/forces/DampingParams.h"

// ✅ Add these:
#include "core/config/types/DomainConfig.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/config/types/FieldsConfig.h"
```

### Step 2: Load Config (if not already done)

```cpp
#include "core/config/loader/ConfigLoader.h"

// Load full config from JSON
config::FullConfig full_config = config::load_config("config.json");

// Access domain config
const auto& domain = full_config.domains[0];
```

### Step 3: Update Force Construction

```cpp
// Create force registry
ForceRegistry registry;

// ✅ Pass config references (not parameter structs!)
registry.add_force(std::make_unique<ElectricFieldForce>(domain));
registry.add_force(std::make_unique<MagneticFieldForce>(domain.fields.magnetic));
registry.add_force(std::make_unique<DampingForce>(domain.environment, DampingModel::Friction));
```

### Step 4: Ensure Config Lifetime

⚠️ **CRITICAL**: Config objects **must outlive** force objects!

```cpp
// ✅ CORRECT: config outlives forces
{
    config::DomainConfig domain = load_config(...);
    
    ForceRegistry registry;
    registry.add_force(std::make_unique<ElectricFieldForce>(domain));
    
    // ... use registry ...
    
} // domain destructs AFTER registry

// ❌ WRONG: config destructs before forces!
{
    ForceRegistry registry;
    
    {
        config::DomainConfig domain = load_config(...);
        registry.add_force(std::make_unique<ElectricFieldForce>(domain));
    } // domain destructs HERE - forces now have dangling references!
    
    // ... use registry ... ⚠️ UNDEFINED BEHAVIOR!
}
```

### Step 5: Update Tests

Replace parameter struct setup with config setup:

```cpp
// ❌ OLD (v1.0):
TEST_CASE("MyForceTest") {
    MagneticFieldParams params;
    params.uniform_field_T = Vec3{0, 0, 1.0};
    MagneticFieldForce force(params);
    // ...
}

// ✅ NEW (v1.1):
TEST_CASE("MyForceTest") {
    config::DomainConfig domain;
    domain.fields.magnetic.field_strength_T = Vec3{0, 0, 1.0};
    domain.fields.magnetic.enabled = true;
    
    MagneticFieldForce force(domain.fields.magnetic);
    // ...
}
```

---

## ✅ Benefits of SSOT Migration

### Before (v1.0): Parameter Duplication

```cpp
// Config (SSOT)
domain.fields.dc.axial_V = 1000.0;
domain.geometry.length_m = 0.1;

// Force parameters (DUPLICATE!)
AnalyticalFieldParams params;
params.dc_axial_voltage_V = 1000.0;  // ❌ Copied from config!
params.length_m = 0.1;               // ❌ Copied from config!

// Problem: Two sources of truth, must keep in sync!
```

### After (v1.1): Single Source of Truth

```cpp
// Config (SSOT)
domain.fields.dc.axial_V = 1000.0;
domain.geometry.length_m = 0.1;

// Force (reference)
ElectricFieldForce force(domain);  // ✅ References config directly!

// Benefit: One source of truth, always in sync!
```

### Advantages

1. **No duplication**: ~200 lines of parameter code deleted
2. **Type safety**: Use strongly-typed config structs
3. **Automatic propagation**: Config changes visible immediately
4. **Cleaner code**: Fewer includes, simpler constructors
5. **Better maintainability**: Update config in one place

---

## 🧪 Testing

All tests have been updated and pass (100%):

```bash
cd build
ctest -R "Force"

# Results:
# ForceRegistry .................. PASSED
# ElectricFieldForce ............. PASSED
# MagneticDampingForces .......... PASSED
# SpaceChargeForce ............... PASSED
# ForceIntegration ............... PASSED
```

---

## 📚 Updated Documentation

- **ARCHITECTURE.md**: Force System section updated with SSOT pattern
- **DEVELOPERS_GUIDE.md**: "Adding New Force Types" rewritten for SSOT
- **This document**: Migration guide and breaking changes

---

## 🔍 Example: Complete Migration

**Before (v1.0):**

```cpp
#include "core/physics/forces/MagneticFieldParams.h"
#include "core/physics/forces/AnalyticalFieldParams.h"
#include "core/physics/forces/DampingParams.h"

// Manual parameter setup
MagneticFieldParams mag_params;
mag_params.uniform_field_T = Vec3{0, 0, 1.5};
mag_params.enabled = true;

AnalyticalFieldParams field_params;
field_params.instrument_type = InstrumentType::IMS;
field_params.dc_axial_voltage_V = 1000.0;
field_params.length_m = 0.1;

DampingParams damp_params;
damp_params.pressure_Pa = 101325.0;
damp_params.temperature_K = 300.0;
damp_params.gas_mass_kg = 28.0 * AMU_TO_KG;

// Create forces
ForceRegistry registry;
registry.add_force(std::make_unique<MagneticFieldForce>(mag_params));
registry.add_force(std::make_unique<ElectricFieldForce>(field_params));
registry.add_force(std::make_unique<DampingForce>(damp_params, DampingModel::HardSphere));
```

**After (v1.1):**

```cpp
#include "core/config/loader/ConfigLoader.h"

// Load config from JSON (SSOT)
config::FullConfig full_config = config::load_config("config.json");
const auto& domain = full_config.domains[0];

// Create forces (reference config directly)
ForceRegistry registry;
registry.add_force(std::make_unique<ElectricFieldForce>(domain));
registry.add_force(std::make_unique<MagneticFieldForce>(domain.fields.magnetic));
registry.add_force(std::make_unique<DampingForce>(domain.environment, DampingModel::HardSphere));
```

**Result:**
- ✅ Fewer includes
- ✅ No manual parameter setup
- ✅ Config-driven (no hardcoded values)
- ✅ SSOT compliant

---

## ⚠️ Common Pitfalls

### Pitfall 1: Dangling References

```cpp
// ❌ WRONG: Config destructs before forces
ForceRegistry create_forces() {
    config::DomainConfig domain;
    domain.fields.magnetic.field_strength_T = Vec3{0, 0, 1.0};
    
    ForceRegistry registry;
    registry.add_force(std::make_unique<MagneticFieldForce>(domain.fields.magnetic));
    
    return registry;  // ⚠️ domain destructs, forces have dangling references!
}

// ✅ CORRECT: Config outlives forces
class Simulation {
    config::FullConfig config_;
    ForceRegistry registry_;
    
    void setup() {
        config_ = load_config("config.json");
        const auto& domain = config_.domains[0];
        
        registry_.add_force(std::make_unique<ElectricFieldForce>(domain));
        // ✅ config_ outlives registry_
    }
};
```

### Pitfall 2: Modifying Config After Force Creation

```cpp
config::DomainConfig domain;
domain.fields.magnetic.field_strength_T = Vec3{0, 0, 1.0};

MagneticFieldForce force(domain.fields.magnetic);

// ✅ This works - force sees new value immediately (SSOT!)
domain.fields.magnetic.field_strength_T = Vec3{0, 0, 2.0};

Vec3 F = force.compute(ion, t, ctx);  // Uses B = 2.0 T (updated value)
```

**Note:** This is a **feature**, not a bug! Forces always see current config values.

### Pitfall 3: Missing Includes

```cpp
// ❌ WRONG: Forgot to include config headers
#include "core/physics/forces/ElectricFieldForce.h"

ElectricFieldForce force(domain);  // ❌ Compiler error: 'domain' undefined

// ✅ CORRECT: Include config types
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/config/types/DomainConfig.h"

config::DomainConfig domain;
ElectricFieldForce force(domain);  // ✅ Compiles
```

---

## 🚀 Next Steps

After migrating to SSOT force system:

1. **Steps 5-8** (in progress): Migrate integrator and main.cpp
2. **Step 9**: Delete `LegacyAdapter` and `paramUtils`
3. **v1.2**: Migrate collision handlers to SSOT
4. **v1.3**: Migrate reaction system to SSOT

---

## 📞 Support

If you encounter issues during migration:

1. Check this migration guide
2. Review updated docs (`ARCHITECTURE.md`, `DEVELOPERS_GUIDE.md`)
3. Check unit tests for examples (`tests/physics/forces/`)
4. Open an issue on GitHub

---

**Last Updated:** 2025-11-22  
**Applies to:** ICARION v1.1+

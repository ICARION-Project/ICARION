# ICARION Developer's Guide

**Version:** 1.0++ (Post Phase 1: Force System Refactoring)  
**Last Updated:** November 21, 2025

This guide provides practical instructions for extending ICARION with new features.

---

## Table of Contents

1. [Adding New Force Types](#adding-new-force-types)
2. [Adding New Collision Models](#adding-new-collision-models)
3. [Adding New Instrument Types](#adding-new-instrument-types)
4. [Testing Guidelines](#testing-guidelines)
5. [Code Style and Conventions](#code-style-and-conventions)

---

## Adding New Force Types

### Overview

ICARION's force system follows a plugin architecture using the **IForce interface**. All forces implement `IForce::compute()` and are managed by `ForceRegistry`.

### Design Principle: Single Source of Truth (SSOT)

**✅ CORRECT PATTERN**: Forces receive **const references** to config objects, NOT parameter structs.

**❌ WRONG PATTERN**: Creating separate `YourForceParams` structs duplicates configuration data.

**Why?** Parameter duplication violates SSOT and creates maintenance burden. Forces should read directly from `config::DomainConfig` or similar config structs.

**Example**:
```cpp
// ✅ CORRECT: Direct config reference
class ElectricFieldForce : public IForce {
    ElectricFieldForce(const config::DomainConfig& domain);
private:
    const config::DomainConfig& domain_;  // Reference, not copy!
};

// ❌ WRONG: Parameter struct (duplicates config!)
struct ElectricFieldParams { /* ... */ };
class ElectricFieldForce : public IForce {
    ElectricFieldForce(const ElectricFieldParams& params);
private:
    ElectricFieldParams params_;  // Duplication!
};
```

---

### Step-by-Step Guide

#### 1. Create Force Header (`src/core/physics/forces/YourForce.h`)

**✅ CORRECT PATTERN**: Use **const reference** to config, not parameter struct!

```cpp
#pragma once

#include "IForce.h"
#include "ForceContext.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "core/config/DomainConfig.h"  // ✅ Direct config reference

namespace ICARION {
namespace physics {

/**
 * @brief Your force description
 * 
 * Physics: F = ... (formula)
 * 
 * ✅ Uses const reference to DomainConfig (SSOT pattern)
 */
class YourForce : public IForce {
public:
    /**
     * @brief Constructor
     * 
     * @param domain Reference to simulation domain config (SSOT)
     * @param additional_param Any non-config parameters
     * 
     * ⚠️ Config reference must outlive this object!
     */
    explicit YourForce(const config::DomainConfig& domain, double additional_param = 0.0);
    
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override;
    
private:
    const config::DomainConfig& domain_;  ///< ✅ Reference, not copy!
    double additional_param_;
};

} // namespace physics
} // namespace ICARION
```

#### 2. Implement Force (`src/core/physics/forces/YourForce.cpp`)

```cpp
#include "YourForce.h"
#include "core/utils/mathUtils.h"

namespace ICARION {
namespace physics {

YourForce::YourForce(const config::DomainConfig& domain, double additional_param)
    : domain_(domain)  // ✅ Store reference
    , additional_param_(additional_param) {
}

Vec3 YourForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    // ✅ Read config directly from domain_
    if (!domain_.enable_some_feature) {
        return Vec3{0, 0, 0};
    }
    
    // Example: Use domain config parameters
    double param = domain_.some_config_value;
    Vec3 force = {0, 0, 0};
    
    // ... implement physics using domain_ and ctx ...
    
    return force;
}

} // namespace physics
} // namespace ICARION
```

#### 3. Create Unit Tests (`tests/physics/forces/test_your_force.cpp`)

```cpp
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/physics/forces/YourForce.h"
#include "core/config/DomainConfig.h"

using namespace ICARION::physics;
using namespace ICARION::config;
using Catch::Matchers::WithinAbs;

TEST_CASE("YourForce - Basic functionality", "[forces][yourforce]") {
    // ✅ Create config (SSOT)
    DomainConfig domain;
    domain.enable_some_feature = true;
    domain.some_config_value = 1.0;
    
    YourForce force(domain, 0.5);  // ✅ Pass by reference
    
    ForceContext ctx;
    ctx.domain = domain;  // ✅ Same config in context
    ctx.temperature_K = 300.0;
    
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{100, 0, 0};
    ion.ion_charge_C = 1.602e-19;
    ion.mass_kg = 100 * 1.66e-27;
    
    Vec3 F = force.compute(ion, 0.0, ctx);
    
    // Add your assertions
    REQUIRE(std::isfinite(F.x));
    REQUIRE(std::isfinite(F.y));
    REQUIRE(std::isfinite(F.z));
}
```

#### 4. Add to CMake (`tests/physics/forces/CMakeLists.txt`)

```cmake
add_executable(test_your_force
    test_your_force.cpp
)

target_link_libraries(test_your_force
    PRIVATE
        icarion_core
        Catch2::Catch2WithMain
)

add_test(NAME YourForce COMMAND test_your_force)
```

#### 5. Register in ForceRegistry (User Code)

```cpp
// ✅ Config lives in scope
config::DomainConfig domain = load_config("config.json");

ForceRegistry registry;

// ✅ Force stores reference to domain (no duplication!)
registry.add_force(std::make_unique<YourForce>(domain, 123.45));

// ⚠️ domain must outlive registry!
```

### Best Practices

✅ **Always provide parameter structs** (even if deprecated for Phase 2)  
✅ **Add deprecation warnings** to parameter structs pointing to SSOT violation  
✅ **Write comprehensive unit tests** (aim for >90% coverage)  
✅ **Document physics equations** in class docstrings  
✅ **Use const correctness** (`compute()` must be const)  
✅ **Check for NaN/Inf** in force output  

❌ **Don't mutate ion state** in `compute()`  
❌ **Don't store mutable state** in force objects  
❌ **Don't allocate in hot loops** (pre-allocate in constructor)

---

## Adding New Collision Models

### Overview

Collision damping is implemented in `DampingForce` with multiple models:
- **Friction**: Simple γ·m·v drag
- **HSD**: Elastic collisions with momentum transfer
- **Langevin**: Ion-neutral polarization interactions

### Adding a New Model

#### 1. Add to DampingModel Enum (`src/core/physics/forces/DampingForce.h`)

```cpp
enum class DampingModel {
    None,
    HSD,
    Langevin,
    Friction,
    YourNewModel  // Add here
};
```

#### 2. Add Parameters to DampingParams

```cpp
struct DampingParams {
    DampingModel model = DampingModel::None;
    
    // ... existing parameters ...
    
    // --- YourNewModel parameters ---
    double your_parameter = 0.0;
};
```

#### 3. Implement in DampingForce::compute()

```cpp
Vec3 DampingForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    // ... existing code ...
    
    switch (params_.model) {
        case DampingModel::YourNewModel: {
            gamma = compute_your_model_gamma(ion, ctx);
            break;
        }
        // ... other cases ...
    }
    
    return Vec3(-gamma * ion.mass_kg * ion.vel.x,
                -gamma * ion.mass_kg * ion.vel.y,
                -gamma * ion.mass_kg * ion.vel.z);
}
```

#### 4. Add Helper Method

```cpp
double DampingForce::compute_your_model_gamma(const IonState& ion, const ForceContext& ctx) const {
    // Your damping coefficient calculation
    return some_value;
}
```

#### 5. Add Tests

```cpp
TEST_CASE("DampingForce - YourNewModel", "[forces][damping]") {
    DampingParams params;
    params.model = DampingModel::YourNewModel;
    params.your_parameter = 123.0;
    
    DampingForce force(params);
    
    // Test your model
}
```

---

## Adding New Instrument Types

### Overview

Instrument-specific electric field calculations are in `ElectricFieldForce`. Each instrument type has analytical field formulas.

### Step-by-Step Guide

#### 1. Add Instrument to Enum (`src/instrument/InstrumentTypes.h`)

```cpp
enum class InstrumentType {
    LQIT = 0,
    IMS = 1,
    Orbitrap = 2,
    QuadrupoleRF = 3,
    TOF = 4,
    FTICR = 5,
    YourInstrument = 8,  // Add here
    // ...
};
```

⚠️ **Note:** This file violates SSOT and will be moved to `core/config/types/` in Phase 2.

#### 2. Add Parameters to AnalyticalFieldParams

```cpp
struct AnalyticalFieldParams {
    InstrumentType instrument_type = InstrumentType::UnknownInstrument;
    
    // ... existing parameters ...
    
    // YourInstrument-specific
    double your_voltage_V = 0.0;
    double your_frequency_Hz = 0.0;
};
```

#### 3. Implement Field Calculation in ElectricFieldForce

```cpp
Vec3 ElectricFieldForce::compute_analytical(const IonState& ion, double t) const {
    switch (params_.instrument_type) {
        case InstrumentType::YourInstrument:
            return compute_your_instrument_field(ion, t);
        // ... other cases ...
    }
}

Vec3 ElectricFieldForce::compute_your_instrument_field(const IonState& ion, double t) const {
    // Your electric field formula: E(x, y, z, t)
    
    double Ex = /* ... */;
    double Ey = /* ... */;
    double Ez = /* ... */;
    
    return Vec3{Ex, Ey, Ez};
}
```

#### 4. Add Tests

```cpp
TEST_CASE("ElectricFieldForce - YourInstrument", "[forces][electric]") {
    AnalyticalFieldParams params;
    params.instrument_type = InstrumentType::YourInstrument;
    params.your_voltage_V = 1000.0;
    params.your_frequency_Hz = 1e6;
    
    ElectricFieldForce force(params);
    
    // Test field at various positions
    IonState ion;
    ion.pos = Vec3{0.001, 0, 0};
    
    Vec3 E = force.compute(ion, 0.0, ForceContext{});
    
    // Verify field properties (symmetry, magnitude, direction)
}
```

#### 5. Document Physics

Add to class docstring or separate documentation:

```cpp
/**
 * YourInstrument Field:
 * 
 * Electric field: E(r, t) = ...
 * 
 * Parameters:
 * - your_voltage_V: Applied voltage [V]
 * - your_frequency_Hz: Drive frequency [Hz]
 * 
 * Assumptions:
 * - Ideal geometry (no fringe fields)
 * - Linear approximation near axis
 */
```

---

## Testing Guidelines

### Test Organization

```
tests/physics/forces/
├── test_force_registry.cpp          # ForceRegistry functionality
├── test_electric_field_force.cpp    # All instrument types
├── test_magnetic_damping_forces.cpp # Magnetic + damping
├── test_space_charge_force.cpp      # Space charge
├── test_force_integration.cpp       # Multi-force scenarios
└── minimal/                         # Standalone debug tests
    ├── minimal_electric_force.cpp
    ├── minimal_magnetic_force.cpp
    ├── minimal_damping_force.cpp
    └── minimal_space_charge_force.cpp
```

### Test Categories

1. **Unit Tests** (Catch2): Fine-grained component testing
2. **Integration Tests**: Multi-force scenarios
3. **Minimal Tests**: Standalone physics validation (no framework)

### Running Tests

```bash
# All tests
cd build && ctest

# Force tests only
ctest -R "Force"

# Specific test
./tests/physics/forces/test_electric_field_force

# Minimal tests (no CMake needed)
cd tests/physics/forces/minimal
./minimal_electric_force
```

### Writing Good Tests

```cpp
TEST_CASE("Descriptive name explaining what is tested", "[tag1][tag2]") {
    SECTION("Specific scenario") {
        // Arrange
        YourForce force(params);
        IonState ion = create_test_ion(/*...*/);
        
        // Act
        Vec3 F = force.compute(ion, 0.0, ctx);
        
        // Assert
        REQUIRE_THAT(F.x, WithinAbs(expected_Fx, 1e-10));
        REQUIRE(F.y < 0.0);  // Qualitative check
    }
}
```

**Key Principles:**
- Test one thing per section
- Use descriptive names
- Include edge cases (zero, infinity, boundary conditions)
- Verify both magnitude and direction
- Check for NaN/Inf
- Test symmetries and conservation laws

---

## Code Style and Conventions

### Naming

- **Classes:** `PascalCase` (e.g., `ElectricFieldForce`)
- **Functions:** `snake_case` (e.g., `compute_field()`)
- **Variables:** `snake_case` (e.g., `ion_charge_C`)
- **Constants:** `UPPER_SNAKE_CASE` (e.g., `ELEM_CHARGE_C`)
- **Members:** `trailing_underscore_` (e.g., `params_`)

### Units

Always include units in variable names:
```cpp
double voltage_V;           // Volts
double frequency_Hz;        // Hertz
double mass_kg;             // Kilograms
double position_m;          // Meters
double CCS_m2;              // m² (collision cross-section)
double time_s;              // Seconds
```

### Documentation

Use Doxygen-style comments:

```cpp
/**
 * @brief Brief description (one line)
 * 
 * @param ion The ion state
 * @param t Simulation time [s]
 * @param ctx Force computation context
 * @return Force vector [N]
 * 
 * Detailed description with physics equations and assumptions.
 */
Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const;
```

### SSOT Violations

When creating parameter structs (Phase 1 pattern), always add deprecation warning:

```cpp
// ============================================================================
// ⚠️ DEPRECATED: YourParams violates SSOT principle!
// ============================================================================
// This struct duplicates parameters from FullConfig → ...
// 
// **TODO (Phase 2):** Replace with direct config reference.
// **KEPT FOR NOW:** To avoid breaking changes during Phase 1.
// ============================================================================

struct YourParams {
    // ...
};
```

### Error Handling

```cpp
// Constructor validation
if (invalid_condition) {
    throw std::invalid_argument("ClassName: Clear error message explaining what's wrong");
}

// Runtime checks (hot path)
assert(ion.mass_kg > 0 && "Ion mass must be positive");

// Defensive checks
if (!std::isfinite(force.x) || !std::isfinite(force.y) || !std::isfinite(force.z)) {
    throw std::runtime_error("Force calculation produced non-finite result");
}
```

---

## Future Enhancements (Phase 2+)

### Planned Changes

1. **SSOT Refactor**: Remove all parameter structs, use `FullConfig` directly
2. **Move InstrumentType**: `instrument/InstrumentTypes.h` → `core/config/types/`
3. **Stochastic Forces**: Separate random kicks from deterministic forces
4. **GPU Acceleration**: CUDA kernels for space charge and field evaluation
5. **Field Caching**: Pre-compute and interpolate fields for performance

### Migration Strategy

When Phase 2 begins:
1. Update force constructors to accept config references
2. Deprecate parameter structs with compile warnings
3. Update all call sites incrementally
4. Remove parameter structs after full migration
5. Update this guide with new patterns

---

## Getting Help

- **Architecture Questions**: See `ARCHITECTURE.md`
- **API Reference**: See `docs/CLI_USAGE.md`
- **Config Schema**: See `schema/icarion-config.schema.json`
- **Bug Reports**: Create GitHub issue with minimal reproducer

---

**Document Status:** Living document, updated with each major refactoring phase.

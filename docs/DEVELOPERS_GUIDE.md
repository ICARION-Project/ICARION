# ICARION Developer's Guide

**Version:** 1.0 
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

**Version History:**
- **Pre-v1.0**: Forces used parameter structs (MagneticFieldParams, AnalyticalFieldParams, etc.) - DEPRECATED
- **v1.0+**: Forces use **const config references** (SSOT pattern)

### Design Principle: Single Source of Truth (SSOT)

****IMPORTANT (v1.0)**: Forces now use **const config references**, not parameter structs.

**MODERN (v1.0)**: Forces receive **const references** to config objects

**Why SSOT?**

- **No data duplication**: Config changes automatically propagate
- **Type safety**: Use strongly-typed config structs
- **Maintainability**: Single place to update parameters
- **Performance**: No copying of config data

---

### Step-by-Step Guide

#### 1. Create Force Header (`src/core/physics/forces/YourForce.h`)

**CORRECT PATTERN**: Use **const reference** to config, not parameter struct!

```cpp
#pragma once

#include "IForce.h"
#include "ForceContext.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "core/config/DomainConfig.h"  // Direct config reference

namespace ICARION {
namespace physics {

/**
 * @brief Your force description
 * 
 * Physics: F = ... (formula)
 * 
 * Uses const reference to DomainConfig (SSOT pattern)
 */
class YourForce : public IForce {
public:
    /**
     * @brief Constructor
     * 
     * @param domain Reference to simulation domain config (SSOT)
     * @param additional_param Any non-config parameters
     * 
     * **Config reference must outlive this object!
     */
    explicit YourForce(const config::DomainConfig& domain, double additional_param = 0.0);
    
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override;
    
private:
    const config::DomainConfig& domain_;  ///< Reference, not copy!
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
    : domain_(domain)  // Store reference
    , additional_param_(additional_param) {
}

Vec3 YourForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    // Read config directly from domain_
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
    // Create config (SSOT)
    DomainConfig domain;
    domain.enable_some_feature = true;
    domain.some_config_value = 1.0;
    
    YourForce force(domain, 0.5);  // Pass by reference
    
    ForceContext ctx;
    ctx.domain = domain;  // Same config in context
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
// Config lives in scope
config::DomainConfig domain = load_config("config.json");

ForceRegistry registry;

// Force stores reference to domain (no duplication!)
registry.add_force(std::make_unique<YourForce>(domain, 123.45));

// **domain must outlive registry!
```

### Best Practices (v1.0)

**DO:**

- **Use const config references**, not parameter structs
- **Store references as members**: `const config::DomainConfig& domain_;`
- **Read config on-demand**: `double V = domain_.fields.dc.axial_V;`
- **Write comprehensive unit tests** (aim for >90% coverage)
- **Document physics equations** in class docstrings
- **Use const correctness** (`compute()` must be const)
- **Check for NaN/Inf** in force output
- **Validate config in constructor** (throw if invalid)

****DON'T:**

- **Don't create parameter structs** (violates SSOT!)
- **Don't copy config data** (use references!)
- **Don't mutate ion state** in `compute()`
- **Don't store mutable state** in force objects
- **Don't allocate in hot loops** (pre-allocate in constructor)
- **Don't use raw pointers** (use const references or shared_ptr)


---

## Adding New Collision Models

### Overview

ICARION supports two types of collision models:

1. **Deterministic:** Continuous damping forces (use `DampingForce`)
2. **Stochastic:** Discrete collision events (use `ICollisionHandler`)

**IMPORTANT:** This guide covers stochastic models. For deterministic models, see "Adding New Force Types".

### Design Principle: SSOT Compliance

**DO:** Read directly from `EnvironmentConfig`  
**DO:** Store `const config::EnvironmentConfig&` references  
****DON'T:** Create parameter structs (SSOT violation!)  
****DON'T:** Copy config data into handler objects

### Step-by-Step Guide

#### 1. Implement ICollisionHandler

**File:** `src/core/physics/collisions/YourCollisionHandler.h`

```cpp
#pragma once

#include "ICollisionHandler.h"

namespace ICARION::physics {

/**
 * @brief Your collision model description
 * 
 * Physics: [Describe collision mechanism]
 * 
 * SSOT Design: Reads parameters directly from EnvironmentConfig.
 */
class YourCollisionHandler : public ICollisionHandler {
public:
    /**
     * @brief Construct handler
     * @param your_param Model-specific parameter (if needed)
     */
    explicit YourCollisionHandler(double your_param = 1.0);
    
    bool handle_collision(
        IonState& ion,
        double dt,
        EhssRng& rng,
        const config::EnvironmentConfig& env  // SSOT!
    ) override;
    
    std::string name() const override { return "YourModel"; }
    
private:
    double your_param_;
};

} // namespace ICARION::physics
```

#### 2. Implement Logic

**File:** `src/core/physics/collisions/YourCollisionHandler.cpp`

```cpp
#include "YourCollisionHandler.h"
#include "utils/constants.h"

namespace ICARION::physics {

YourCollisionHandler::YourCollisionHandler(double your_param)
    : your_param_(your_param)
{}

bool YourCollisionHandler::handle_collision(
    IonState& ion,
    double dt,
    EhssRng& rng,
    const config::EnvironmentConfig& env
) {
    // Read parameters directly from env (SSOT!)
    const double T_K = env.temperature_K;
    const double n = env.particle_density_m_3;
    const double m_neutral = env.neutral_mass_kg;
    const Vec3 v_gas = env.gas_velocity_m_s;
    
    // Your collision physics here
    // ...
    
    // Modify ion velocity if collision occurs
    if (collision_occurred) {
        ion.vel = new_velocity;
        return true;
    }
    
    return false;
}

} // namespace ICARION::physics
```

#### 3. Add to CollisionHandlerFactory

**File:** `src/core/physics/collisions/CollisionHandlerFactory.cpp`

```cpp
#include "YourCollisionHandler.h"

std::unique_ptr<ICollisionHandler> CollisionHandlerFactory::create(
    const config::PhysicsConfig& config,
    const GeometryMap* geometry_map,
    double gamma_for_ou
) {
    switch (config.collision_model) {
        // ... existing cases ...
        
        case CM::YourModel:  // Add your model
            return std::make_unique<YourCollisionHandler>(config.your_param);
        
        default:
            throw std::runtime_error("Unknown collision model");
    }
}
```

#### 4. Add Enum Value

**File:** `src/core/config/types/PhysicsEnums.h`

```cpp
enum class CollisionModel {
    // ... existing values ...
    YourModel,         // Add your model
    UnknownCollisionModel
};
```

**File:** `src/core/config/conversion/EnumMapper.h`

```cpp
inline CollisionModel collision_model_from_string(std::string_view str) {
    // ... existing mappings ...
    if (str_lower == "yourmodel") return CollisionModel::YourModel;
    // ...
}

inline std::string collision_model_to_string(CollisionModel model) {
    // ... existing cases ...
    case CollisionModel::YourModel: return "YourModel";
    // ...
}
```

#### 5. Write Unit Tests

**File:** `tests/physics/collisions/test_your_collision_handler.cpp`

```cpp
#include <catch2/catch_test_macros.hpp>
#include "YourCollisionHandler.h"

TEST_CASE("YourCollisionHandler: Basic functionality", "[collision][yourmodel]") {
    // Setup
    physics::YourCollisionHandler handler(1.0);
    
    IonState ion;
    ion.mass_kg = 50.0 * 1.66e-27;
    ion.vel = Vec3{100, 0, 0};
    
    config::EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    // ...
    
    EhssRng rng(42);
    double dt = 1e-9;
    
    // Act
    bool collision = handler.handle_collision(ion, dt, rng, env);
    
    // Assert
    if (collision) {
        REQUIRE(ion.vel.magnitude() < 100.0);  // Velocity changed
    }
}
```

#### 6. Update CMakeLists.txt

```cmake
# tests/physics/collisions/CMakeLists.txt

add_executable(test_your_collision_handler
    test_your_collision_handler.cpp
)
target_link_libraries(test_your_collision_handler PRIVATE
    icarion_core
    Catch2::Catch2WithMain
)
add_test(NAME YourCollisionHandler COMMAND test_your_collision_handler)
```

### Best Practices

1. **SSOT Compliance:** Always read from `EnvironmentConfig` directly
2. **Return Value:** Return `true` if collision occurred, `false` otherwise
3. **In-Place Modification:** Modify `ion.vel` directly (don't return new velocity)
4. **Validation:** Check for invalid parameters in constructor
5. **Logging:** Use `io::debug_log()` for warnings/diagnostics
6. **Testing:** Test edge cases (zero velocity, zero density, etc.)

### Common Patterns

#### Pattern 1: Collision Probability

```cpp
// Mean free path approach
double collision_rate = n * sigma * v_rel_mag;
double P = 1.0 - std::exp(-collision_rate * dt);

if (rng.uniform01() < P) {
    // Apply collision
}
```

#### Pattern 2: Post-Collision Velocity

```cpp
// Thermal scattering
double v_thermal = std::sqrt(BOLTZMANN_CONSTANT * T_K / reduced_mass);

// Random direction
double theta = std::acos(2.0 * rng.uniform01() - 1.0);
double phi = 2.0 * M_PI * rng.uniform01();

Vec3 v_new;
v_new.x = v_thermal * std::sin(theta) * std::cos(phi);
v_new.y = v_thermal * std::sin(theta) * std::sin(phi);
v_new.z = v_thermal * std::cos(theta);

ion.vel = v_new + v_gas;  // Add gas velocity
```

### Troubleshooting

**Problem: Ions equilibrate to ~70-80% of thermal energy instead of 100%**

**Error:** Kinetic energy does not equilibrate correctly to gas temperature. This effect is more or less independent of temperature, and persists across different random seeds and time steps, as well as neutral gas molecules (with different masses).

**Cause:** Collision probability calculated with bulk gas velocity (`v_rel = ion.vel - v_gas`), but collision applied with different Maxwell-Boltzmann sampled neutral. This velocity inconsistency creates systematic bias.

**Solution:** Sample neutral velocity FIRST, then calculate collision probability using the SAME sampled velocity. See commit 92d29c1 for implementation.

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

****Note:** This file violates SSOT and will be moved to `core/config/types/` in the future.

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

### SSOT Compliance

All new code must follow the Single Source of Truth principle:

```cpp
// CORRECT: Use const config reference
class MyForce : public IForce {
public:
    MyForce(const config::DomainConfig& domain);
    
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override {
        // Read config on-demand
        double axial_V = domain_.fields.dc.axial_V;
        // ...
    }
    
private:
    const config::DomainConfig& domain_;  // Reference (SSOT)
};
```

```cpp
// INCORRECT: Parameter struct (SSOT violation)
struct MyForceParams {
    double axial_V;  // Duplicates domain.fields.dc.axial_V
};

class MyForce : public IForce {
public:
    MyForce(const MyForceParams& params);  // Will be rejected
    
private:
    MyForceParams params_;  // SSOT violation
};
```

Rationale:

- Single Source of Truth: Config changes propagate automatically
- Type Safety: Compiler enforces correct config structure
- No Duplication: Zero-copy, references only
- Maintainability: One place to change parameters

Migration Note: All legacy parameter structs (MagneticFieldParams, AnalyticalFieldParams, DampingParams) were deleted in ICARION v1.0. New code must use direct config references.

See: src/core/physics/forces/ElectricFieldForce.h (reference implementation)

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

## Future Enhancements

### In Progress

1. **Reaction System Database Unification** (Phase 3D, approximately 4-6h):
   - Unify species types (remove ICARION::io::Species, reactionUtils::Species)
   - Wire reaction_handler directly into integrator
   - Delete legacy reaction loading code
   - Blocker: Type mismatch between species databases

2. **SimulationEngine Integration** (Phase 5A, approximately 12h):
   - Main simulation loop using integration strategies
   - Orchestration of ForceRegistry + CollisionHandler + ReactionHandler
   - Boundary condition handling
   - Output management (HDF5Writer v2)
   - Status: Ready to start (Phase 4 complete)

### Planned

1. **InstrumentType Location** (Low priority, approximately 30min):
   - Move instrument/InstrumentTypes.h to core/config/types/InstrumentType.h
   - Eliminates dependency cycle config to instrument

2. **Stochastic Forces** (Design phase):
   - Separate random kicks from deterministic forces
   - New IStochasticForce interface (similar to collision handlers)

3. **Space Charge** ✅ (COMPLETE - November 2025):
   - Automatic method selection: N<1000→Direct (O(N²)), N≥1000→Grid (O(N log N))
   - CIC charge deposition with O(h²) convergence
   - Poisson solver with 5 methods (Gauss-Seidel, SOR, CG, Multigrid, FFT)
   - Files: `src/core/physics/spacecharge/*`, `src/core/physics/forces/SpaceCharge{Direct,Grid}.{h,cpp}`
   - Tests: 17 unit tests, 3 integration tests (all passing)
   - **When to use Direct:** N<1000, exact results needed
   - **When to use Grid:** N≥1000, fast approximation
   - **Configuration:** Set `physics.enable_space_charge = true` in config
   - **Performance tuning:** Adjust grid size (default 64³), update frequency
   - **Limitations:** CPU-only (GPU planned for v1.1), Dirichlet BC causes errors near boundaries

4. **GPU Acceleration** (Planned for v1.1):
   - CUDA kernels for space charge Poisson solver
   - Grid-based field evaluation on GPU
   - Target: 10-100x speedup for N≥10,000

5. **Field Caching** (Optimization phase):
   - Pre-compute field on regular grid
   - Trilinear interpolation for fast evaluation
   - Benefit: Approximately 10x speedup for analytical fields

### Completed

- **Force System SSOT** (Phase 1, Steps 1-4 complete)
  - IForce interface with ForceRegistry
  - ElectricFieldForce, MagneticFieldForce, DampingForce, SpaceChargeForce
  - Full unit test coverage (27/27 tests passing)

- **Collision System SSOT** (Phase 2C complete)
  - ICollisionHandler interface with factory
  - EHSS, HSS, OU collision handlers
  - Energy conservation validation

- **Reaction System Handlers** (Phase 3C complete)
  - IReactionHandler interface with factory
  - StochasticReactionHandler implementation
  - Database-driven reaction loading

- **Integration Strategies** (Phase 4A/4B complete, November 2025)
  - IIntegrationStrategy interface
  - RK4Strategy (4th-order Runge-Kutta, fixed-step)
  - RK45Strategy (Dormand-Prince 5(4), adaptive timestep with FSAL)
  - BorisStrategy (symplectic pusher for electromagnetic fields)
  - IntegrationStrategyFactory for runtime selection
  - Comprehensive test suites (27/27 tests passing, 100%)
  - SSOT-compliant (uses DomainConfig pointer, callback-based acceleration)
  - **Strategies tested standalone, NOT yet integrated into main.cpp**
  - Files: `src/core/integrator/strategies/*`

**Next Step:** Phase 5 (SimulationEngine) will integrate strategies into main simulation loop

---

## Getting Help

- **Architecture Questions**: See `ARCHITECTURE.md`
- **API Reference**: See `docs/CLI_USAGE.md`
- **Config Schema**: See `schema/icarion-config.schema.json`
- **Bug Reports**: Create GitHub issue with minimal reproducer

---

**Document Status:** Living document, updated with each major refactoring phase.

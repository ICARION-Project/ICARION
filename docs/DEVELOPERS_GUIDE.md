# ICARION Developer's Guide

**Version:** 1.0 
**Last Updated:** December 2025

This guide provides practical instructions for extending ICARION with new features.  
**Baseline:** Core CPU paths are SoA-first (`IonEnsemble`). AoS helpers exist only for tests/legacy entry points and should not be used in hot loops.

**Separation of concerns:** The engine orchestrates while pluggable strategies/registries do the work: forces (`IForce` + `ForceRegistry`), collisions (`ICollisionHandler` + factory), reactions (`IReactionHandler` + factory), fields/geometry (`IFieldModel`/`IDomainGeometry` via DomainManager), integrators (`IIntegrationStrategy` + factory). Swap components without changing the engine loop.

---

## Table of Contents

1. [Adding New Force Types](#adding-new-force-types)
2. [Adding New Collision Models](#adding-new-collision-models)
3. [Adding New Instrument Types](#adding-new-instrument-types)
4. [GPU Development Guide](#gpu-development-guide)
5. [Testing Guidelines](#testing-guidelines)
6. [Code Style and Conventions](#code-style-and-conventions)

---

## Adding New Force Types

### Overview

ICARION's force system follows a plugin architecture using the **IForce interface**. All forces implement a single SoA entry point `IForce::compute(...)` (ensemble + ion index + time + context) and are managed by `ForceRegistry`. AoS hooks have been removed from the hot path; tests wrap single ions into a scratch SoA when needed.  

**New in v1.0 SoA path:**  
- `IForce::compute_soa(const ForceState&, t, ctx)` is required; it must use the provided state (no ensemble fetches).  
- `ForceRegistry::compute_total_force_soa` aggregates via `compute_soa` and uses `ForceState::ensemble_index` for space-charge (no AoS staging).  
- RK45 calls only the SoA path; no scratch ensembles are built in integration stages.

**Version:** 1.0 uses **const config references** (Single Source of Truth pattern)

### Design Principle: Single Source of Truth (SSOT)

**IMPORTANT**: Forces use **const config references**, not parameter structs.

**Pattern**: Forces receive **const references** to config objects

**Why SSOT?**

- **No data duplication**: Config changes automatically propagate
- **Type safety**: Use strongly-typed config structs
- **Maintainability**: Single place to update parameters
- **Performance**: No copying of config data

---

### Step-by-Step Guide (SoA-first)

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
    
Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override;  // AoS (tests)
Vec3 compute(const core::IonEnsemble& ensemble, size_t i, double t, const ForceContext& ctx) const override; // SoA hot path
    
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
    domain.name = "test_domain";
    domain.instrument = InstrumentType::NoFixedInstrument;
    
    YourForce force(domain, 0.5);  // Pass by reference
    
    ForceContext ctx;
    ctx.domain = &domain;  // Pointer to config
    
    IonState ion;
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{100, 0, 0};
    ion.ion_charge_C = 1.602e-19;  // Elementary charge
    ion.mass_kg = 100.0 * 1.66e-27;  // 100 amu
    
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

- Implement `compute` (SoA-only); tests may wrap single ions into a scratch ensemble if needed.
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
        PhysicsRng& rng,
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
    PhysicsRng& rng,
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
    ion.pos = Vec3{0, 0, 0};
    ion.vel = Vec3{100, 0, 0};
    ion.mass_kg = 50.0 * 1.66e-27;  // 50 amu
    ion.ion_charge_C = 1.602e-19;
    
    config::EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    // ...
    
    PhysicsRng rng(42);
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

**Problem: Ions starting at domain boundaries are immediately deactivated**

**Error:** Ions placed exactly at or very close to entrance boundaries (e.g., z=0 for domain starting at z=0) drift backwards after first collision and are instantly eliminated. This is especially problematic at low pressures where random thermal velocities from collisions can exceed drift velocities.

**Symptoms:**
- Ions with position z ≈ 0 and negative velocity (vz < 0) after first collision
- Ion marked inactive after ~1-2 collisions
- Problem worse at low pressures (< 1000 Pa) where thermal diffusion dominates
- Example: IMS drift tube with ion at z=0, domain origin at (0,0,0)

**Cause:** 
1. Ion starts exactly at entrance boundary (z = 0)
2. Boundary check: ion is active if `z >= -EPSILON` (typically -1e-12)
3. First collision randomizes velocity with thermal component
4. Even small backward velocity component → `z < -EPSILON` → ion deactivated
5. Effect amplified at low pressure: fewer collisions → larger thermal kicks

**Solution:** **Shift domain origin backward** to create entrance buffer zone:

```cpp
// For a 10cm IMS drift tube with ions starting at z=0:
config::DomainConfig domain;
domain.geometry.origin_m = Vec3{0.0, 0.0, -0.02};  // Shift 2cm backward
domain.geometry.length_m = 0.10 + 0.04;             // Add 4cm buffer (2cm entrance + 2cm exit)
domain.geometry.radius_m = 0.5;                     // Wide radius to prevent radial losses

// Ion starting position remains z=0 (relative to lab frame)
// But now ion is 2cm INSIDE the domain (domain extends from z=-0.02 to z=0.12)
```

**Why this works:**
- Ion at z=0 is now safely inside domain (domain starts at z=-0.02)
- Backward collisions at z=0 → z=-0.001 are still within domain bounds
- Provides safety margin for thermal diffusion

**Best Practices:**
- **Rule of thumb**: Add buffer distance ≥ 2× mean free path
- **Low pressure (< 100 Pa)**: Use buffer ≥ 2-5 cm
- **High pressure (> 10 kPa)**: Buffer ≥ 0.5-1 cm usually sufficient
- **Wide radius**: Use `radius_m ≥ 5 × length_m` to prevent radial losses from diffusion

**Future Development:**
- **v1.1+**: Planned boundary types beyond current absorbing boundaries:
  - `BoundaryType::Reflecting` - elastic reflection at boundaries
  - `BoundaryType::Emitting` - continuous ion source at entrance
  - `BoundaryType::Periodic` - wraparound for bulk simulations
- **Current (v1.0)**: Only absorbing boundaries implemented

**Related Issues:**
- See `tests/instruments/test_ims_drift.cpp` for working example
- Commits: 3868379 (boundary fix), e00e8b0 (EHSS test), b9a57f4 (CTest fixes)

---

## Adding New Instrument Types

### Overview

Instrument-specific electric field calculations live in FieldModels and are consumed by `ElectricFieldForce`. Analytical formulas are implemented in `AnalyticalFieldModel`; grid/BEM/FEM fields use `FieldProviderModel` (wraps `IFieldProvider`). Field models are injected via `PhysicsSetup` into `ForceRegistry` (SSOT); `ElectricFieldForce` falls back only if none is provided.

Multi-domain geometry handling lives in `IDomainGeometry` strategies (e.g., `CylindricalGeometry`, `OrbitrapGeometry`) used by `DomainManager` and `DomainContext` for transforms and boundary checks. DomainManager now only orchestrates these strategies (no AoS boundary helpers); geometry classes encapsulate containment/intersection logic.

### Space-Charge Architecture (v1.0)

- `ISpaceChargeModel` exposes `update_fields()` + `sample_electric_field()`. ForceRegistry owns an optional instance and adds Coulomb force directly in the SoA loop (no fallback AoS conversions).
- Models:
  - `SpaceChargeDirectModel` – exact O(N²), shared across domains for small ion counts.
  - `SpaceChargeGridModel` – geometry-driven Poisson solver (Dirichlet masks + bounding boxes from `IDomainGeometry`).
  - `SpaceChargeGPUModel` – wraps `gpu::GPUSpaceChargeP3M`; compiles as a stub in CPU-only builds.
- `SpaceChargeModelFactory` decides per-domain: try GPU if `physics.enable_space_charge_gpu` and CUDA build, else grid, else direct. Logging records fallbacks automatically.
- Configuration overrides / CLI: `physics.enable_space_charge` toggles feature, `physics.enable_space_charge_gpu` requests GPU acceleration (safe to enable even on CPU because the factory degrades gracefully).

### Step-by-Step Guide

#### 1. Add Instrument to Enum (`src/core/config/types/InstrumentTypes.h`)

```cpp
enum class InstrumentType {
    LQIT = 0,
    IMS = 1,
    Orbitrap = 2,
    QuadrupoleRF = 3,
    TOF = 4,
    FTICR = 5,
    NoFixedInstrument = 6,
    UnknownInstrument = 7,
    YourInstrument = 8,  // Add here
    // ...
};
```

#### 2. Add Field Parameters to DomainConfig

Add instrument-specific field parameters to the appropriate config structures in `src/core/config/types/`.

#### 3. Implement Field Calculation in a FieldModel

Add your instrument field to `AnalyticalFieldModel` (or a new `IFieldModel` implementation):

```cpp
class YourFieldModel : public IFieldModel {
public:
    explicit YourFieldModel(const DomainConfig& dom) : dom_(&dom) {}
    Vec3 E(const Vec3& pos, double t) const override {
        // Your instrument-specific analytical field: E(x, y, z, t)
        double Ex = /* ... */;
        double Ey = /* ... */;
        double Ez = /* ... */;
        return Vec3{Ex, Ey, Ez};
    }
private:
    const DomainConfig* dom_;
};
```

#### 4. Add Tests

Add a parity or direct test for your FieldModel and ElectricFieldForce:

```cpp
TEST_CASE("YourFieldModel parity", "[forces][field]") {
    DomainConfig domain;
    domain.instrument = InstrumentType::YourInstrument;
    domain.fields.dc.axial_V = ValueOrWaveform(1000.0);
    // instrument-specific params...
    domain.finalize();

    YourFieldModel model(domain);
    ElectricFieldForce force(domain);

    IonState ion;
    ion.pos = Vec3{0.001, 0, 0};
    ion.mass_kg = ion.ion_charge_C = 1.0;

    ForceContext ctx;
    ctx.field_model = &model;
    ctx.domain = &domain;

    Vec3 F = force.compute(ion, 0.0, ctx);
    REQUIRE(std::isfinite(F.x));
    REQUIRE(std::isfinite(F.y));
    REQUIRE(std::isfinite(F.z));
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

## GPU Development Guide

**Added:** November 2025 (v1.0)
**Status:** Integration/space-charge/collision paths wired via factories; boundary helper remains experimental.

### Overview

ICARION's GPU acceleration is designed for **easy extensibility**. This guide shows how to add new GPU-accelerated features.

**Current GPU Features (v1.0):**
- RK4/RK45/Boris batch integrators via `GPUIntegrationStrategy` (factory-selected wrapper; automatic fallback + grid/E-field checks)
- HSS/EHSS collision helper (active-ion threshold default 5000; EHSS geometry upload TODO)
- Reactions via `GPUReactionBackend` (constant / Arrhenius / modified Arrhenius, mixtures). Flattens per-domain tables, launches XORWOW kernel, updates SoA species/mass/charge/CCS/mobility. Device buffers are pooled (RNG, species/domain/flags/reaction tables); CPU stochastic handler is the fallback.
- Field-provider upload for integration when ElectricFieldForce is present
- Space charge P³M helper wired through `SpaceChargeGPUModel` (opt-in via `physics.enable_space_charge_gpu`, CPU fallback guaranteed)
- Boundary check helper supports absorption/cylindrical only and is not wired into the main loop

`IntegrationStrategyFactory` emits the GPU wrapper automatically whenever `simulation.enable_gpu` is true (and CUDA is available); the wrapper calls `IIntegrationStrategy::step_batch()` internally and falls back to the CPU strategy if the batch conditions are not met.

**Note:** In v1.0 the runtime GPU path is disabled (CPU-only), even if `enable_gpu` is set. CUDA builds still compile GPU code for development; CPU fallback always wins.

### Prerequisites

**Build Requirements:**
- CUDA Toolkit 11.0+ (tested with CUDA 12.0)
- GPU with Compute Capability 7.5+ (Turing, Ampere, Ada, Hopper)
- CMake 3.16+

**Build GPU-enabled ICARION:**
```bash
cmake -B build -DUSE_GPU_ACCEL=ON
make -C build -j$(nproc)
```

**Check GPU availability at runtime:**
```cpp
#include "core/gpu/GPUContext.h"

if (icarion::gpu::GPUContext::is_cuda_available()) {
    auto context = icarion::gpu::GPUContext::create(0);  // Device 0
    std::cout << "GPU: " << context->get_properties().name << "\n";
}
```

### Adding a New GPU Kernel

#### Example: Adding GPU Collision Kernel

**Step 1: Create Kernel File (`src/core/gpu/collision_batch.cu`)**

```cuda
#include "collision_batch.cuh"
#include "core/gpu/GPUContext.h"

namespace icarion {
namespace gpu {

// Device function: Compute collision for one ion
__device__ bool check_collision_gpu(
    const Vec3& pos,
    const Vec3& vel,
    double mass,
    double charge,
    double dt,
    // ... collision parameters
) {
    // Your collision logic here
    return collision_occurred;
}

// Kernel: Process all ions in parallel
__global__ void collision_batch_kernel(
    // Input/output ion state arrays (SoA)
    const double* x, const double* y, const double* z,
    const double* vx, const double* vy, const double* vz,
    double* vx_out, double* vy_out, double* vz_out,
    const double* mass, const double* charge,
    const bool* active, bool* collision_flags,
    // Collision parameters
    double dt, int N,
    // ... collision model params
) {
    // Grid-stride loop (handles any N)
    for (int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < N;
         i += gridDim.x * blockDim.x) {
        
        if (!active[i]) continue;
        
        Vec3 pos = {x[i], y[i], z[i]};
        Vec3 vel = {vx[i], vy[i], vz[i]};
        
        // Check collision
        bool collided = check_collision_gpu(
            pos, vel, mass[i], charge[i], dt, /* params */
        );
        
        if (collided) {
            collision_flags[i] = true;
            // Update velocity
            vx_out[i] = /* new vx */;
            vy_out[i] = /* new vy */;
            vz_out[i] = /* new vz */;
        }
    }
}

// Host function: Launch kernel
void collision_batch(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    /* collision params */,
    cudaStream_t stream
) {
    int N = ions_in.count;
    int threads = 256;
    int blocks = std::min((N + threads - 1) / threads, 2048);
    
    collision_batch_kernel<<<blocks, threads, 0, stream>>>(
        ions_in.pos_x, ions_in.pos_y, ions_in.pos_z,
        ions_in.vel_x, ions_in.vel_y, ions_in.vel_z,
        ions_out.vel_x, ions_out.vel_y, ions_out.vel_z,
        ions_in.mass_kg, ions_in.charge_C,
        ions_in.active, /* collision_flags */,
        dt, N, /* params */
    );
    
    CUDA_CHECK(cudaGetLastError());
}

} // namespace gpu
} // namespace icarion
```

**Collision Path Wiring (as implemented)**

- `ICollisionHandler` exposes optional batch hooks:
  - `supports_batch()` advertises accelerator support.
  - `handle_batch(...)` receives the SoA ensemble plus an index list for the current domain and either processes the entire batch (returning `true`) or requests CPU fallback (`false`).
- `GPUCollisionHandler` wraps an EHSS/HSS CPU handler, wires up `GPUCollisionHelper`, and implements the batch hook. It logs and falls back automatically if CUDA prerequisites are not met or if the batch is below the configured threshold.
- `SimulationEngine::perform_collisions(...)` groups indices per-domain after birth/domain detection, then:
  1. Calls `handle_batch(...)` when the handler advertises support.
  2. Falls back to `handle_collisions_cpu(...)` (per-ion RNG, SoA views) whenever the batch hook returns `false` or is unavailable.


- **Multi-gas support:** The CUDA kernels now sample neutral velocities per mixture component (matching the CPU algorithm), accumulate component-specific collision rates, and select the actual collision partner using the same Monte-Carlo scheme. For cache efficiency the device struct keeps a fixed-size mixture buffer (`MAX_GPU_GAS_COMPONENTS = 8`); additional components are truncated with a runtime warning.
```cpp
void SimulationEngine::perform_collisions(IonEnsemble& ensemble,
                                          double dt,
                                          const std::vector<int>& domains) {
    if (!collision_handler_) return;
    const bool has_batch = collision_handler_->supports_batch();
    for (size_t dom = 0; dom < config_.domains.size(); ++dom) {
        auto& indices = domain_buckets[dom];
        if (indices.empty()) continue;
        const auto& env = config_.domains[dom].environment;
        bool handled = has_batch &&
            collision_handler_->handle_batch(ensemble, indices, dt,
                                             env, rng_by_ion_);
        if (!handled) {
            handle_collisions_cpu(ensemble, dt, indices, env);
        }
    }
}
```

This keeps the SimulationEngine agnostic of GPU specifics—the handler encapsulates accelerator logic, and SoA data never round-trips through temporary AoS buffers.

### GPU Development Best Practices

#### 1. Memory Access Patterns

**GOOD: Coalesced Access (SoA)**
```cuda
// Threads access consecutive memory locations
for (int i = threadIdx.x; i < N; i += blockDim.x) {
    x[i] = compute(x[i]);  // Coalesced: x[0], x[1], x[2], ...
}
```

**BAD: Strided Access (AoS)**
```cuda
// Threads access non-consecutive memory
struct Ion { double x, y, z, vx, vy, vz; };
Ion* ions;
for (int i = threadIdx.x; i < N; i += blockDim.x) {
    ions[i].x = compute(ions[i].x);  // Strided: ions[0].x, ions[1].x, ...
}
```

#### 2. Kernel Launch Configuration

**Grid-Stride Loop Pattern:**
```cuda
__global__ void my_kernel(double* data, int N) {
    // Each thread processes multiple elements
    for (int i = blockIdx.x * blockDim.x + threadIdx.x;
         i < N;
         i += gridDim.x * blockDim.x) {
        
        data[i] = compute(data[i]);
    }
}

// Launch with limited blocks (avoids scheduler overhead)
int threads = 256;  // Good occupancy on most GPUs
int blocks = std::min((N + threads - 1) / threads, 2048);
my_kernel<<<blocks, threads>>>(data, N);
```

**Why Grid-Stride?**
- Works for any N (1 to 10M+)
- Optimal occupancy across GPU architectures
- Reduces scheduling overhead
- Better instruction cache utilization

#### 3. Error Handling

**Always check CUDA errors:**
```cuda
// After kernel launch
CUDA_CHECK(cudaGetLastError());

// After synchronization
CUDA_CHECK(cudaStreamSynchronize(stream));

// Asynchronous error check
cudaError_t err = cudaGetLastError();
if (err != cudaSuccess) {
    // Log error and fall back to CPU
    fprintf(stderr, "CUDA error: %s\n", cudaGetErrorString(err));
    return false;
}
```

#### 4. Async Pipeline

**Overlap CPU and GPU work:**
```cpp
// Upload (async)
ion_state_conversion::upload_ions(ions, ions_gpu_in, stream);

// Compute (async, queued after upload)
my_kernel<<<blocks, threads, 0, stream>>>(ions_gpu_in, ions_gpu_out);

// Download (async, queued after compute)
ion_state_conversion::download_ions(ions_gpu_out, ions, stream);

// Synchronize once at end
context.synchronize();
```

#### 5. Performance Profiling

**Use NVIDIA Nsight Systems:**
```bash
# Profile GPU-enabled run
nsys profile -o report ./icarion_main config.json

# View timeline
nsys-ui report.nsys-rep
```

**Key metrics to check:**
- Kernel duration
- Memory transfer time
- GPU occupancy
- Warp execution efficiency

### Testing GPU Code

#### Unit Tests

**Test GPU context:**
```cpp
TEST(GPUContext, Creation) {
    if (!GPUContext::is_cuda_available()) {
        GTEST_SKIP() << "CUDA not available";
    }
    
    auto context = GPUContext::create(0);
    ASSERT_NE(context, nullptr);
    EXPECT_GT(context->get_properties().total_memory, 0);
}
```

`GPUIntegrationStrategy` exercises the helper indirectly inside the main loop, but the helper can still be unit-tested in isolation:

**Test kernel correctness:**
```cpp
TEST(IntegrationKernel, SingleIon) {
    auto context = GPUContext::create(0);
    
    // CPU reference
    IonState ion_cpu = /* initial state */;
    rk4_step_cpu(ion_cpu, dt);
    
    // GPU result
    std::vector<IonState> ions = {ion_cpu};
    GPUIntegrationHelper helper(context, 1);  // Threshold = 1
    helper.integrate_batch_rk4(ions, dt, 0.0);
    
    // Compare (tolerance: 1e-12)
    EXPECT_NEAR(ions[0].pos.x, ion_cpu.pos.x, 1e-12);
}
```

#### Performance Tests

**Measure speedup:**
```cpp
TEST(GPUPerformance, LargeScale) {
    auto context = GPUContext::create(0);
    std::vector<IonState> ions(100000);
    
    // CPU baseline
    auto t0 = now();
    cpu_integrate(ions, dt);
    auto cpu_time = elapsed(t0);
    
    // GPU time
    t0 = now();
    gpu_helper->integrate_batch_rk4(ions, dt, 0.0);
    auto gpu_time = elapsed(t0);
    
    double speedup = cpu_time / gpu_time;
    EXPECT_GT(speedup, 1.0);  // Verify GPU is faster
}
```

### Common Pitfalls

#### Don't: Mix Host and Device Pointers
```cpp
double* host_ptr = new double[N];
double* device_ptr;
cudaMalloc(&device_ptr, N * sizeof(double));

// ERROR: Can't dereference device pointer on CPU
double x = device_ptr[0];  // Segfault!

// ERROR: Can't pass host pointer to kernel
my_kernel<<<...>>>(host_ptr, N);  // Invalid memory access!
```

#### Do: Use Explicit Transfer
```cpp
cudaMemcpy(device_ptr, host_ptr, N * sizeof(double), cudaMemcpyHostToDevice);
my_kernel<<<...>>>(device_ptr, N);
cudaMemcpy(host_ptr, device_ptr, N * sizeof(double), cudaMemcpyDeviceToHost);
```

#### Don't: Forget to Synchronize
```cpp
my_kernel<<<...>>>(data, N);
// Kernel is async! Results not ready yet!
use_results(data);  // Race condition!
```

#### Do: Synchronize Before Using Results
```cpp
my_kernel<<<...>>>(data, N);
cudaStreamSynchronize(stream);  // Wait for kernel
use_results(data);  // Safe now
```

### GPU Code Review Checklist

Before submitting GPU code, verify:

- [ ] Error checking after every CUDA call
- [ ] Grid-stride loop for scalability
- [ ] Coalesced memory access (SoA layout)
- [ ] Proper synchronization
- [ ] CPU fallback on GPU errors
- [ ] Unit tests with CPU reference
- [ ] Performance tests (speedup measurement)
- [ ] Documentation of kernel algorithm
- [ ] Conditional compilation (`#ifdef ICARION_USE_GPU`)

### Adding a New GPU Integrator

**Example: Adding Velocity Verlet GPU support**

#### Step 1: Create CUDA Kernel (`src/core/gpu/integrate_verlet_batch.cu`)

```cuda
#include "integrate_verlet_batch.cuh"
#include "core/gpu/compute_acceleration.cuh"  // Shared force eval

namespace icarion {
namespace gpu {

__global__ void integrate_verlet_batch_kernel(
    const IonStateGPU ions_in,
    IonStateGPU ions_out,
    Vec3 E_field, Vec3 B_field,
    double dt, int N
) {
    // Grid-stride loop
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
         i < N; 
         i += gridDim.x * blockDim.x) {
        
        if (!ions_in.active[i]) continue;
        
        // Load ion state
        Vec3 pos = {ions_in.pos_x[i], ions_in.pos_y[i], ions_in.pos_z[i]};
        Vec3 vel = {ions_in.vel_x[i], ions_in.vel_y[i], ions_in.vel_z[i]};
        double mass = ions_in.mass_kg[i];
        double charge = ions_in.charge_C[i];
        
        // Velocity Verlet: x_new = x + v*dt + 0.5*a*dt^2
        Vec3 acc = compute_acceleration(pos, vel, mass, charge, E_field, B_field);
        Vec3 pos_new = pos + vel * dt + acc * (0.5 * dt * dt);
        
        // v_new = v + 0.5*(a_old + a_new)*dt
        Vec3 acc_new = compute_acceleration(pos_new, vel, mass, charge, E_field, B_field);
        Vec3 vel_new = vel + (acc + acc_new) * (0.5 * dt);
        
        // Write outputs
        ions_out.pos_x[i] = pos_new.x;
        ions_out.pos_y[i] = pos_new.y;
        ions_out.pos_z[i] = pos_new.z;
        ions_out.vel_x[i] = vel_new.x;
        ions_out.vel_y[i] = vel_new.y;
        ions_out.vel_z[i] = vel_new.z;
    }
}

void integrate_verlet_batch(
    const IonStateGPU& ions_in,
    IonStateGPU& ions_out,
    Vec3 E_field, Vec3 B_field,
    double dt,
    cudaStream_t stream
) {
    int N = ions_in.count;
    int threads = 256;
    int blocks = std::min((N + threads - 1) / threads, 2048);
    
    integrate_verlet_batch_kernel<<<blocks, threads, 0, stream>>>(
        ions_in, ions_out, E_field, B_field, dt, N
    );
    
    CUDA_CHECK(cudaGetLastError());
}

} // namespace gpu
} // namespace icarion
```

#### Step 2: Extend GPUIntegrationHelper

**In `GPUIntegrationHelper.h`:**
```cpp
class GPUIntegrationHelper {
public:
    // ... existing methods ...
    
    /**
     * @brief Velocity Verlet integration (2nd-order symplectic)
     */
    bool integrate_batch_verlet(
        std::vector<IonState>& ions,
        double dt,
        double t,
        const IFieldProvider* field_provider = nullptr
    );
};
```

**In `GPUIntegrationHelper.cpp`:**
```cpp
bool GPUIntegrationHelper::integrate_batch_verlet(
    std::vector<IonState>& ions,
    double dt, double t,
    const IFieldProvider* field_provider
) {
    // Same pattern as RK4/Boris
    auto ions_in = pool_->get_device_buffer<IonStateGPU>(1);
    auto ions_out = pool_->get_device_buffer<IonStateGPU>(1);
    
    // Upload
    ion_state_conversion::upload_ions(ions, *ions_in, context_.get_stream());
    
    // Compute
    Vec3 E{0,0,0}, B{0,0,0};  // Or extract from field_provider
    integrate_verlet_batch(*ions_in, *ions_out, E, B, dt, context_.get_stream());
    
    // Download
    ion_state_conversion::download_ions(*ions_out, ions, context_.get_stream());
    context_.synchronize();
    
    stats_.gpu_integrations++;
    return true;
}
```

#### Step 3: Update SimulationEngine Dispatch

**In `SimulationEngine.h`:**
```cpp
#ifdef ICARION_USE_GPU
    enum class IntegratorType { RK4, RK45, Boris, Verlet, Unknown };
#endif
```

**In `SimulationEngine.cpp`:**
```cpp
// Cache integrator type
if (dynamic_cast<VerletStrategy*>(integrator_.get())) {
    integrator_type_ = IntegratorType::Verlet;
}

// Dispatch
case IntegratorType::Verlet:
    threshold /= 2;  // Verlet: 2 force evals, same as Boris
    return gpu_helper_->integrate_batch_verlet(ions, dt, t, field_provider);
```

#### Step 4: Validate CPU/GPU Parity

**Create `tests/integrator/test_verlet_parity.cpp`:**
```cpp
TEST_CASE("Verlet GPU/CPU Parity - Free Particle", "[gpu][verlet][parity]") {
    // Create identical ions
    std::vector<IonState> ions_cpu(1000), ions_gpu(1000);
    for (int i = 0; i < 1000; ++i) {
        ions_cpu[i] = ions_gpu[i] = create_test_ion();
    }
    
    // CPU reference
    VerletStrategy verlet_cpu;
    for (auto& ion : ions_cpu) {
        verlet_cpu.step(ion, 0.0, 1e-9, force_registry, ions_cpu);
    }
    
    // GPU
    auto context = GPUContext::create(0);
    auto helper = GPUIntegrationHelper::create(*context, 1);
    helper->integrate_batch_verlet(ions_gpu, 1e-9, 0.0);
    
    // Compare (tolerance: 1e-12)
    for (int i = 0; i < 1000; ++i) {
        REQUIRE_THAT(ions_gpu[i].pos.x, WithinAbs(ions_cpu[i].pos.x, 1e-12));
        REQUIRE_THAT(ions_gpu[i].vel.x, WithinAbs(ions_cpu[i].vel.x, 1e-12));
    }
}
```

#### Step 5: Performance Benchmarking

```cpp
// Measure threshold point where GPU becomes beneficial
for (int N : {1000, 2000, 3000, 5000, 10000}) {
    auto t_cpu = benchmark_cpu(N, verlet_cpu);
    auto t_gpu = benchmark_gpu(N, verlet_gpu);
    std::cout << N << " ions: " << (t_cpu / t_gpu) << "× speedup\n";
}
```

**Note:** GPU performance varies significantly with hardware, kernel complexity, and data transfer patterns. Benchmark your specific use case.

### GPU Integrator Design Guidelines

**Smart Threshold Selection:**
```cpp
// Rule of thumb: threshold ∝ 1 / (force_evaluations_per_step)
// Default base threshold: 5000 ions (configurable via gpu_threshold)
Boris:   1 eval  → threshold / 2
Verlet:  2 evals → threshold / 2
RK4:     4 evals → threshold
RK45:    6 evals → threshold
```

**Reuse Shared Code:**
```cpp
// All integrators use same force evaluation
#include "core/gpu/compute_acceleration.cuh"

__device__ Vec3 compute_acceleration(
    Vec3 pos, Vec3 vel, 
    double mass, double charge,
    Vec3 E_field, Vec3 B_field
) {
    // Lorentz force: F = q(E + v×B)
    Vec3 v_cross_B = {
        vel.y * B_field.z - vel.z * B_field.y,
        vel.z * B_field.x - vel.x * B_field.z,
        vel.x * B_field.y - vel.y * B_field.x
    };
    Vec3 F = E_field * charge + v_cross_B * charge;
    return F / mass;
}
```

**Memory Efficiency:**
```cpp
// Minimize register usage for high occupancy
// Target: <64 registers per thread
// Check with: nvcc --ptxas-options=-v
```

**Error Propagation:**
```cpp
// Multi-step integrators: accumulate error carefully
// Use Kahan summation for long simulations
__device__ double kahan_add(double sum, double x, double& c) {
    double y = x - c;
    double t = sum + y;
    c = (t - sum) - y;
    return t;
}
```

### Additional Resources

**CUDA Programming Guide:**
- https://docs.nvidia.com/cuda/cuda-c-programming-guide/

**CUDA Best Practices:**
- https://docs.nvidia.com/cuda/cuda-c-best-practices-guide/

**ICARION GPU Reference (v1.0):**
- GPU integration available for RK4, RK45, Boris integrators
- Automatic dispatch based on particle count (default threshold: 5000, Boris: 2500)
- Implementation details in `src/core/integrator/` and integration strategy classes
- Validation tests: `tests/integrator/test_rk45_boris_parity.cpp` (407 lines)
- Space charge GPU: `src/core/physics/spacecharge/SpaceChargeGPUModel.{h,cpp}` + `GPUSpaceChargeP3M.{h,cu}`

**GPU Performance:**
GPU performance varies significantly with hardware, simulation complexity, and data transfer patterns. Benchmark your specific use case to determine optimal threshold values.

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

### Common Testing Issues

#### Orbitrap Boundary Checking

**Problem:** Ions get deactivated immediately in Orbitrap simulations despite starting inside boundaries.

**Root Cause:** Hyperlogarithmic electrode surfaces require special boundary checking. Common mistakes:

1. **Using cylindrical radius check** instead of `find_domain_index()`:
   ```cpp
   // WRONG: Simple radius check fails for hyperlogarithmic geometry
   bool inside = (r <= domain.geometry.radius_m);
   
   // CORRECT: Use DomainManager's comprehensive check
   int domain_idx = domain_manager->find_domain_index(ion.pos);
   bool inside = (domain_idx >= 0);
   ```

2. **Wrong bisection bracket direction** for solving hyperlogarithmic equation:
   ```cpp
   // WRONG: Searches outward (diverges at large z)
   double r_lo = 0.3 * R;
   double r_hi = 3.0 * R;
   
   // CORRECT: Hyperlogarithmic surfaces curve INWARD
   double r_lo = 0.1 * R;
   double r_hi = R;
   ```

**Solution:** For Orbitrap instruments in `SimulationEngine`, use:
```cpp
if (domain_config.instrument == config::Instrument::Orbitrap) {
    Vec3 pos_global = domain_manager_->local_to_global_pos(pos_after, domain_idx);
    int check_domain = domain_manager_->find_domain_index(pos_global);
    still_inside = (check_domain == domain_idx);
}
```

**Affected Files:**
- `src/core/integrator/SimulationEngine.cpp` (boundary checks after integration)
- `src/core/integrator/DomainManager.cpp` (bisection solver for r(z))
- `tests/integrator/test_domain_manager.cpp` (test bisection implementation)

**See:** Commit `b4f358a` - "fix: Correct Orbitrap boundary checking"

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

**Note:** ICARION uses direct config references (const DomainConfig&) instead of parameter structs. This ensures Single Source of Truth (SSOT) compliance.

See: `src/core/physics/forces/ElectricFieldForce.h` (reference implementation)

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

1. **Reaction System Database Unification**:
   - Unify species types across modules
   - Wire reaction_handler directly into integrator
   - Delete legacy reaction loading code

2. **SimulationEngine Integration**:
   - Main simulation loop using integration strategies
   - Orchestration of ForceRegistry + CollisionHandler + ReactionHandler
   - Boundary condition handling
   - Output management (HDF5Writer v2)

### Planned

1. **InstrumentType Location** (Low priority, approximately 30min):
   - Move instrument/InstrumentTypes.h to core/config/types/InstrumentType.h
   - Eliminates dependency cycle config to instrument

2. **Stochastic Forces** (Design phase):
   - Separate random kicks from deterministic forces
   - New IStochasticForce interface (similar to collision handlers)

2. **Space Charge** [AVAILABLE]:
   - Automatic method selection: N<1000→Direct (O(N²)), N≥1000→Grid (O(N log N))
   - CIC charge deposition with O(h²) convergence
   - Poisson solver with 5 methods (Gauss-Seidel, SOR, CG, Multigrid, FFT)
   - Files: `src/core/physics/spacecharge/*`, `src/core/physics/forces/SpaceCharge{Direct,Grid}.{h,cpp}`
   - **When to use Direct:** N<1000, exact results needed
   - **When to use Grid:** N≥1000, fast approximation
   - **Configuration:** Set `physics.enable_space_charge = true` in config
   - **Performance tuning:** Adjust grid size (default 64³), update frequency
   - **Limitations:** CPU-only (GPU planned for future release), Dirichlet BC causes errors near boundaries

3. **GPU Space Charge Acceleration** (Planned for future release):
   - CUDA kernels for space charge Poisson solver
   - Grid-based field evaluation on GPU

4. **Field Caching**:
   - Pre-compute field on regular grid
   - Trilinear interpolation for fast evaluation
   - Benefit: Reduces repeated analytical field evaluations

### v1.0 Features

- **Force System SSOT**:
  - IForce interface with ForceRegistry
  - ElectricFieldForce, MagneticFieldForce, DampingForce, SpaceChargeForce
  - Full unit test coverage

- **Collision System SSOT**:
  - ICollisionHandler interface with factory
  - EHSS, HSS, OU collision handlers (SoA overrides implemented)
  - Energy conservation validation
  - Parity check: `tests/physics/collisions/test_collision_soa_parity.cpp`

- **Reaction System Handlers**:
  - IReactionHandler interface with factory
  - StochasticReactionHandler implementation (SoA path implemented)
  - GPUReactionHandler wrapper + `GPUReactionBackend` stub (factory returns the wrapper when `simulation.enable_gpu` is true; currently logs and falls back to CPU until kernels land)
  - Parity check: `tests/physics/reactions/test_reaction_soa_parity.cpp`
  - Database-driven reaction loading

- **Integration Strategies**:
  - IIntegrationStrategy interface
  - RK4Strategy (4th-order Runge-Kutta, fixed-step)
  - RK45Strategy (Dormand-Prince 5(4), adaptive timestep with FSAL)
  - BorisStrategy (symplectic pusher for electromagnetic fields)
  - IntegrationStrategyFactory for runtime selection
  - SSOT-compliant (uses DomainConfig pointer, callback-based acceleration)
  - Files: `src/core/integrator/strategies/*`

---

## Getting Help

- **Architecture Questions**: See `ARCHITECTURE.md`
- **API Reference**: See `docs/CLI_USAGE.md`
- **Config Schema**: See `schema/icarion-config.schema.json`
- **Bug Reports**: Create GitHub issue with minimal reproducer

---

**Document Status:** Living document, updated with each major refactoring phase.

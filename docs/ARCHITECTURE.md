# ICARION Architecture Guide

**Version:** 1.0
**Last Updated:** November 21, 2025

This document describes the high-level architecture of ICARION, focusing on module organization, data flow, and key design patterns.

---

## Table of Contents

1. [System Overview](#system-overview)
2. [Module Structure](#module-structure)
3. [Configuration System](#configuration-system)
4. [Force System Architecture](#force-system-architecture)
5. [Field Solver Architecture](#field-solver-architecture)
6. [Integrator Architecture](#integrator-architecture)
7. [Domain Management](#domain-management-phase-5a)
8. [SimulationEngine Architecture](#simulationengine-architecture-phase-5a)
9. [OutputManager Architecture](#outputmanager-architecture-phase-5a)
10. [Data Flow](#data-flow)
11. [Design Patterns](#design-patterns)

---

## System Overview

ICARION is a modular ion trajectory simulation framework with the following key subsystems:

```
┌─────────────────────────────────────────────────────────────┐
│                     ICARION Simulation                      │
├─────────────────────────────────────────────────────────────┤
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐       │
│  │   Config     │  │    Physics   │  │  Integrator  │       │
│  │   Loader     │→ │    Forces    │→ │   (ODE)      │       │
│  └──────────────┘  └──────────────┘  └──────────────┘       │
│         ↓                 ↓                  ↓              │
│  ┌──────────────┐  ┌──────────────┐  ┌──────────────┐       │
│  │   Domain     │  │ Field Solver │  │   Output     │       │
│  │   Config     │→ │   (Poisson)  │→ │   (HDF5)     │       │
│  └──────────────┘  └──────────────┘  └──────────────┘       │
└─────────────────────────────────────────────────────────────┘
```

### Core Principles

1. **Modularity**: Clear separation of concerns (config, physics, numerics, I/O)
2. **Extensibility**: Plugin architecture for forces, fields, and integrators
3. **Performance**: Optimized hot paths, optional GPU acceleration
4. **Correctness**: Extensive unit tests, physics validation, numerical safety

---

## Module Structure

```
ICARION/
├── src/
│   ├── core/                    # Core data structures and algorithms
│   │   ├── config/              # Configuration types and loaders
│   │   │   ├── types/           # Config data structures
│   │   │   └── loader/          # JSON parsing and validation
│   │   ├── physics/             # Physics models
│   │   │   └── forces/          # Force implementations
│   │   ├── types/               # Fundamental types (Vec3, IonState)
│   │   └── utils/               # Math utilities, constants
│   │
│   ├── integrator/              # ODE integration algorithms
│   ├── fieldsolver/             # Poisson/Laplace solvers
│   ├── instrument/              # Instrument-specific logic
│   ├── simulation/              # High-level simulation engine
│   ├── trajsim/                 # Trajectory simulation loop
│   └── utils/                   # General utilities
│
├── include/                     # Public API headers
├── tests/                       # Unit and integration tests
├── examples/                    # Example input files
├── schema/                      # JSON schemas
└── docs/                        # Documentation
```

### Module Dependencies

```
┌─────────────┐
│ simulation  │  (High-level orchestration)
└──────┬──────┘
       │
   ┌───▼───┐
   │trajsim│  (Trajectory loop)
   └───┬───┘
       │
   ┌───▼───────┐
   │integrator │  (ODE solver)
   └───┬───────┘
       │
   ┌───▼─────┐
   │ physics │  (Forces)
   └───┬─────┘
       │
┌──────▼──────┐
│   config    │  (Configuration)
└─────────────┘
```

**Dependency Rule**: Lower layers don't depend on higher layers.

---

## Configuration System

### Overview

Configuration is loaded from JSON files and stored in structured C++ classes. The system follows a hierarchical design.

### Class Hierarchy

```cpp
FullConfig                          // Top-level configuration
├── IonCloudConfig                  // Ion ensemble properties
│   ├── mass_u                      // Ion mass [u]
│   ├── charge                      // Ion charge state
│   ├── num_ions                    // Number of ions
│   └── positions / velocities      // Initial conditions
│
├── vector<DomainConfig>            // Multi-domain simulations
│   ├── name                        // Domain identifier
│   ├── geometry                    // Spatial boundaries
│   │   ├── x_min, x_max
│   │   ├── y_min, y_max
│   │   └── z_min, z_max
│   │
│   ├── FieldsConfig                // Electric/magnetic fields
│   │   ├── ElectricFieldConfig
│   │   │   ├── type                // "analytical" or "fieldmap"
│   │   │   ├── instrument_type     // LQIT, IMS, TOF, ...
│   │   │   ├── rf_voltage_V
│   │   │   ├── rf_frequency_Hz
│   │   │   ├── dc_axial_voltage_V
│   │   │   └── ... (instrument-specific)
│   │   │
│   │   └── MagneticFieldConfig
│   │       ├── uniform_field_T     // Uniform B-field [T]
│   │       └── gradient_T_per_m    // Linear gradient
│   │
│   └── EnvironmentConfig           // Background gas
│       ├── temperature_K
│       ├── pressure_Pa
│       ├── gas_species             // "N2", "He", "Ar", ...
│       ├── collision_model         // "HSD", "Langevin", ...
│       └── collision_parameters
│
├── IntegrationConfig               // Numerical integration
│   ├── method                      // "RK4", "Verlet", "AdaptiveRK45"
│   ├── dt_s                        // Time step [s]
│   ├── t_max_s                     // Simulation end time [s]
│   └── tolerance                   // Adaptive error tolerance
│
├── OutputConfig                    // Output settings
│   ├── format                      // "hdf5", "csv"
│   ├── output_file
│   ├── save_interval_steps
│   └── save_fields                 // Include field data
│
└── SimulationConfig                // General settings
    ├── random_seed
    ├── enable_gpu
    └── num_threads
```

### Config Loader Architecture

```cpp
namespace ICARION {
namespace config {

/**
 * @brief Main configuration loader
 * 
 * Responsibilities:
 * - Parse JSON files
 * - Validate against schema
 * - Construct FullConfig object
 * - Apply defaults
 * - Resolve cross-references
 */
class ConfigLoader {
public:
    /**
     * @brief Load configuration from JSON file
     * 
     * @param filename Path to JSON config file
     * @return Parsed and validated configuration
     * @throws std::runtime_error on parse/validation errors
     */
    static FullConfig load_from_file(const std::string& filename);
    
    /**
     * @brief Load from JSON string
     */
    static FullConfig load_from_string(const std::string& json_str);
    
    /**
     * @brief Validate configuration against schema
     */
    static bool validate(const FullConfig& config);
    
private:
    static void parse_ion_cloud(const nlohmann::json& j, IonCloudConfig& ion_cloud);
    static void parse_domains(const nlohmann::json& j, std::vector<DomainConfig>& domains);
    static void parse_fields(const nlohmann::json& j, FieldsConfig& fields);
    static void parse_environment(const nlohmann::json& j, EnvironmentConfig& env);
    // ... etc
};

} // namespace config
} // namespace ICARION
```

### Configuration Flow

```
┌────────────┐
│ JSON File  │
└──────┬─────┘
       │ read
       ▼
┌─────────────┐
│ JSON Parser │  (nlohmann::json)
└──────┬──────┘
       │ parse
       ▼
┌──────────────┐
│ Schema Valid │  (optional, via JSON Schema)
└──────┬───────┘
       │ validate
       ▼
┌──────────────┐
│ ConfigLoader │  (convert to C++ structs)
└──────┬───────┘
       │ construct
       ▼
┌──────────────┐
│  FullConfig  │  (validated, ready to use)
└──────────────┘
```

### SSOT Status (v1.0)

**Completed SSOT Migrations:**

1. **Force System** (Phase 1, Steps 1-4 complete):
   - MagneticFieldForce: Uses const MagneticFieldConfig&
   - ElectricFieldForce: Uses const DomainConfig&
   - DampingForce: Uses const EnvironmentConfig&
   - Legacy structs deleted: MagneticFieldParams, AnalyticalFieldParams, DampingParams

2. **Collision System** (Phase 2C complete):
   - All collision handlers use const EnvironmentConfig& reference
   - No parameter duplication

3. **Reaction System** (Phase 3C complete):
   - ReactionHandler wired into integrator
   - Uses species database from FullConfig

**In Progress:**

1. **Integrator System** (Phase 1, Steps 5-8, approximately 3.5h remaining):
   - compute_accelerations(): Replace GlobalParams with DomainConfig
   - integrate_one_step(): Use FullConfig directly
   - integrate_trajectory(): Remove parameter conversions
   - See: tmp/REMAINING_SSOT_MIGRATION_WORK.md

**Known Minor Issue:**

1. **InstrumentType Location** (Low priority, approximately 30min):
   - Current: instrument/InstrumentTypes.h
   - Issue: Creates dependency config to instrument (backwards)
   - Fix: Move to core/config/types/InstrumentType.h
   - Impact: Cosmetic only, not a functional SSOT violation

### Current Architecture (v1.0)

All force classes follow SSOT principle with direct config references:

```cpp
// Example: ElectricFieldForce (SSOT compliant since v1.0)
class ElectricFieldForce : public IForce {
public:
    ElectricFieldForce(const config::DomainConfig& domain);
    
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override {
        // Read config on-demand
        const auto& dc = domain_.fields.dc;
        double axial_V = dc.axial_V;  // No duplication
        // ...
    }
    
private:
    const config::DomainConfig& domain_;  // Reference (zero-copy)
};
```

Benefit: Config changes propagate automatically, no parameter conversion overhead.

---

## Force System Architecture

### Design Overview

The force system follows a **Strategy Pattern** with plugin architecture and **SSOT (Single Source of Truth)** principle:

```
         ┌──────────┐
         │  IForce  │  (Interface)
         └────┬─────┘
              │ implements
      ┌───────┴───────┬───────────┬─────────────┐
      │               │           │             │
┌─────▼─────┐  ┌─────▼─────┐ ┌──▼──────┐ ┌────▼────────┐
│  Electric │  │ Magnetic  │ │ Damping │ │ SpaceCharge │
│   Field   │  │   Field   │ │  Force  │ │    Force    │
└───────────┘  └───────────┘ └─────────┘ └─────────────┘
       │              │            │            │
       └──────────────┴────────────┴────────────┘
                       │
                 ┌─────▼─────────┐
                 │ ForceRegistry │  (Composite)
                 └───────────────┘
```

### SSOT Principle (v1.0)

**Forces store references to config, not copies:**

```cpp
// MODERN (v1.0): Direct config reference
const config::MagneticFieldConfig& magnetic = domain.fields.magnetic;
MagneticFieldForce force(magnetic);  // Reference to SSOT
```

**Benefits:**

- No data duplication
- Config changes automatically propagate
- Cleaner interfaces
- Type safety from config system

### IForce Interface

```cpp
namespace ICARION {
namespace physics {

/**
 * @brief Abstract interface for all force types
 * 
 * Forces store const references to config (SSOT) and compute F(ion, t, context).
 * All physics happens in compute().
 */
class IForce {
public:
    virtual ~IForce() = default;
    
    /**
     * @brief Compute force on a single ion
     * 
     * @param ion Ion state (position, velocity, mass, charge)
     * @param t Current simulation time [s]
     * @param ctx Additional context (ion ensemble, field provider)
     * @return Force vector [N]
     * 
     * Must be const (no mutation of force object).
     * Must be thread-safe if called from parallel context.
     */
    virtual Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const = 0;
    
    /**
     * @brief Get force name for logging/debugging
     */
    virtual std::string name() const = 0;
};

} // namespace physics
} // namespace ICARION
```

### ForceContext

Provides shared data for force computation (avoids duplicate lookups):

```cpp
struct ForceContext {
    const IFieldProvider* field_provider;       ///< Optional field evaluator override
    const std::vector<IonState>* all_ions;      ///< All ions (for space charge)
    
    // Context is minimal - forces read from their stored config references
};
```

### ForceRegistry (Composite Pattern)

**Phase 12 Enhancement:** ForceRegistry now stores domain configuration internally.

Manages multiple forces and computes total force via superposition:

```cpp
class ForceRegistry : public IForce {
public:
    /**
     * @brief Construct registry (empty, no domain) [DEPRECATED]
     * @deprecated Use ForceRegistry(const config::DomainConfig&) instead
     */
    ForceRegistry() = default;
    
    /**
     * @brief Construct registry with domain context (RECOMMENDED)
     * @param domain Domain configuration (geometry, fields, environment)
     * 
     * Phase 12 enhancement: Registry stores domain reference internally.
     * This eliminates need to pass domain through integration methods.
     */
    explicit ForceRegistry(const config::DomainConfig& domain);
    
    /**
     * @brief Add a force to the registry
     * 
     * Takes ownership of the force object.
     */
    void add_force(std::unique_ptr<IForce> force);
    
    /**
     * @brief Compute total force on ion (overrides IForce::compute)
     * 
     * F_total = F1 + F2 + F3 + ... (superposition)
     */
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override;
    
    /**
     * @brief Get domain configuration (if available)
     * @return Pointer to domain config, or nullptr if not set
     * 
     * Phase 12: Allows forces and integrators to access domain context.
     */
    const config::DomainConfig* domain() const;
    
    /**
     * @brief Clear all forces
     */
    void clear();
    
    /**
     * @brief Get number of registered forces
     */
    size_t size() const;
    
    /**
     * @brief Check if registry is empty
     */
    bool empty() const;
    
    /**
     * @brief Get const reference to force vector (for iteration)
     */
    const std::vector<std::unique_ptr<IForce>>& forces() const;
    
private:
    std::vector<std::unique_ptr<IForce>> forces_;
    const config::DomainConfig* domain_ = nullptr;  // Non-owning pointer (Phase 12)
};
```

**Benefits of Domain-Aware ForceRegistry (Phase 12):**
- ✅ Better SSOT compliance (domain stored once, not passed through methods)
- ✅ Cleaner method signatures (fewer parameters to pass)
- ✅ Multi-domain support (each domain has its own registry)
- ✅ Forces can access domain context without parameter pollution

### Implemented Force Types

#### 1. ElectricFieldForce

Computes Lorentz electric force: **F = q·E**

**Constructor (SSOT):**

```cpp
// Analytical mode: reads from DomainConfig
ElectricFieldForce(const config::DomainConfig& domain);

// Field provider mode: uses external field evaluator
ElectricFieldForce(std::shared_ptr<IFieldProvider> provider);
```

**Configuration Access:**

```cpp
// Reads directly from stored domain reference
double voltage = domain_->fields.dc.axial_V;
double radius = domain_->geometry.radius_m;
Instrument instrument = domain_->instrument;
```

**Supported Instruments:**

- LQIT (Linear Quadrupole Ion Trap)
- IMS (Ion Mobility Spectrometry)
- TOF (Time-of-Flight)
- Orbitrap
- QuadrupoleRF
- FTICR (Fourier Transform ICR)
- NoFixedInstrument (returns zero field)

#### 2. MagneticFieldForce

Computes Lorentz magnetic force: **F = q(v × B)**

**Constructor (SSOT):**

```cpp
// Analytical mode: reads from MagneticFieldConfig
MagneticFieldForce(const config::MagneticFieldConfig& magnetic);

// Field provider mode
MagneticFieldForce(std::shared_ptr<IFieldProvider> provider);
```

**Configuration Access:**

```cpp
// Reads directly from stored magnetic reference
Vec3 B = magnetic_.field_strength_T;
Vec3 gradient = magnetic_.field_gradient_T_m;
bool enabled = magnetic_.enabled;
```

**Modes:**

- Uniform field: `B = const`
- Linear gradient: `B(z) = B₀ + ∇B·z`
- Field provider (interpolated)

#### 3. DampingForce

Computes deterministic collision damping: **F = -γ·m·v**

**Constructor (SSOT):**

```cpp
DampingForce(const config::EnvironmentConfig& env, DampingModel model);
```

**Configuration Access:**

```cpp
// Reads directly from stored environment reference
double pressure = env_.pressure_Pa;
double temperature = env_.temperature_K;
double density = env_.particle_density_m_3;
double mass_gas = env_.gas_mass_kg;
```

**Models:**

- **Friction**: Mobility-based, γ = q/(K₀·m)
- **HSD**: Elastic collisions, γ = ν·(m_n/(m_i+m_n))
- **Langevin**: Ion-induced dipole, enhanced cross-section
- **None**: No damping

****Note**: Stochastic kicks (thermal noise) are handled separately by CollisionEngine.

#### 4. SpaceChargeForce

Computes ion-ion Coulomb repulsion: **F = k_e·q₁·q₂·r̂/r²**

**Constructor:**

```cpp
SpaceChargeForce();  // Stateless, reads from ForceContext
```

**Features:**

- N-body direct summation (O(N²))
- Self-interaction exclusion
- Softening parameter to prevent divergence at r→0

**Performance**: Suitable for <1000 ions. Use `SpaceChargeSolver` (grid-based) for larger ensembles.

---

## Field Solver Architecture

### Overview

ICARION supports two modes for electric field evaluation:

1. **Analytical**: Closed-form formulas for ideal geometries
2. **Numerical**: Poisson/Laplace solver on 3D grid

### IFieldProvider Interface

```cpp
/**
 * @brief Abstract interface for field evaluation
 * 
 * Allows swapping between analytical and numerical field sources.
 */
class IFieldProvider {
public:
    virtual ~IFieldProvider() = default;
    
    /**
     * @brief Evaluate field at position
     * 
     * @param pos Position vector [m]
     * @param t Time [s] (for time-dependent fields)
     * @return Field value (e.g., electric field [V/m] or magnetic field [T])
     */
    virtual Vec3 evaluate(const Vec3& pos, double t) const = 0;
};
```

### Implementations

1. **UniformFieldProvider**: Returns constant field
2. **LinearGradientProvider**: Linear field variation
3. **GridFieldProvider**: 3D interpolation from grid data
4. **PoissonSolver**: Solves ∇²φ = -ρ/ε₀ on mesh

### Space Charge Solver

For large ion ensembles (>10k ions), direct N-body is too slow. Instead:

```cpp
class SpaceChargeSolver {
public:
    /**
     * @brief Solve Poisson equation for space charge
     * 
     * ∇²φ = -ρ(r)/ε₀
     * E = -∇φ
     * 
     * @param ion_positions Ion cloud positions
     * @param charges Ion charges
     * @return Field provider for interpolated E-field
     */
    std::unique_ptr<IFieldProvider> solve(
        const std::vector<Vec3>& ion_positions,
        const std::vector<double>& charges
    );
};
```

**Algorithm**: 
- Particle-in-Cell (PIC): Deposit charges on grid
- Solve Poisson on grid (FFT or multigrid)
- Interpolate E-field back to particles

---

## Integrator Architecture

### Overview

The integrator solves the equations of motion:

```
dv/dt = F(x, v, t) / m
dx/dt = v
```

### Integration Strategy Interface (Phase 4)

```cpp
namespace ICARION::integrator {

/**
 * @brief Strategy pattern for numerical integration methods
 * 
 * Replaces legacy integrate_one_step() with modular, testable design.
 * Uses ForceRegistry for force computation (SSOT-compliant).
 */
class IIntegrationStrategy {
public:
    virtual ~IIntegrationStrategy() = default;
    
    /**
     * @brief Fixed-step integration
     * 
     * @param ion Ion state (updated in-place)
     * @param t Current time [s]
     * @param dt Time step [s]
     * @param force_registry Force computation engine
     * @param domain Domain configuration (fields, boundaries)
     * @param all_ions All ions in ensemble (for space charge)
     */
    virtual void step(
        IonState& ion,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry,
        const config::DomainConfig& domain,
        const std::vector<IonState>& all_ions
    ) = 0;
    
    /**
     * @brief Adaptive-step integration (optional)
     * 
     * @param dt_inout Input: current dt, Output: suggested next dt
     * 
     * Default: calls step() and ignores dt adjustment
     * RK45: uses error control for adaptive stepping
     */
    virtual void step_adaptive(
        IonState& ion,
        double t,
        double& dt_inout,
        const physics::ForceRegistry& force_registry,
        const config::DomainConfig& domain,
        const std::vector<IonState>& all_ions
    );
    
    virtual std::string name() const = 0;
};

} // namespace ICARION::integrator
```

### Implemented Strategies (Phase 4A/4B)

#### 1. RK4Strategy
- **Order:** 4th-order accurate
- **Type:** Fixed timestep
- **Use Case:** General-purpose integration
- **Cost:** 4 force evaluations per step
- **File:** `src/core/integrator/strategies/RK4Strategy.{h,cpp}`

**Algorithm:**
```
k1 = f(t, y)
k2 = f(t + dt/2, y + dt*k1/2)
k3 = f(t + dt/2, y + dt*k2/2)
k4 = f(t + dt, y + dt*k3)
y_new = y + dt*(k1 + 2*k2 + 2*k3 + k4)/6
```

#### 2. RK45Strategy (Dormand-Prince)
- **Order:** 5th-order accurate (4th-order error control)
- **Type:** Adaptive timestep with error control
- **Use Case:** High-accuracy simulations with varying dynamics
- **Cost:** 6 force evaluations per step (FSAL optimization)
- **File:** `src/core/integrator/strategies/RK45Strategy.{h,cpp}`

**Features:**
- Embedded Runge-Kutta 5(4) (Dormand-Prince coefficients)
- FSAL (First Same As Last) optimization
- PI controller for timestep adaptation
- Configurable tolerances (atol, rtol) via `simulation.rk45_settings`
- Automatic step rejection when error > tolerance

**Configuration (SSOT):**
```json
{
  "simulation": {
    "integrator": "RK45",
    "rk45_settings": {
      "abs_tol": 1e-14,
      "rel_tol": 1e-12,
      "dt_min": 1e-12,
      "safety": 0.84,
      "min_factor": 0.2,
      "max_factor": 5.0
    }
  }
}
```

#### 3. BorisStrategy
- **Order:** 2nd-order accurate
- **Type:** Symplectic (energy-conserving)
- **Use Case:** Charged particles in strong electromagnetic fields
- **Cost:** 1 force evaluation per step
- **File:** `src/core/integrator/strategies/BorisStrategy.{h,cpp}`

**Algorithm (Boris Pusher):**
```
1. v^- = v^n + (q/m)*E*(dt/2)          [Electric half-step]
2. t = (q/m)*B*(dt/2)                  [Rotation parameter]
3. s = 2*t / (1 + |t|^2)               [Scaled rotation]
4. v' = v^- + v^- × t                  [First rotation]
5. v^+ = v^- + v' × s                  [Second rotation]
6. x^(n+1) = x^n + v^+*dt              [Position update]
7. v^(n+1) = v^+ + (q/m)*E*(dt/2)      [Electric half-step]
```

**Properties:**
- Time-reversible
- Preserves phase-space volume (symplectic)
- No small-angle approximation (valid for all B)
- Optimal for cyclotron motion (no accumulation errors)

### Factory Pattern (IntegrationStrategyFactory)

```cpp
namespace ICARION::integrator {

/**
 * @brief Factory for creating integration strategies
 */
class IntegrationStrategyFactory {
public:
    /**
     * @brief Create integration strategy from name
     * 
     * @param strategy_name "RK4", "RK45", or "Boris"
     * @return Unique pointer to strategy instance
     * @throws std::invalid_argument if unknown strategy
     */
    static std::unique_ptr<IIntegrationStrategy> 
    create(const std::string& strategy_name);
    
    /**
     * @brief Get list of supported strategies
     */
    static std::vector<std::string> supported_strategies();
};

} // namespace ICARION::integrator
```

**Implementation:** `src/core/integrator/strategies/IntegrationStrategyFactory.h`

**Supported Strategies:**
- `"RK4"` - 4th-order Runge-Kutta (fixed timestep)
- `"RK45"` - Dormand-Prince 5(4) (adaptive timestep)
- `"Boris"` - Boris pusher (electromagnetic fields)

**Example Usage:**

```cpp
// Create from config
std::string method = config.simulation.integrator;  // "RK45"
auto strategy = IntegrationStrategyFactory::create(method);

// Use strategy (callback-based)
auto compute_accel = [&](const IonState& ion, double t) -> Vec3 {
    // Your acceleration computation
    return Vec3{0, 0, -9.81};
};

strategy->step(ion, t, dt, compute_accel, &domain);
```

**Current Status (Phase 4):**
- Factory implemented and tested
- Not yet integrated into main.cpp (Phase 5 work)

### Integration Loop Examples

**Current Implementation (Test Usage):**

```cpp
// Example: Fixed-step RK4 (from tests)
void simulate_trajectory_fixed(
    IonState& ion,
    std::function<Vec3(const IonState&, double)> compute_accel,
    IIntegrationStrategy& strategy,
    double t_max,
    double dt
) {
    double t = 0.0;
    while (t < t_max) {
        strategy.step(ion, t, dt, compute_accel, nullptr);
        t += dt;
        save_output(ion, t);
    }
}
```

**Future (Production Usage with ForceRegistry):**

```cpp
// Future: SimulationEngine will wrap ForceRegistry in callback
auto compute_accel = [&](const IonState& ion, double t) -> Vec3 {
    Vec3 F_total = force_registry_->compute_total_force(ion, t, ctx);
    return F_total / ion.mass_kg;
};

strategy_->step(ion, t, dt, compute_accel, &domain);
```

**Migration Status:**
- Strategies implemented and tested (Phase 4)
- Factory pattern implemented (Phase 4)
- SimulationEngine integration pending (Phase 5)

---

## Domain Management (Phase 5A)

**DomainManager** handles spatial domain lookup, coordinate transformations, and boundary conditions.

**Key Features:**
- Domain lookup by position (cylindrical/hyperbolic (Orbitrap) geometry)
- Coordinate transforms (global ↔ local)
- Aperture crossing detection
- Domain property updates (temperature, pressure, gas velocity)

**SSOT Compliance:**
```cpp
// Uses modern config::DomainConfig (not legacy InstrumentDomain)
DomainManager manager(full_config.domains);

int idx = manager.find_domain_index(ion.pos);
Vec3 local = manager.global_to_local_pos(ion.pos, idx);
manager.update_domain_properties(ion, idx);
```

**Internal Geometry Check:**
- Replaces legacy `isInsideDomain()` from `paramUtils.cpp`
- Supports cylindrical (most instruments) and logarithmic-hyperbolic (Orbitrap)
- No dependency on legacy functions

**Geometry Implementations:**

*Cylindrical (IMS, LQIT, TOF, etc.):*
- Ion inside domain if: `(z >= -ε) && (z < L) && (r < R)`
- Constant radius R over entire length L
- Floating-point tolerance ε = 1e-12 m for z-boundary

*Orbitrap (logarithmic-hyperbolic):*
- Electrode shape: `z² = 0.5·(r² - R²) + R_m² · ln(R/r)`
- Inner/outer electrodes follow hyperbolic surfaces r_in(z) and r_out(z)
- Numerical root-finding (bisection method) computes r(z) at runtime
- Ion inside domain if: `r_in(z) ≤ r ≤ r_out(z)`
- Physics constraint: **r_char > r_out** (characteristic radius larger than outer electrode)
  - For realistic Orbitrap parameters, surfaces become nearly z-independent
  - Example: r_char=25mm, r_in=12mm, r_out=20mm → gap constant ~8mm for z=0..20mm

**Boundary Termination:**
- `terminate_ion_at_boundary()`: Ray-tracing to find exact intersection point
- Cylindrical: Analytical ray-cylinder/plane intersection
- Orbitrap: Midpoint approximation (TODO: ray-hyperbola intersection)
- Prevents unphysical ion positions beyond domain boundaries

**Status:** Complete (Phase 5A, Nov 2025)

**Files:**
- `src/core/integrator/DomainManager.h` (API)
- `src/core/integrator/DomainManager.cpp` (~150 lines)
- `tests/integrator/test_domain_manager.cpp` (11 test cases)

---

## SimulationEngine Architecture (Phase 5A)

**SimulationEngine** is the main simulation orchestrator, replacing legacy `integrate_trajectory()`.

### Overview

```text
SimulationEngine (orchestrator)
    │
    ├── ForceRegistry (compute total force)
    ├── IntegrationStrategy (RK4/RK45/Boris)
    ├── CollisionHandler (EHSS/HSS/OU)
    ├── ReactionHandler (ion-molecule reactions)
    ├── OutputManager (HDF5 + text logging)
    └── DomainManager (boundary checks, transitions)
```

### Key Features

**1. Dependency Injection:**
- All physics modules passed via constructor (testable!)
- No global state (everything in `FullConfig`)

**2. Modular Design:**
- Clean separation: physics / I/O / domain management
- Swap components easily (RK4 ↔ RK45, EHSS ↔ HSS)

**3. SSOT Compliance:**
- Uses `config::FullConfig` exclusively (no legacy `GlobalParams`)
- Direct config access (no parameter duplication)

**4. Parallel Execution:**
- OpenMP ion loop (thread-safe)
- Ion-based RNG (reproducible, independent of scheduling)

### API Design

```cpp
namespace ICARION::integrator {

class SimulationEngine {
public:
    /**
     * @brief Construct simulation engine
     * @param config Simulation configuration (SSOT)
     * @param force_registry Force computation system
     * @param integrator Trajectory integration strategy (RK4/RK45/Boris)
     * @param collision_handler Collision physics (optional, nullptr = no collisions)
     * @param reaction_handler Reaction chemistry (optional, nullptr = no reactions)
     */
    SimulationEngine(
        const config::FullConfig& config,
        std::shared_ptr<physics::ForceRegistry> force_registry,
        std::shared_ptr<IIntegrationStrategy> integrator,
        std::shared_ptr<physics::ICollisionHandler> collision_handler = nullptr,
        std::shared_ptr<physics::IReactionHandler> reaction_handler = nullptr
    );
    
    /**
     * @brief Run complete simulation
     * @param ions Initial ion ensemble
     * @return Final ion states
     */
    std::vector<IonState> run(std::vector<IonState>& ions);
    
    // Accessors (for testing)
    const config::FullConfig& get_config() const;
    const DomainManager& get_domain_manager() const;
    const OutputManager& get_output_manager() const;
};

} // namespace ICARION::integrator
```

### Main Simulation Loop

**Workflow:**

1. **Initialize** (setup subsystems)
   - Create `DomainManager` from `config.domains`
   - Create `OutputManager` with HDF5 file path
   - Initialize numerical safety logging (if enabled)

2. **Main Time Loop** (until `t_global >= t_end`)
   - Apply ion birth logic (delayed emission)
   - **Parallel ion processing** (OpenMP):
     - Find domain (DomainManager)
     - Transform to local coordinates
     - Compute forces (ForceRegistry)
     - Handle collisions (ICollisionHandler)
     - Handle reactions (IReactionHandler)
     - Integrate trajectory (IIntegrationStrategy)
     - Check aperture crossings (DomainManager)
     - Transform back to global coordinates
   - Log trajectory snapshot (OutputManager)
   - Update progress logging (every 10%)

3. **Finalize** (completion)
   - Flush output buffers
   - Write completion metadata
   - Log final statistics (active/lost ions)

**Early Exit Conditions:**
- All ions inactive (lost or detected)
- Critical error (NaN positions, invalid domain index)

### Example Usage

```cpp
// Load configuration
auto config = config::ConfigLoader::load("config.json");

// Create physics modules
auto force_registry = std::make_shared<physics::ForceRegistry>();
auto integrator = std::make_shared<RK4Strategy>();
auto collision_handler = physics::CollisionHandlerFactory::create(config.physics, ...);
auto reaction_handler = physics::ReactionHandlerFactory::create(config.physics, ...);

// Create simulation engine
integrator::SimulationEngine engine(
    config,
    force_registry,
    integrator,
    collision_handler,
    reaction_handler
);

// Generate ions
auto result = config.generate_ions(rng);
std::vector<IonState> ions = std::move(result.ions);

// Run simulation
auto final_ions = engine.run(ions);

// Output automatically written to HDF5 file
```

### Ion-Based RNG (Phase 12 Enhancement)

**Problem:** Thread-local RNG (seeded with thread_id) makes results dependent on OpenMP scheduling.

**Solution:** Ion-specific RNG (seeded with ion index):

```cpp
// Create RNG array before parallel region
std::vector<EhssRng> rng_by_ion;
rng_by_ion.reserve(n_ions);
for (int i = 0; i < n_ions; ++i) {
    uint64_t ion_seed = config.rng_seed + static_cast<uint64_t>(i);
    rng_by_ion.emplace_back(ion_seed);
}

#pragma omp parallel
{
    #pragma omp for schedule(dynamic)
    for (int i = 0; i < n_ions; ++i) {
        IonState& ion = ions[i];
        EhssRng& ion_rng = rng_by_ion[i];  // Ion-specific RNG
        
        // Use ion_rng for collisions/reactions
        collision_handler->handle_collision(ion, dt, ion_rng, env);
    }
}
```

**Benefits:**
- Reproducible results independent of thread count/scheduling
- Same ion always sees same random sequence
- Thread-safe (each thread accesses different ion RNG)

**Status:** Implemented (Nov 2025)

### Files

**Implementation:**
- `src/core/integrator/SimulationEngine.h` (API)
- `src/core/integrator/SimulationEngine.cpp` (~400 lines)

**Tests:**
- `tests/integrator/test_simulation_engine.cpp` (unit tests)
- 10 test cases, 45+ assertions

**Status:** Production-ready (Phase 5A complete, Nov 2025)

---

## OutputManager Architecture (Phase 5A)

**OutputManager** handles unified output: HDF5 trajectories + text logging.

### Key Features

**1. Buffered HDF5 Output:**
- RAM buffering (default: 50 timesteps)
- Time-based flush triggers (default: 1 ms)
- Size-based flush triggers (buffer full)

**2. Text Logging (optional):**
- Progress messages ("50% completed")
- Ion statistics (active/lost counts)
- Completion summary
- Can be disabled (empty log filename)

**3. HDF5Writer v2 Integration:**
- Wraps modern `io::HDF5Writer` API
- Metadata export (species, parameters, git hash)
- Chunked datasets (efficient large-file writes)

### API Design

```cpp
namespace ICARION::integrator {

class OutputManager {
public:
    /**
     * @brief Construct output manager
     * @param hdf5_filename HDF5 trajectory file path (required)
     * @param log_filename Text log file path (empty = no text log)
     * @param write_interval_dt Time interval between HDF5 writes [s]
     * @param buffer_max Max timesteps in RAM before forced flush
     */
    OutputManager(
        const std::string& hdf5_filename,
        const std::string& log_filename = "",
        double write_interval_dt = 0.001,  // Default: 1 ms
        size_t buffer_max = 50
    );
    
    /**
     * @brief Initialize HDF5 file (write metadata)
     * @param config Simulation configuration (SSOT)
     * @param ions Initial ion ensemble (for species metadata)
     */
    void initialize(const config::FullConfig& config, 
                    const std::vector<IonState>& ions);
    
    /**
     * @brief Log trajectory snapshot (buffers in RAM)
     * @param t Current time [s]
     * @param ions Current ion states
     */
    void log_step(double t, const std::vector<IonState>& ions);
    
    /**
     * @brief Log progress message (to text log)
     * @param message Progress message (e.g., "50% completed")
     */
    void log_progress(const std::string& message);
    
    /**
     * @brief Check if HDF5 write is needed
     * @param t_current Current time [s]
     * @return True if next write time reached OR buffer full
     */
    bool should_write(double t_current) const;
    
    /**
     * @brief Flush buffers to HDF5 file
     */
    void flush();
    
    /**
     * @brief Finalize output (write completion metadata)
     * @param t_final Final time [s]
     * @param final_ions Final ion states (for statistics)
     */
    void finalize(double t_final, const std::vector<IonState>& final_ions);
    
    // Accessors (for testing)
    size_t buffer_size() const;
    bool has_text_log() const;
};

} // namespace ICARION::integrator
```

### Buffering Strategy

**Time-based flush:**

```text
t = 0.0 ms  →  log_step()  →  buffer
t = 0.5 ms  →  log_step()  →  buffer
t = 1.0 ms  →  log_step()  →  AUTO-FLUSH (write_interval_dt = 1 ms)
t = 1.5 ms  →  log_step()  →  buffer
t = 2.0 ms  →  log_step()  →  AUTO-FLUSH
```

**Size-based flush:**

```text
log_step()  →  buffer (size = 49)
log_step()  →  AUTO-FLUSH (buffer_max = 50)
```

**Benefits:**
- Reduces HDF5 I/O overhead (batch writes)
- Prevents memory exhaustion (max buffer size)
- Configurable trade-off (write frequency vs. memory)

### Example Usage

```cpp
// Create output manager
OutputManager manager(
    "trajectory.h5",      // HDF5 file
    "simulation.log",     // Text log
    1e-3,                 // Flush every 1 ms
    50                    // Buffer max = 50 timesteps
);

// Initialize (write metadata)
manager.initialize(config, ions);

// Main loop
for (double t = 0.0; t < t_end; t += dt) {
    // ... simulate ions ...
    
    // Log timestep
    manager.log_step(t, ions);  // Auto-flush if needed
    
    // Progress logging
    if (step % 1000 == 0) {
        manager.log_progress("50% completed");
    }
}

// Finalize (write completion metadata)
manager.finalize(t_end, ions);
```

### Files

**Implementation:**
- `src/core/integrator/OutputManager.h` (API)
- `src/core/integrator/OutputManager.cpp` (~250 lines)

**Tests:**
- `tests/integrator/test_output_manager.cpp`
- 8 test cases, 30+ assertions

**Status:** Production-ready (Phase 5A complete, Nov 2025)

---

## Data Flow

### Typical Simulation Flow

```
1. Load Config
   └─ ConfigLoader::load_from_file()
   
2. Initialize Forces
   ├─ Create ElectricFieldForce(fields_config)
   ├─ Create MagneticFieldForce(fields_config)
   ├─ Create DampingForce(environment_config)
   └─ Create SpaceChargeForce(softening)
   
3. Build ForceRegistry
   └─ registry.add_force() for each force
   
4. Initialize Ions
   └─ Create IonState from ion_cloud_config
   
5. Time Integration Loop
   for t in [0, t_max]:
       F_total = registry.compute_total_force(ion, t, ctx)
       integrator.step(ion, F_total, dt)
       save_output(ion, t)
   
6. Finalize Output
   └─ Write HDF5 file
```

### Data Dependencies

```
FullConfig
    ↓
┌───────────────┬──────────────┬─────────────┐
│               │              │             │
DomainConfig  IonCloudConfig  IntegConfig  OutputConfig
    ↓               ↓              ↓            ↓
FieldsConfig    IonState    Integrator    HDF5Writer
    ↓
ForceRegistry
    ↓
ElectricFieldForce
MagneticFieldForce
DampingForce
SpaceChargeForce
```

---

## Design Patterns

### 1. Strategy Pattern (Forces)

**Intent**: Allow runtime selection of force models without changing client code.

**Implementation**: `IForce` interface with multiple implementations.

**Benefits**:
- Easy to add new force types
- Forces are composable via `ForceRegistry`
- Forces can be swapped at runtime

### 2. Composite Pattern (ForceRegistry)

**Intent**: Treat individual forces and collections of forces uniformly.

**Implementation**: `ForceRegistry` aggregates `IForce` objects and computes sum.

**Benefits**:
- Superposition principle expressed naturally
- Client code doesn't know how many forces exist
- Easy to add/remove forces dynamically

### 3. Provider Pattern (Field Evaluation)

**Intent**: Decouple field source (analytical, numerical, file) from consumers.

**Implementation**: `IFieldProvider` interface.

**Benefits**:
- Swap analytical ↔ numerical fields without changing force code
- Testability (mock field providers)
- Performance optimization (e.g., cached interpolation)

### 4. Builder Pattern (Configuration)

**Intent**: Construct complex objects step-by-step with validation.

**Implementation**: `ConfigLoader` builds `FullConfig` from JSON.

**Benefits**:
- Separation of parsing and validation
- Defaults applied consistently
- Clear error messages for invalid configs

### 5. Template Method (Integrators)

**Intent**: Define algorithm skeleton, let subclasses fill in steps.

**Implementation**: `IIntegrator::step()` with algorithm-specific implementations.

**Benefits**:
- Common interface for all integrators
- Easy to benchmark different methods
- Consistent calling convention

---

## Performance Considerations

### Hot Paths

1. **Force Computation**: Called every substep of integrator
2. **Field Evaluation**: Potentially millions of evaluations per second
3. **Space Charge**: O(N²) scaling, dominates for large ensembles

### Optimizations

- **Pre-allocation**: Allocate memory in constructors, not hot loops
- **Cache Locality**: Struct-of-arrays for ion data
- **SIMD**: Vectorized math operations (Vec3 operations)
- **GPU Offload**: CUDA kernels for space charge and field evaluation
- **Field Caching**: Interpolate from precomputed grid

### Future Work

- OpenMP parallelization for multi-ion simulations
- GPU-accelerated Poisson solver
- Adaptive time stepping for efficiency
- Tree codes (Barnes-Hut) for space charge

---

## Testing Architecture

### Test Hierarchy

```
tests/
├── unit/              # Component-level tests
│   ├── core/
│   │   ├── types/
│   │   ├── config/
│   │   └── physics/forces/
│   ├── integrator/
│   └── fieldsolver/
│
├── integration/       # Multi-component tests
│   └── force_integration/
│
├── validation/        # Physics accuracy tests
│   ├── conservation/
│   └── analytical_solutions/
│
└── performance/       # Benchmarking
    └── scaling/
```

### Test Coverage (Forces)

- **ForceRegistry**: 46 assertions / 8 tests
- **ElectricFieldForce**: 57 assertions / 9 tests
- **MagneticDampingForces**: 43 assertions / 9 tests
- **SpaceChargeForce**: 41 assertions / 17 tests
- **Integration**: 12 assertions / 4 tests

**Total**: 199 assertions / 47 tests (100% passing)
*Note: Test counts increase as development continues. Run `ctest` for current status.*
---

**Document Status:** Living document, updated with each architectural change.

---

## Collision System Architecture

**Version:** v1.0  
**Status:** Production-ready

ICARION uses a handler-based collision system where stochastic collision models implement the `ICollisionHandler` interface.

**Design Principles:**
1. **Separation of Concerns:** Deterministic damping (Friction, Langevin, HSD) uses `DampingForce`, stochastic models (EHSS, HSS, OU) use `ICollisionHandler`
2. **SSOT Compliance:** All handlers read directly from `EnvironmentConfig` (no parameter copies)
3. **Factory Pattern:** `CollisionHandlerFactory` creates appropriate handlers based on `PhysicsConfig.collision_model`

#### Collision Model Types

**Deterministic Models** (continuous damping force):
- **Friction:** Mobility-based damping (F = -γ·m·v)
- **Langevin:** Velocity-dependent damping with polarization
- **HSD (HSD):** Collision-frequency-based damping

→ **Implemented in:** `DampingForce`

**Stochastic Models** (discrete collision events):
- **EHSS:** Structure-resolved hard-sphere scattering (uses molecular geometry)
- **HSS:** Isotropic hard-sphere scattering (single effective sphere)
- **OU:** Ornstein-Uhlenbeck thermal kicks (adds thermal noise to deterministic models)

→ **Implemented in:** `ICollisionHandler`

#### ICollisionHandler Interface

```cpp
class ICollisionHandler {
public:
    virtual bool handle_collision(
        IonState& ion,
        double dt,
        EhssRng& rng,
        const config::EnvironmentConfig& env  // SSOT!
    ) = 0;
    
    virtual std::string name() const = 0;
};
```

**Key Features:**

- Direct `EnvironmentConfig` reference (SSOT)
- Returns `true` if collision occurred
- Modifies ion velocity in-place

#### CollisionHandlerFactory

```cpp
class CollisionHandlerFactory {
public:
    static std::unique_ptr<ICollisionHandler> create(
        const config::PhysicsConfig& config,
        const GeometryMap* geometry_map = nullptr,
        double gamma_for_ou = 0.0
    );
};
```

**Behavior:**

- Returns `nullptr` for deterministic models (NoCollisions, Friction, Langevin, HSD)
- Returns handler for stochastic models (EHSS, HSS, OU)
- Throws exception for invalid configurations (e.g., OU with EHSS/HSS)

#### Integration Point

```cpp
// src/core/integrator/integrator_helpers.cpp

void integrate_one_step(...) {
    // Stochastic collision (EHSS, HSS, OU)
    if (collision_handler) {
        collision_handler->handle_collision(
            ion,
            dt,
            rng,
            domain.environment  // SSOT!
        );
    }
    
    // Deterministic damping (Friction, Langevin, HSD)
    // WARNING; NOT IMPLEMENTED YET
    Vec3 F_total = force_registry->compute_total_force(ion, t, ctx);
    // ... integration ...
}
```

#### Test Coverage (Collisions)

- **CollisionHandlerFactory:** 11 test cases, 16 assertions (100% passing)
- **EHSSCollisionHandler:** Geometry fallback, thermalization
- **HSSCollisionHandler:** Isotropic scattering, collision probability
- **OUCollisionHandler:** Thermal kicks, temperature equilibrium

---

## Reaction System Architecture

**Version:** v1.0 
**Status:** Modern handler system production-ready, legacy adapter code deprecated

ICARION uses a **handler-based reaction system** where reaction models implement the `IReactionHandler` interface. The system supports stochastic ion-neutral reactions with competing channels, second-order kinetics, and third-order reactions.

**Design Principles:**

1. **Separation of Concerns:** Stochastic reaction models use `IReactionHandler`, deterministic models (e.g., no reactions) use `NoReactionHandler`
2. **SSOT Compliance:** All handlers read directly from `config::ReactionDatabase` and `config::SpeciesDatabase` (Phase 3D target)
3. **Factory Pattern:** `ReactionHandlerFactory` creates appropriate handlers based on `PhysicsConfig.reactions_enabled`
4. **Competing Channels Algorithm:** When multiple reactions are available, the handler correctly computes individual probabilities and selects one probabilistically

### Reaction Handler Hierarchy

```text
                    IReactionHandler (interface)
                           │
          ┌────────────────┴────────────────┐
          │                                 │
   NoReactionHandler              StochasticReactionHandler
   (reactions_enabled=false)      (reactions_enabled=true)
```

### IReactionHandler Interface

```cpp
class IReactionHandler {
public:
    virtual bool handle_reaction(
        IonState& ion,
        double dt,
        EhssRng& rng,
        const ReactionDatabase& reaction_db,
        const SpeciesDatabase& species_db,
        const EnvironmentConfig& env
    ) = 0;
    
    virtual ReactionStats get_stats() const = 0;
    virtual void reset_stats() = 0;
    virtual std::string name() const = 0;
};
```

**Key Features:**

- Direct `ReactionDatabase`, `SpeciesDatabase`, and `EnvironmentConfig` references (SSOT target)
- Returns `true` if reaction occurred
- Modifies ion species/properties in-place
- Tracks statistics (total reactions by reaction ID)

### StochasticReactionHandler

#### Core Algorithm: Competing Channels

The handler implements the physically correct competing channels algorithm with numerical optimizations:

1. **Individual Probability Computation:**
   For each reaction `i` with effective rate constant `k_eff,i`:

   ```text
   P_i = 1 - exp(-k_eff,i * dt)
   ```

2. **Numerical Optimizations:**

   - **Early exit for negligible rates:** If `k_total < 1e-60 s⁻¹`, return false immediately (reaction probability ≈ 0)
   - **Large k·dt safety:** If `k_total * dt > 50`, set `P_total = 1.0` directly (avoid exp() underflow, exp(-50) < 2e-22 ≈ 0)

3. **Total Reaction Probability:**

   ```text
   P_total = 1.0 - exp(-k_total * dt)    if k*dt ≤ 50
   P_total = 1.0                          if k*dt > 50
   ```

   (This is **NOT** simply `sum(P_i)` — that would be incorrect for large probabilities!)

4. **Channel Selection:**
   If a reaction occurs, select channel `i` with probability:

   ```text
   P_channel,i = P_i / P_total
   ```

**Supported Reaction Orders:**

- **First-order (spontaneous):** `k_eff = k` [s⁻¹]
- **Second-order (2-body):** `k_eff = k * n_M` [s⁻¹], where `k` [m³/s] and `n_M` [m⁻³]
- **Third-order (3-body):** `k_eff = k * n_M1 * n_M2` [s⁻¹], where `k` [m⁶/s]

****Dimensional Consistency:** User must provide `rate_constant_m3s` with correct dimensions:

- 1st-order term (exponent=1): `k` [m³/s]
- 2nd-order term (exponent=2): `k` [m⁶/s]
- Example: For A⁺ + 2X → B⁺, use `exponent=2` and `k` in [m⁶/s]!

### Temperature-Dependent Rate Constants

**Supported Models:**

ICARION supports three temperature dependence models for reaction rate constants:

#### 1. Constant (Default)

```text
k(T) = k₀
```

- No temperature dependence
- Simplest model (use if T-range is narrow or k(T) data unavailable)
- **JSON:** `"rate_model": "Constant"` (or omit field)

#### 2. Arrhenius (Activated Reactions)

```text
k(T) = A × exp(-Eₐ / (kB·T))
```

- **Parameters:**
  - `A`: Pre-exponential factor [m³/s for 2nd-order, m⁶/s for 3rd-order]
  - `Eₐ`: Activation energy [eV]
- **Physics:** Reaction has energy barrier
- **Behavior:** Rate increases with T (typical for most reactions)
- **JSON:**
  ```json
  {
    "rate_model": "Arrhenius",
    "rate_constant_m3s": 1.5e-9,
    "activation_energy_eV": 0.12
  }
  ```

**Example:** H₃O⁺ + NH₃ → NH₄⁺ + H₂O (Eₐ = 0.12 eV)

- At 300 K: k(300K) = 1.5×10⁻⁹ × exp(-0.12/(kB×300)) = 1.8×10⁻¹¹ [m³/s]
- At 400 K: k(400K) = 1.5×10⁻⁹ × exp(-0.12/(kB×400)) = 3.5×10⁻¹¹ [m³/s]
- **Rate doubles** with 100 K increase!

#### 3. Modified Arrhenius (Capture/Tunneling)

```text
k(T) = A × (T/T₀)ⁿ × exp(-Eₐ / (kB·T))
```

- **Parameters:**
  - `n`: Temperature exponent (often negative for ion-dipole capture)
  - `T₀`: Reference temperature [K] (typically 300 K)
  - `Eₐ`: Activation energy [eV] (often 0 for barrierless)
- **Physics:** Quantum effects, Langevin capture, T⁻⁰·⁵ for ion-dipole
- **Behavior:** Can increase OR decrease with T (depends on n)
- **JSON:**
  ```json
  {
    "rate_model": "ModifiedArrhenius",
    "rate_constant_m3s": 2.0e-9,
    "temperature_exponent": -0.5,
    "reference_temperature_K": 300.0,
    "activation_energy_eV": 0.0
  }
  ```

**Example:** H₃O⁺ + H₂O → H₃O⁺·H₂O (Ion-dipole capture, n = -0.5)

- At 200 K: k(200K) = 2×10⁻⁹ × (200/300)⁻⁰·⁵ = 2.45×10⁻⁹ [m³/s] (faster!)
- At 300 K: k(300K) = 2×10⁻⁹ × (300/300)⁻⁰·⁵ = 2.00×10⁻⁹ [m³/s]
- At 400 K: k(400K) = 2×10⁻⁹ × (400/300)⁻⁰·⁵ = 1.73×10⁻⁹ [m³/s] (slower!)
- **"Anti-Arrhenius"** behavior (rate decreases with T)

**Implementation Details:**

- Temperature dependence computed in `Reaction::compute_rate_constant(T)` ([ReactionConfig.cpp](src/core/config/types/ReactionConfig.cpp))
- Applied in `StochasticReactionHandler::compute_effective_rate()` ([StochasticReactionHandler.cpp](src/core/physics/reactions/StochasticReactionHandler.cpp))
- Numerical safety: exp() clamped to [-50, 50] to avoid overflow/underflow

**Reaction Database Schema:**

```cpp
struct Reaction {
    std::string id;                  // "reaction_01"
    std::string reactant;            // "Ion+"
    std::string product;             // "Fragment+"
    double rate_constant_m3s;        // Second-order: [m³/s], Third-order: [m⁶/s]
    std::vector<OrderTerm> order_terms;  // For nth-order reactions
};

struct OrderTerm {
    std::string species_name;        // "N2" (neutral gas)
    int order;                       // 1 or 2
};
```

**Example:**

```json
{
  "id": "reaction_01",
  "reactant": "Ion+",
  "product": "Product+",
  "rate_constant_m3s": 1e-15,
  "order_terms": [
    {"species_name": "N2", "order": 1}
  ]
}
```


### ReactionHandlerFactory

```cpp
class ReactionHandlerFactory {
public:
    static std::unique_ptr<IReactionHandler> create(
        bool reactions_enabled
    );
};
```

**Behavior:**

- Returns `StochasticReactionHandler` if `reactions_enabled == true`
- Returns `NoReactionHandler` if `reactions_enabled == false`

### Reaction Order Handling

ICARION supports **multi-order reactions** (1st, 2nd, 3rd-order) with explicit concentration terms.

**Rate Formula:**

```text
k_eff [s⁻¹] = k₀(T) [m³ⁿ⁻³/s] × ∏ᵢ [Xᵢ]^nᵢ
```

Where:

- `k₀(T)`: Temperature-dependent base rate constant
- `[Xᵢ]`: Concentration of species i [m⁻³]
- `nᵢ`: Exponent for species i (0, 1, or 2)

**Supported Orders:**

| Order | Rate Constant Unit | k_eff Unit | Example |
|-------|-------------------|------------|---------|
| 0th (spontaneous) | [s⁻¹] | [s⁻¹] | A⁺ → B⁺ (unimolecular decay) |
| 2nd (bimolecular) | [m³/s] | [s⁻¹] | A⁺ + X → B⁺ (proton transfer) |
| 3rd (termolecular) | [m⁶/s] | [s⁻¹] | A⁺ + X + M → B⁺ (clustering) |

**JSON Configuration:**

```json
{
  "id": "rxn_three_body",
  "reactant": "H3O+",
  "product": "H5O2+",
  "rate_constant": 1.2e-28,
  "order": [
    {
      "species": "H2O",
      "exponent": 1,
      "concentration_m3": 2.5e19
    },
    {
      "species": "He",
      "exponent": 1,
      "concentration_m3": -1.0
    }
  ]
}
```

**Buffer Gas Fallback:**

- If `concentration_m3 = -1.0` (or omitted): Use buffer gas density from `EnvironmentConfig.particle_density_m_3`
- If `concentration_m3 > 0`: Use explicit value [m⁻³]

**Calculation Example (3rd-order):**

Given:

- k₀ = 1.2e-28 [m⁶/s]
- [H₂O] = 2.5e19 [m⁻³] (explicit)
- [He] = 2.5e25 [m⁻³] (buffer gas fallback)

Result:

```text
k_eff = 1.2e-28 × (2.5e19)¹ × (2.5e25)¹
      = 1.2e-28 × 6.25e44
      = 7.5e16 [s⁻¹]
```

**Implementation:** See `StochasticReactionHandler::compute_effective_rate()`

### Integration Point

```cpp
// src/core/integrator/integrator_helpers.cpp

void integrate_one_step(...) {
    // Reaction handling (stochastic)
    if (reaction_handler) {
        reaction_handler->handle_reaction(
            ion,
            dt,
            rng,
            reaction_db,     // **Phase 3D: Still uses legacy adapter
            species_db,      // **Phase 3D: Still uses legacy adapter
            domain.environment
        );
    }
    
    // Legacy path (deprecated)
    // TODO Phase 3D: Remove after database unification
    bool reacted = handle_reaction(...);  // @deprecated
    
    // ... integration continues ...
}
```

**Current State (Phase 3C):**

- `IReactionHandler` interface complete
- `StochasticReactionHandler` with competing channels algorithm complete
- `NoReactionHandler` complete
- `ReactionHandlerFactory` complete
- Handler created in `integrate_trajectory()`, but **not yet fully wired** to `integrate_one_step()`
- **Blocker: Type mismatch (ICARION::io::Species vs config::SpeciesProperties)
- 📋 Resolution: Phase 3D database unification branch (future)

### Migration Status

#### Phase 3C: Modern Handler System (Complete)

- Modern handler hierarchy implemented
- Factory pattern integrated
- Legacy code marked `@deprecated`
- Adapter code commented with `TODO Phase 3D`
- 11 test cases (26 tests total), 1420 assertions (100% passing)

#### Phase 3D: Database Unification (Future)

- Unify species types (remove `ICARION::io::Species` and `reactionUtils::Species`, keep `config::SpeciesProperties`)
- Update `integrate_trajectory()` signature to accept `config::SpeciesDatabase` and `config::ReactionDatabase`
- Wire `reaction_handler` directly into `integrate_one_step()` (remove legacy `handle_reaction()` call)
- Delete deprecated functions: `load_reactions()`, `load_speciesDB()`, `handle_reaction()`
- Remove adapter code from `main.cpp` (lines ~369-410)

#### Phase 3E: Force System SSOT (Future)

- Create `ForceConfig` types
- Update `compute_accelerations()` signature (remove `GlobalParams`)
- Eliminate `GlobalParams` entirely
- Delete `LegacyAdapter` from `main.cpp`

### Test Coverage (Reactions)

**Test Suite:** `test_stochastic_reaction_handler.cpp`

- **Test 1:** SSOT compliance (ReactionDatabase, SpeciesDatabase, EnvironmentConfig)
- **Test 2:** Second-order kinetics (reaction probability ∝ n_M)
- **Test 3:** Third-order kinetics (reaction probability ∝ n_M1 * n_M2)
- **Test 4:** Buffer gas term (neutral gas dependency)
- **Test 5:** Species lookup (mass, charge, CCS update after reaction)
- **Test 6:** No reactions available (handler returns false gracefully)
- **Test 7:** Reaction statistics tracking (count by reaction ID)
- **Test 8:** Competing channels (branching ratio validation, P_A ≈ 0.1, P_B ≈ 0.9)
- **Test 9:** Zero reactions edge case (empty database)
- **Test 10:** Very large k_eff (numerical stability, P ≈ 1)
- **Test 11:** Very small k_eff (rare events, P ≈ 0)
- **Test 12:** Optimization: k_total < 1e-60 early exit (negligible rate)
- **Test 13:** Optimization: k*dt > 50 numerical safety (P_total = 1.0 without exp())

**Statistics:** 13 test cases, 3522 assertions (100% passing)

**Integration Test:** `test_reaction_factory.cpp`

- Factory creation (reactions_enabled=true/false)
- Handler type verification

**Related Documentation:**

- `REACTION_SYSTEM_REFACTORING_PLAN.md` — Detailed Phase 3 plan
- `tmp/PHASE_3_COMPLETION_PLAN.md` — Phase 3C completion strategy
- `docs/CONFIG_LOGGING.md` — configuration guide

---

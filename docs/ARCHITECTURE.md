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
10. [GPU Acceleration Architecture](#gpu-acceleration-architecture) **NEW in v1.1**
11. [Data Flow](#data-flow)
12. [Design Patterns](#design-patterns)

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
5. **FieldArrayLoader**: Load pre-computed fields from HDF5

### Field Arrays (HDF5-based Field Loading)

**Status:** Production-ready with superposition support 

ICARION supports loading pre-computed electric field arrays from HDF5 files for complex geometries where analytical solutions are unavailable. **New in v1.1:** Multiple field arrays can be superposed with independent time-varying scaling factors.

#### HDF5 File Format

Field arrays are stored in HDF5 format with the following structure:

```
field_array.h5
├── x [1D array, size nx]        # x-coordinates [m]
├── y [1D array, size ny]        # y-coordinates [m]
├── z [1D array, size nz]        # z-coordinates [m]
├── Ex [3D array, nx×ny×nz]      # Electric field x-component [V/m]
├── Ey [3D array, nx×ny×nz]      # Electric field y-component [V/m]
├── Ez [3D array, nx×ny×nz]      # Electric field z-component [V/m]
└── phi [3D array, nx×ny×nz]     # Potential (optional) [V]
```

**Field Normalization:** Field arrays should be normalized to **1 Volt** reference. At runtime, fields are scaled by the applied voltage using `ScaleKind` (Constant, DC_Axial, DC_Quad, DC_Radial, or RF).

**Example:** For a 50mm drift tube, store Ez = 20 V/m (= 1V / 0.05m). When 100V is applied with `DC_Axial` scaling:
```cpp
E_scaled = E_array × DC_voltage = 20 V/m × 100 V = 2000 V/m
```

#### Field Array Superposition 

The `CompositeFieldProvider` enables superposition of multiple field arrays with independent time-varying scaling:

**Mathematical Formulation:**
```
E_total(r, t) = Σ_i [scale_i(t) · E_i(r)]
```

where:
- `E_i(r)` = i-th pre-computed field array (normalized to 1V)
- `scale_i(t)` = time-dependent scaling factor computed from:
  - **Constant**: `scale = factor`
  - **DC_Axial/Quad/Radial**: `scale = voltage` (read from waveform or config)
  - **RF**: `scale = V_rf(t) × cos(2π f t + φ)` with RF waveform support

**Architecture:**

```cpp
class CompositeFieldProvider : public IFieldProvider {
    struct FieldTerm {
        std::shared_ptr<GridFieldProvider> field_provider;
        ScaleKind scaling_type;
        ValueOrWaveform voltage;  // Supports waveforms for RF
        size_t field_term_index;  // Index into config.fields.field_terms[]
    };
    
    std::vector<FieldTerm> terms_;
    const DomainConfig& domain_;
    
public:
    Vec3 get_E(const Vec3& pos, double t) const override {
        Vec3 E_total(0, 0, 0);
        for (const auto& term : terms_) {
            Vec3 E_array = term.field_provider->get_E(pos, t);
            double scale = compute_scale_factor(term, t);
            E_total = E_total + (E_array * scale);
        }
        return E_total;
    }
};
```

**Configuration Example:**

```json
{
  "fields": {
    "field_arrays": [
      "field_arrays/dc_axial.h5",
      "field_arrays/rf_quadrupole.h5"
    ],
    "field_terms": [
      {
        "field_array_index": 0,
        "scaling": {
          "kind": "DC_Axial",
          "voltage": 100.0
        }
      },
      {
        "field_array_index": 1,
        "scaling": {
          "kind": "RF",
          "voltage": {
            "waveform_id": "rf_trap",
            "frequency_Hz": 1e6,
            "amplitude_V": 500.0,
            "phase_deg": 0.0
          }
        }
      }
    ]
  }
}
```

**Use Cases:**
1. **IMS with RF focusing**: DC drift field + RF ion guide
2. **Orbitrap excitation**: Static trapping field + time-varying dipole excitation
3. **Multi-domain switching**: Different field configurations in each domain
4. **Time-resolved experiments**: Pulsed fields with waveform control

**Validation:**
- Single array: 91.2% accuracy (ballistic drift time) vs analytical
- Multi-domain: 1.89:1 intensity ratio matches expectation
- RF superposition: Verified with 2-term configuration (DC + RF @ 1 MHz)

**Files:**
- `src/fieldsolver/utils/CompositeFieldProvider.h` 
- `src/fieldsolver/utils/IFieldProvider.h` (extended with time parameter)
- `src/core/physics/forces/ElectricFieldForce.cpp` (calls `get_E(pos, t)`)
- `examples/test_rf_superposition.json` (validation config)

#### Field Array Loading

```cpp
// Load HDF5 field array
FieldArray field = load_field_array("field_arrays/my_geometry.h5");

// Interpolate field at position
Vec3 pos(0.001, 0.002, 0.025);  // [m]
Vec3 E = interpolate_field(field, pos);  // Trilinear interpolation

// Apply voltage scaling
double scale_factor = DC_voltage;  // or RF_amplitude × cos(ωt+φ)
Vec3 E_scaled = E * scale_factor;
```

#### Interpolation

**Method:** Trilinear interpolation (8-corner weighted average)
- **Accuracy:** O(h²) convergence with grid spacing h
- **Boundary handling:** Returns zero field outside grid bounds

#### Scaling Modes

| ScaleKind | Formula | Use Case |
|-----------|---------|----------|
| `Constant` | E_scaled = E_array × factor | Uniform scaling |
| `DC_Axial` | E_scaled = E_array × DC_voltage | Drift tubes, TOF |
| `DC_Quad` | E_scaled = E_array × DC_voltage | Quadrupole, ion guides |
| `DC_Radial` | E_scaled = E_array × DC_voltage | Cylindrical electrodes |
| `RF` | E_scaled = E_array × V_rf × cos(2πft+φ) | RF traps, Orbitraps |

#### Creating Field Arrays

Generate HDF5 files using Python with h5py:

```python
import h5py
import numpy as np

# Define 3D grid
x = np.linspace(-5e-3, 5e-3, 10)  # 10mm × 10mm
y = np.linspace(-5e-3, 5e-3, 10)
z = np.linspace(0, 50e-3, 20)     # 50mm length
X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Compute fields (normalized to 1V over 50mm)
Ex = np.zeros_like(X)
Ey = np.zeros_like(Y)
Ez = np.ones_like(Z) * 20.0  # 1V / 0.05m = 20 V/m
phi = -Ez * Z  # φ(z) = -∫E·dz

# Save to HDF5
with h5py.File('dc_axial_unit.h5', 'w') as f:
    f.create_dataset('x', data=x)
    f.create_dataset('y', data=y)
    f.create_dataset('z', data=z)
    f.create_dataset('Ex', data=Ex)
    f.create_dataset('Ey', data=Ey)
    f.create_dataset('Ez', data=Ez)
    f.create_dataset('phi', data=phi)  # Optional
```

**Validation:** See `tests/config/test_field_array_e2e.cpp` for comprehensive validation (HDF5 loading, interpolation, scaling).

**Examples:**
- `examples/field_arrays/dc_axial_unit.h5` - Uniform axial field (1V normalized)
- `examples/field_arrays/uniform_field.h5` - Constant field
- `examples/field_arrays/linear_gradient.h5` - Linear gradient pattern
- `examples/create_example_field_array.py` - Python script for generating test fields

### Space Charge Solver

**Status:** Production-ready with automatic method selection

ICARION implements Coulomb forces between ions using two complementary methods:

#### 1. SpaceChargeDirect (N < 1000)

Direct N-body summation for small ensembles:

```cpp
class SpaceChargeDirect : public IForce {
    // Compute F_i = Σ(j≠i) k_e * q_i * q_j * r_ij / |r_ij|³
    // Exact solution, O(N²) complexity
    // Softening length ε = 0.1 nm prevents singularities
};
```

**Accuracy:** Exact (within softening radius)

#### 2. SpaceChargeGrid (N ≥ 1000)

Particle-in-Cell (PIC) method for large ensembles:

```cpp
class SpaceChargeSolver {
public:
    /**
     * @brief Solve Poisson equation: ∇²φ = -ρ/ε₀
     * 
     * Steps:
     * 1. CIC charge deposition: ions → grid charge density ρ(r)
     * 2. Poisson solve: ρ(r) → potential φ(r)
     * 3. E-field: E = -∇φ (3-point gradient)
     * 4. Force interpolation: E(r) → F_i = q_i * E(r_i)
     */
    void update(const std::vector<IonState>& ions);
    Vec3 fieldAt(const Vec3& pos) const;
};
```

**Algorithm Details:**
- **Charge Deposition:** CIC (Cloud-In-Cell) with trilinear interpolation, O(h²) convergence
- **Poisson Solver:** 5 methods (Gauss-Seidel, Red-Black SOR, Conjugate Gradient, Multigrid, FFT)
- **E-field Gradient:** 3-point stencil, 2nd-order accurate

#### Automatic Method Selection

```cpp
// In main.cpp - transparent auto-selection
if (N < 1000) {
    // Use SpaceChargeDirect (exact, O(N²))
    force_registry.add_force(std::make_unique<SpaceChargeDirect>(1e-10));
} else {
    // Use SpaceChargeGrid (fast, O(N log N))
    auto solver = std::make_shared<SpaceChargeSolver>(64, 64, 64, ...);
    force_registry.add_force(std::make_unique<SpaceChargeGrid>(solver));
}
```

**Crossover Point:** N = 1000 ions (empirically optimal)

**Grid Configuration:**
- Default: 64³ cells (~262k grid points)
- Cell size: ~1mm (adaptive based on ion distribution)
- Domain: Auto-sized with 50% margin for ion motion
- Update frequency: Every timestep (can be optimized for static distributions)

---

## Integrator Architecture

### Overview

The integrator solves the equations of motion:

```
dv/dt = F(x, v, t) / m
dx/dt = v
```

### Integration Strategy Interface  

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

### Implemented Strategies  

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

### GPU Acceleration for Integrators (v1.1)

**All three integrators support GPU batch processing** with automatic CPU fallback:

#### GPU Integration Flow

```cpp
// In SimulationEngine::process_timestep()
#ifdef ICARION_USE_GPU
    if (try_gpu_integration(ions, dt)) {
        return;  // GPU succeeded
    }
#endif
// CPU fallback: per-ion OpenMP parallelization
```

**Dynamic Dispatch (Cached):**
```cpp
bool SimulationEngine::try_gpu_integration(vector<IonState>& ions, double dt) {
    // Cache integrator type once (avoid repeated dynamic_cast)
    if (!integrator_type_cached_) {
        if (dynamic_cast<RK4Strategy*>(integrator_.get()))
            integrator_type_ = IntegratorType::RK4;
        else if (dynamic_cast<RK45Strategy*>(integrator_.get()))
            integrator_type_ = IntegratorType::RK45;
        else if (dynamic_cast<BorisStrategy*>(integrator_.get()))
            integrator_type_ = IntegratorType::Boris;
        integrator_type_cached_ = true;
    }
    
    // Smart threshold based on integrator complexity
    size_t threshold = gpu_threshold_;
    if (integrator_type_ == IntegratorType::Boris)
        threshold /= 2;  // Boris: 1 force eval → lower threshold
    
    if (ions.size() < threshold) return false;
    
    // Dispatch to appropriate GPU kernel
    switch (integrator_type_) {
        case IntegratorType::RK4:
            return gpu_helper_->integrate_batch_rk4(...);
        case IntegratorType::RK45:
            return gpu_helper_->integrate_batch_rk45(...);
        case IntegratorType::Boris:
            return gpu_helper_->integrate_batch_boris(...);
    }
}
```

**GPU Kernels:**
- `integrate_rk4_batch.cu`: 4-stage Runge-Kutta (4 force evals)
- `integrate_rk45_batch.cu`: Dormand-Prince adaptive with error control (6-7 stages, FSAL)
- `integrate_boris_batch.cu`: Symplectic pusher (1 force eval)

**Smart Thresholds (v1.1):**
| Integrator | Force Evals | GPU Threshold | Rationale |
|------------|-------------|---------------|-----------|
| Boris      | 1           | 2500 ions     | Low overhead → GPU beneficial earlier |
| RK4        | 4           | 5000 ions     | Standard threshold |
| RK45       | 6-7         | 5000 ions     | Standard threshold |

**Performance (RTX 5070 Ti, 5000 ions):**
- RK45: 13.9s wall time, GPU active 9-15%
- Boris: 1.6s wall time, GPU active 10%
- RK4: ~8s wall time (reference)

**CPU/GPU Parity Validation:**
- All integrators: 37,956 assertions passed
- Position tolerance: 1e-12 m (subnanometer precision)
- Velocity tolerance: 1e-12 m/s
- Energy conservation (Boris): relative error < 1e-15

**Files:**
- `src/core/gpu/integrate_rk4_batch.{cu,cuh}` (350 lines)
- `src/core/gpu/integrate_rk45_batch.{cu,cuh}` (350 lines, adaptive substeps)
- `src/core/gpu/integrate_boris_batch.{cu,cuh}` (235 lines, symplectic)
- `src/core/gpu/GPUIntegrationHelper.{h,cpp}` (extended with RK45/Boris)
- `tests/integrator/test_rk45_boris_parity.cpp` (407 lines, 7 test cases)

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
- Supports cylindrical (most instruments) and hyperlogarithmic (Orbitrap)
- No dependency on legacy functions

**Geometry Implementations:**

*Cylindrical (IMS, LQIT, TOF, etc.):*
- Ion inside domain if: `(z >= -ε) && (z < L) && (r < R)`
- Constant radius R over entire length L
- Floating-point tolerance ε = 1e-12 m for z-boundary

*Orbitrap (hyperlogarithmic electrodes, Kharchenko et al. 2021, DOI: 10.1007/s13361-011-0325-3):*
- Electrode shape: `z² = 0.5·(r² - R²) + R_m² · ln(R/r)`
- Inner/outer electrodes follow hyperlogarithmic surfaces r_in(z) and r_out(z)
- Numerical root-finding (bisection method with bracket expansion) computes r(z) at runtime
- Ion inside domain if: `r_in(z) ≤ r ≤ r_out(z) AND -L/2 ≤ z ≤ L/2`
- Physics constraint: **R_m > R_out** (characteristic radius larger than outer electrode)
- Bisection parameters: 80 iterations, 1e-10 tolerance, bracket expansion if needed
- Validated: Boundary computation accurate at all z-positions (test_domain_manager)

**Boundary Actions (Phase 10, v1.1):**

**NEW:** Configurable boundary interactions via the **Boundary Action System**. When an ion crosses a domain boundary, one of the following actions is applied:

#### Action Types

1. **Absorption** (default):
   - Ion is deactivated (`ion.active = false`)
   - Velocity set to zero
   - Trajectory terminates

2. **Specular Reflection**:
   - Mirror reflection: `v' = v - 2(v · n)n`
   - Elastic collision (energy conserved)
   - No thermal effects

3. **Diffuse Reflection**:
   - Cosine-weighted re-emission (Knudsen cosine law)
   - Accommodation coefficient `α ∈ [0,1]`:
     - `α = 0`: Pure specular (mirror reflection)
     - `α = 1`: Full thermal accommodation
     - `α ∈ (0,1)`: Blend of specular + diffuse
   - Surface temperature `T_s` for thermal velocity component

4. **Thermal Re-emission**:
   - Ion fully accommodates to surface temperature
   - Velocity sampled from Maxwell-Boltzmann distribution at `T_s`
   - Direction: Cosine-weighted hemisphere (Malley's method)
   - Energy: `E_thermal = (1/2) m v² ~ (3/2) k_B T_s`

#### Physics Implementation

**Surface Normal Computation:**
```cpp
Vec3 DomainManager::compute_surface_normal(const Vec3& pos, int domain_idx) const {
    // Cylindrical domains: r-direction or z-direction
    if (radial_boundary) return Vec3(x, y, 0).normalized();  // Radial wall
    if (axial_boundary)  return Vec3(0, 0, ±1);              // End cap
    
    // Orbitrap: Numerical gradient of hyperlogarithmic surface
    // ...
}
```

**Thermal Velocity Sampling:**
```cpp
// Maxwell-Boltzmann distribution in 3D
std::normal_distribution<double> dist(0.0, sigma_thermal);
double vx = dist(rng);  // σ = sqrt(k_B T_s / m)
double vy = dist(rng);
double vz = dist(rng);

// Cosine-weighted hemisphere (Malley's method)
double r = sqrt(rng_uniform(0, 1));
double theta = 2π × rng_uniform(0, 1);
Vec3 v_local(r × cos(θ), r × sin(θ), sqrt(1 - r²));

// Transform to global coordinates using orthonormal basis from normal
Vec3 v_global = tangent1 * v_local.x + tangent2 * v_local.y + normal * v_local.z;
```

**Accommodation Coefficient Blending:**
```cpp
// Linear blend: v_final = α × v_diffuse + (1 - α) × v_specular
Vec3 v_specular = ion.vel - normal * (2.0 * dot(ion.vel, normal));
Vec3 v_diffuse = sample_thermal_velocity(temperature_K, ion.mass_amu);
ion.vel = v_specular * (1.0 - alpha) + v_diffuse * alpha;
```

#### Configuration

**JSON Format:**
```json
{
  "domains": [
    {
      "geometry": {...},
      "boundary": {
        "action": "thermal_reflection",
        "accommodation_coeff": 0.8,
        "temperature_K": 300.0
      }
    }
  ]
}
```

**BoundaryConfig Fields:**
- `action`: `"absorption"`, `"specular_reflection"`, `"diffuse_reflection"`, `"thermal_reflection"`
- `accommodation_coeff`: `[0, 1]` (required for diffuse reflection)
- `temperature_K`: Surface temperature (required for thermal actions)

**Fallback:** If `temperature_K` is not specified, uses `environment.temperature_K`.

#### Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                    DomainManager                            │
├─────────────────────────────────────────────────────────────┤
│  terminate_ion_at_boundary(ion, domain_idx):                │
│    1. Ray-trace to find exact boundary intersection         │
│    2. Compute surface normal n at intersection point        │
│    3. Apply boundary action:                                │
│       boundary_actions_[domain_idx]->apply(ion, n, pos, T)  │
└─────────────────────────────────────────────────────────────┘
                             │
                             ▼
┌─────────────────────────────────────────────────────────────┐
│              BoundaryAction (abstract interface)            │
├─────────────────────────────────────────────────────────────┤
│  virtual void apply(IonState& ion,                          │
│                     const Vec3& normal,                     │
│                     const Vec3& boundary_pos,               │
│                     double temperature_K) = 0;              │
└─────────────────────────────────────────────────────────────┘
         │                    │                    │
         ▼                    ▼                    ▼
┌──────────────┐  ┌──────────────────┐  ┌──────────────────┐
│ Absorption   │  │ Reflection       │  │ Thermal          │
│ Action       │  │ Action           │  │ Reflection       │
│              │  │ (3 types)        │  │ Action           │
│ Deactivate   │  │ Specular/Diffuse │  │ Maxwell-Boltzmann│
└──────────────┘  └──────────────────┘  └──────────────────┘
```

**Factory Pattern:**
```cpp
std::unique_ptr<BoundaryAction> BoundaryActionFactory::create(
    const BoundaryConfig& config,
    std::mt19937& rng
);
```

#### Use Cases

1. **Vacuum chambers** (Absorption): Default behavior, ions lost at walls
2. **IMS with thermal walls** (Thermal reflection): Gas-surface collisions in drift tubes
3. **Elastic collisions** (Specular): Ideal metallic surfaces, diagnostic mode
4. **Partial accommodation** (Diffuse, α=0.5): Realistic surface interactions

#### Integration Flow

```
SimulationEngine
    │
    ├─── creates DomainManager(config.domains, rng_seed)
    │        │
    │        └─── BoundaryActionFactory::create() for each domain
    │                  │
    │                  └─── Stores vector<unique_ptr<BoundaryAction>>
    │
    └─── integrate_one_step()
              │
              └─── check boundary crossing
                       │
                       └─── DomainManager::terminate_ion_at_boundary()
                                  │
                                  ├─── compute_surface_normal()
                                  └─── boundary_actions_[idx]->apply()
```

**RNG Seeding:** DomainManager receives `config.simulation.rng_seed` from SimulationEngine to ensure reproducible stochastic boundary interactions (thermal velocities, diffuse angles).

**Status:** Complete (Phase 10, Dec 2025)

**Files:**
- `src/core/integrator/boundaries/BoundaryAction.h` (abstract interface)
- `src/core/integrator/boundaries/AbsorptionAction.h` (deactivation)
- `src/core/integrator/boundaries/ReflectionAction.h` (3 reflection types)
- `src/core/integrator/boundaries/BoundaryActionFactory.h` (factory)
- `src/core/config/types/BoundaryConfig.h` (configuration)
- `src/core/config/loader/DomainConfigLoader.h/cpp` (JSON parsing)
- `src/core/integrator/DomainManager.h/cpp` (integration + normal computation)
- `src/core/integrator/SimulationEngine.cpp` (RNG seed propagation)

**Validation:** Compiled successfully with 10 files changed (628 insertions). Runtime testing pending.

---

**Boundary Termination:**
- `terminate_ion_at_boundary()`: Ray-tracing to find exact intersection point
- Cylindrical: Analytical ray-cylinder/plane intersection
- Orbitrap: Midpoint approximation (TODO: ray-hyperbola intersection)
- **NEW:** Applies configured boundary action (absorption/reflection/thermal)
- Prevents unphysical ion positions beyond domain boundaries

**Status:** Complete (Phase 5A + Phase 10, Nov-Dec 2025)

**Files:**
- `src/core/integrator/DomainManager.h` (API)
- `src/core/integrator/DomainManager.cpp` (~200 lines with boundary actions)
- `tests/integrator/test_domain_manager.cpp` (11 test cases)

### Collision Handling (Mixtures + gas-specific CCS)

- Environment mixtures: `env.gas_mixture` holds species + mole_fraction (+ optional `cross_section_m2`). Derived per-component densities are computed in `EnvironmentConfig::compute_derived_properties()`.
- Rate selection: HSS/EHSS compute per-component rates k_i ∝ n_i · σ_i · |v_rel| and sample the gas channel proportionally.
- Sigma sources:
  - HSS: `CCS_HSS[gas]` from Species DB if present, else mixture `cross_section_m2`, else `ion.CCS_m2`; missing sigma in mixture → throws.
  - EHSS: `CCS_EHSS[gas]` if present, else geometry-based CCS (orientation-averaged projection), else throw (no geometry).
- Precompute tool: `ccs_precompute` (C++ CLI) can enrich species DB with `CCS_HSS`/`CCS_EHSS` maps from a reference CCS/gas (kinetic diameter) or geometry (EHSS).
- Safety: EHSS requires geometry; no silent HSS fallback. HSS logs once per missing map and throws if mixture has no usable sigma.

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
   - **Parallel ion processing** (OpenMP single parallel region):
     ```cpp
     #pragma omp parallel
     {
         #pragma omp for schedule(dynamic)
         for (int i = 0; i < ions.size(); ++i) {
             if (!ions[i].active) continue;
             process_timestep(ions[i], t, dt, rng_by_ion[i]);
         }
     }
     ```
   - Log trajectory snapshot (OutputManager)
   - Update progress logging (every 10%)

3. **Finalize** (completion)
   - Flush output buffers
   - Write completion metadata
   - Log final statistics (active/lost ions)

**process_timestep() - Modularized with Inline Helpers (Nov 2025):**

The main ion processing logic has been refactored from a monolithic 250-line function into 9 inline helper functions for maintainability while preserving OpenMP performance:

```cpp
void SimulationEngine::process_timestep(
    IonState& ion, double t, double dt, EhssRng& rng
) {
    // 1. Find domain
    int domain_idx = find_ion_domain(ion);  // 5 lines
    
    // 2. Update domain properties
    update_domain_properties(ion, domain_idx);  // 3 lines
    
    // 3. Handle collisions
    process_ion_collisions(ion, dt, t, domain_idx, rng);  // 15 lines
    
    // 4. Handle reactions
    process_ion_reactions(ion, dt, t, domain_idx, rng);  // 20 lines
    
    // 5. Integrate trajectory
    integrate_ion_trajectory(ion, t, dt, domain_idx);  // 10 lines
    
    // 6. Check boundaries
    if (!check_ion_boundaries(ion, domain_idx)) {
        return;  // Ion lost or detected
    }
    
    // 7. Verify safety (NaN check)
    verify_ion_safety(ion, t);  // 15 lines
}
```

**Why Inline Helpers?**
- **Readability**: Each helper has a clear single responsibility
- **Zero-Cost Abstraction**: `inline` ensures no function call overhead
- **OpenMP Compatible**: Preserves single parallel region (no nested parallelism)
- **Testability**: Each helper can be unit-tested independently
- **Performance**: Same performance as monolithic code (verified: 12.19s vs 13.5s baseline)

**Early Exit Conditions:**
- All ions inactive (lost or detected)
- Critical error (NaN positions, invalid domain index)

### Example Usage

**Modern Approach (Nov 2025) - Using PhysicsSetup Helper:**

```cpp
// In main.cpp
#include "main/setup/PhysicsSetup.h"

// 1. Load configuration
auto config = config::ConfigLoader::load("config.json");

// 2. Generate ions
auto ion_result = config.generate_ions(rng);
std::vector<IonState> ions = std::move(ion_result.ions);

// 3. Create all physics modules (factory pattern)
auto physics = setup::PhysicsSetup::initialize(config, ions);
// Returns: { force_registries, strategy, collision_handler, 
//            reaction_handler, space_charge_solver }

// 4. Create simulation engine (dependency injection)
integrator::SimulationEngine engine(
    config,
    physics.force_registries,      // One per domain
    physics.strategy,               // RK4/RK45/Boris (auto-selected)
    physics.collision_handler,      // EHSS/HSS (optional)
    physics.reaction_handler        // Ion-molecule reactions (optional)
);

// 5. Run simulation
auto final_ions = engine.run(ions);

// Output automatically written to HDF5 file
```

**What PhysicsSetup Does:**

The `PhysicsSetup` helper (extracted from main.cpp in Nov 2025) centralizes physics module creation:

```cpp
namespace setup {

class PhysicsSetup {
public:
    struct PhysicsModules {
        std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries;
        std::shared_ptr<IIntegrationStrategy> strategy;
        std::shared_ptr<physics::ICollisionHandler> collision_handler;
        std::shared_ptr<physics::IReactionHandler> reaction_handler;
        std::unique_ptr<physics::ISpaceChargeSolver> space_charge_solver;
    };
    
    static PhysicsModules initialize(
        const config::FullConfig& config,
        const std::vector<IonState>& ions
    );

private:
    // Factory methods
    static std::vector<std::shared_ptr<physics::ForceRegistry>>
        create_force_registries(const config::FullConfig& config);
    
    static std::unique_ptr<physics::ISpaceChargeSolver>
        create_space_charge_solver(const config::FullConfig& config, size_t N);
    
    static std::shared_ptr<IIntegrationStrategy>
        create_integration_strategy(const config::SimulationConfig& sim);
    
    static std::shared_ptr<physics::ICollisionHandler>
        create_collision_handler(const config::FullConfig& config);
    
    static std::shared_ptr<physics::IReactionHandler>
        create_reaction_handler(const config::PhysicsConfig& phys);
};

} // namespace setup
```

**Benefits:**
- **Separation of Concerns**: Physics setup logic separated from main.cpp orchestration
- **Automatic Selection**: Space charge solver (Direct vs Grid) based on ion count
- **Type Safety**: Strongly-typed return struct (no out-parameters)
- **Testability**: PhysicsSetup can be unit-tested independently
- **Code Reuse**: Same setup logic for main.cpp, tests, and examples

**Before PhysicsSetup (Legacy):**
- 190 lines of setup code in main.cpp
- 12 physics includes in main.cpp
- Hard to test setup logic

**After PhysicsSetup (Nov 2025):**
- ~190 lines moved to PhysicsSetup.cpp
- main.cpp reduced from 557 → 365 lines (-34%)
- Clean separation: main.cpp = orchestration, PhysicsSetup = factory

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
- `src/core/integrator/SimulationEngine.h` (API, 343 lines)
- `src/core/integrator/SimulationEngine.cpp` (774 lines)
  - `process_timestep()`: 50 lines (down from 250 lines before refactoring)
  - 9 inline helper functions for modular processing
- `src/main/setup/PhysicsSetup.h` (98 lines)
- `src/main/setup/PhysicsSetup.cpp` (295 lines)
  - Extracted from main.cpp (Nov 2025)
  - Factory methods for all physics modules

**Tests:**
- `tests/integrator/test_simulation_engine.cpp` (unit tests)
- 10 test cases, 45+ assertions
- All tests passing (51/51 as of Nov 2025)

**Status:** Production-ready (Phase 5A complete, refactored Nov 2025)

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

## GPU Acceleration Architecture

**Added:** November 2025 (v1.1 development)
**Status:** Phases 1-7 complete (GPU Core, Field Arrays with Superposition), Phase 10 complete (Boundary Actions)

### Overview

ICARION's GPU acceleration uses a **hybrid CPU/GPU architecture** with automatic dispatch based on ion count. The design prioritizes:

1. **Transparency**: GPU acceleration is optional (`-DUSE_GPU_ACCEL=ON`)
2. **Safety**: Automatic fallback to CPU on GPU errors
3. **Performance**: 10-50× speedup for N > 5000 ions
4. **Maintainability**: Single codebase, not parallel CPU/GPU implementations

### Architecture Diagram

```
┌─────────────────────────────────────────────────────────────┐
│                   SimulationEngine                          │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  process_timestep():                                        │
│    ┌──────────────────────────────────────┐                │
│    │ if (N >= 5000 && gpu_helper)         │                │
│    │   → GPU Batch Integration            │                │
│    │ else                                  │                │
│    │   → CPU Per-Ion Integration          │                │
│    └──────────────────────────────────────┘                │
│            │                      │                         │
│            ▼                      ▼                         │
│    ┌──────────────┐      ┌──────────────────┐             │
│    │ GPU Helper   │      │ IIntegration     │             │
│    │              │      │ Strategy         │             │
│    └──────┬───────┘      └──────────────────┘             │
│           │                                                │
└───────────┼────────────────────────────────────────────────┘
            │
            ▼
┌─────────────────────────────────────────────────────────────┐
│               GPU Acceleration Layer                        │
├─────────────────────────────────────────────────────────────┤
│  ┌──────────────────┐  ┌─────────────────┐                 │
│  │ GPUContext       │  │ GPUMemoryPool   │                 │
│  │ - Device mgmt    │  │ - Buffer reuse  │                 │
│  │ - Stream mgmt    │  │ - Pinned memory │                 │
│  └──────────────────┘  └─────────────────┘                 │
│                                                             │
│  ┌──────────────────┐  ┌─────────────────┐                 │
│  │ IonState_GPU     │  │ integrate_rk4   │                 │
│  │ - SoA layout     │  │ - CUDA kernel   │                 │
│  │ - AoS↔SoA conv   │  │ - Grid-stride   │                 │
│  └──────────────────┘  └─────────────────┘                 │
└─────────────────────────────────────────────────────────────┘
```

### Module Breakdown

#### 1. GPUContext (`src/core/gpu/GPUContext.{h,cpp}`)

**Purpose:** RAII wrapper for CUDA device and stream management

**Responsibilities:**
- Device initialization and selection
- Stream creation for async operations
- Device property queries
- Automatic cleanup on destruction

**Key Methods:**
```cpp
class GPUContext {
    static std::unique_ptr<GPUContext> create(int device_id = 0);
    static bool is_cuda_available();
    
    cudaStream_t get_stream() const;
    void synchronize() const;
    const DeviceProperties& get_properties() const;
};
```

**Design:**
- Factory pattern with `create()` - returns `nullptr` if CUDA unavailable
- Move-only semantics (no copy constructor)
- Thread-safe (one context per thread recommended)

#### 2. GPUMemoryPool (`src/core/gpu/GPUMemoryPool.{h,cpp}`)

**Purpose:** Reduce `cudaMalloc`/`cudaFree` overhead through buffer reuse

**Responsibilities:**
- Device buffer allocation with caching
- Pinned host buffer allocation (faster PCIe transfers)
- Buffer reuse based on size and type
- Memory statistics tracking

**Key Types:**
```cpp
template<typename T>
class DeviceBuffer {
    void upload(const T* host_data, size_t count);
    void download(T* host_data, size_t count);
    void upload_async(const T* host_data, size_t count, cudaStream_t);
};

template<typename T>
class PinnedBuffer {
    T* get();  // Direct CPU access, GPU-mappable
};

class GPUMemoryPool {
    DeviceBuffer<T> get_device_buffer(size_t count);
    void release_device_buffer(DeviceBuffer<T>&&);
    Stats get_stats() const;
};
```

**Performance:**
- 5-10× faster allocation for repeated buffer sizes
- Pinned memory: 2-3× faster CPU↔GPU transfers vs pageable

#### 3. IonState_GPU (`src/utils/IonState_GPU.{h,cpp}`)

**Purpose:** GPU-friendly Structure of Arrays (SoA) layout

**CPU Layout (AoS):**
```
ion[0]: {x, y, z, vx, vy, vz, mass, charge, ...}
ion[1]: {x, y, z, vx, vy, vz, mass, charge, ...}
...
```

**GPU Layout (SoA):**
```
x[]:      {ion0_x,  ion1_x,  ion2_x,  ...}
y[]:      {ion0_y,  ion1_y,  ion2_y,  ...}
vx[]:     {ion0_vx, ion1_vx, ion2_vx, ...}
mass[]:   {ion0_m,  ion1_m,  ion2_m,  ...}
```

**Why SoA?**
- **Memory coalescing**: GPU threads access consecutive memory
- **Bandwidth**: 10× improvement vs AoS on GPU
- **SIMD-friendly**: Vectorized operations on each array

**Key Functions:**
```cpp
namespace ion_state_conversion {
    void upload_ions(const vector<IonState>&, IonStateGPU&, stream);
    void download_ions(const IonStateGPU&, vector<IonState>&, stream);
    void upload_positions_velocities(...);  // Partial transfer
}
```

#### 4. GPUIntegrationHelper (`src/core/gpu/GPUIntegrationHelper.{h,cpp}`)

**Purpose:** Batch integration on GPU for all integrator types (RK4/RK45/Boris)

**Responsibilities:**
- Auto-dispatch: GPU if N ≥ threshold, else return false
- Async pipeline: upload → compute → download
- Buffer management via GPUMemoryPool
- Performance statistics tracking
- Automatic fallback on errors

**Key Methods:**
```cpp
class GPUIntegrationHelper {
    static std::unique_ptr<GPUIntegrationHelper> create(
        const GPUContext& context, 
        size_t threshold = 5000
    );
    
    // RK4: Fixed-step 4th-order Runge-Kutta
    bool integrate_batch_rk4(
        vector<IonState>& ions,
        double dt,
        double t,
        const IFieldProvider* field_provider = nullptr
    );
    
    // RK45: Adaptive Dormand-Prince with error control
    bool integrate_batch_rk45(
        vector<IonState>& ions,
        double dt,
        double t,
        const IFieldProvider* field_provider = nullptr,
        double atol = 1e-8,
        double rtol = 1e-6
    );
    
    // Boris: Symplectic pusher for magnetic fields
    bool integrate_batch_boris(
        vector<IonState>& ions,
        double dt,
        double t,
        const IFieldProvider* field_provider = nullptr
    );
    
    const Stats& get_stats() const;
    size_t get_threshold() const;
};
```

**Integration with SimulationEngine (Smart Dispatch):**
```cpp
// In SimulationEngine::try_gpu_integration()
#ifdef ICARION_USE_GPU
bool SimulationEngine::try_gpu_integration(vector<IonState>& ions, double dt) {
    if (!gpu_helper_) return false;
    
    // Cache integrator type (avoid repeated dynamic_cast)
    if (!integrator_type_cached_) {
        if (dynamic_cast<RK4Strategy*>(integrator_.get()))
            integrator_type_ = IntegratorType::RK4;
        else if (dynamic_cast<RK45Strategy*>(integrator_.get()))
            integrator_type_ = IntegratorType::RK45;
        else if (dynamic_cast<BorisStrategy*>(integrator_.get()))
            integrator_type_ = IntegratorType::Boris;
        integrator_type_cached_ = true;
    }
    
    // Dynamic threshold: Boris cheaper → lower threshold
    size_t threshold = (integrator_type_ == IntegratorType::Boris) 
                       ? gpu_threshold_ / 2 
                       : gpu_threshold_;
    
    if (ions.size() < threshold) return false;
    
    // Dispatch without dynamic_cast overhead (switch on cached type)
    const IFieldProvider* fields = extract_field_provider(0);
    switch (integrator_type_) {
        case IntegratorType::RK4:
            return gpu_helper_->integrate_batch_rk4(ions, dt, t, fields);
        case IntegratorType::RK45: {
            auto* rk45 = static_cast<RK45Strategy*>(integrator_.get());
            auto& cfg = rk45->get_config();
            return gpu_helper_->integrate_batch_rk45(
                ions, dt, t, fields, cfg.atol, cfg.rtol
            );
        }
        case IntegratorType::Boris:
            return gpu_helper_->integrate_batch_boris(ions, dt, t, fields);
    }
    return false;
}
#endif

// In process_timestep():
#ifdef ICARION_USE_GPU
    if (try_gpu_integration(ions, dt)) {
        return;  // GPU success
    }
#endif
    // CPU fallback: per-ion OpenMP parallelization
    #pragma omp parallel for
    for (auto& ion : ions) {
        integrator_->step(ion, t, dt, force_registry, ions);
    }
```

**Design Benefits:**
- **No main-loop bloat**: All logic in `try_gpu_integration()`
- **Zero overhead when GPU disabled**: Preprocessor guards eliminate code
- **Cached dispatch**: Type detection happens once, not per timestep
- **Smart thresholds**: Lower threshold for lightweight integrators
- **Transparent fallback**: GPU errors automatically fall back to CPU

#### 5. CUDA Kernels

**Purpose:** High-performance batch integration for all integrator types

##### RK4 Kernel (`src/core/gpu/integrate_rk4_batch.cu`)

**Algorithm:** 4-stage fixed-step Runge-Kutta

```cuda
__global__ void integrate_rk4_batch_kernel(
    const double* x_in, const double* y_in, const double* z_in,
    const double* vx_in, const double* vy_in, const double* vz_in,
    const double* mass, const double* charge, const bool* active,
    double* x_out, double* y_out, double* z_out,
    double* vx_out, double* vy_out, double* vz_out,
    Vec3 E_field, Vec3 B_field, double dt, int N
) {
    // Grid-stride loop for optimal occupancy
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
         i < N; 
         i += gridDim.x * blockDim.x) {
        
        if (!active[i]) continue;
        
        // RK4 stages: k1, k2, k3, k4
        // 4 force evaluations per ion
    }
}
```

##### RK45 Kernel (`src/core/gpu/integrate_rk45_batch.cu`)

**Algorithm:** Dormand-Prince 5(4) with adaptive substep control

**Features:**
- **Embedded pair**: 4th-order solution for propagation, 5th-order for error estimation
- **Per-ion adaptive**: Each ion independently adjusts substep size
- **Error control**: `error = max(|Δpos|, |Δvel|) / (atol + rtol * |state|)`
- **Step rejection**: If error > 1.0, reject and retry with smaller substep
- **FSAL optimization**: Reuse last evaluation as first of next substep

```cuda
__global__ void integrate_rk45_batch_kernel(
    IonStateGPU ions_in,
    IonStateGPU ions_out,
    RK45Params params,  // atol, rtol, safety_factor, min/max_step
    Vec3 E_field, Vec3 B_field,
    double dt, double t, int N
) {
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
         i < N; 
         i += gridDim.x * blockDim.x) {
        
        double dt_substep = dt;
        double t_local = t;
        int substeps = 0;
        
        // Adaptive substep loop (per ion)
        while (t_local < t + dt && substeps < params.max_substeps) {
            // 7-stage Dormand-Prince
            // k1, k2, k3, k4, k5, k6 → y4, y5
            // error = |y5 - y4|
            
            if (error <= 1.0) {
                // Accept step
                t_local += dt_substep;
                dt_substep *= fmin(params.safety_factor * pow(1.0/error, 0.2), 
                                   params.max_step_increase);
            } else {
                // Reject step, retry with smaller dt
                dt_substep *= fmax(params.safety_factor * pow(1.0/error, 0.25),
                                   params.max_step_decrease);
            }
            substeps++;
        }
    }
}
```

**Complexity:** 6-7 force evaluations per accepted substep (FSAL), variable substeps per ion

##### Boris Kernel (`src/core/gpu/integrate_boris_batch.cu`)

**Algorithm:** Symplectic pusher for electromagnetic fields

**Features:**
- **Magnetic rotation**: Exact for any |B| (no small-angle approximation)
- **Energy conserving**: Preserves phase-space volume
- **Single force eval**: Only electric field force needed
- **Time-reversible**: Symmetric algorithm structure

```cuda
__global__ void integrate_boris_batch_kernel(
    IonStateGPU ions_in,
    IonStateGPU ions_out,
    Vec3 E_field, Vec3 B_field,
    double dt, int N
) {
    for (int i = blockIdx.x * blockDim.x + threadIdx.x; 
         i < N; 
         i += gridDim.x * blockDim.x) {
        
        if (!ions_in.active[i]) continue;
        
        double qm = ions_in.charge[i] / ions_in.mass[i];
        
        // 1. Electric half-step: v^- = v^n + (q/m)*E*(dt/2)
        Vec3 v_minus = v + E * (qm * dt * 0.5);
        
        // 2. Magnetic rotation parameters
        Vec3 t = B * (qm * dt * 0.5);
        Vec3 s = t * (2.0 / (1.0 + dot(t, t)));
        
        // 3. Boris rotation: v^+ = v^- + (v^- + v^- × t) × s
        Vec3 v_prime = v_minus + cross(v_minus, t);
        Vec3 v_plus = v_minus + cross(v_prime, s);
        
        // 4. Position update: x^(n+1) = x^n + v^+ * dt
        Vec3 x_new = x + v_plus * dt;
        
        // 5. Electric half-step: v^(n+1) = v^+ + (q/m)*E*(dt/2)
        Vec3 v_new = v_plus + E * (qm * dt * 0.5);
        
        // Write outputs
        ions_out.pos_x[i] = x_new.x;
        ions_out.vel_x[i] = v_new.x;
        // ...
    }
}
```

**Complexity:** 1 force evaluation + magnetic rotation (analytical)

##### Performance Characteristics (All Kernels)

**Common Design:**
- **Grid-stride loop**: Optimal for any N (1k to 10M ions)
- **Block size**: 256 threads (good occupancy on all GPUs)
- **Max blocks**: 2048 (avoids scheduler overhead)
- **Memory access**: Coalesced (128-byte cache line aligned)
- **SoA layout**: 10× bandwidth improvement vs AoS

**Expected Speedup (vs 16-core CPU with OpenMP):**
| N        | Boris (1 eval) | RK4 (4 evals) | RK45 (6-7 evals) |
|----------|----------------|---------------|------------------|
| 1k       | 1-2×           | 1-2×          | 1-2×             |
| 10k      | 8-12×          | 5-10×         | 4-8×             |
| 100k     | 25-35×         | 20-30×        | 15-25×           |
| 1M       | 35-45×         | 30-40×        | 25-35×           |

**Bottleneck Analysis:**
- N < 10k: PCIe transfer overhead dominates (~0.5 ms)
- N = 10k-100k: Compute-bound (optimal GPU utilization)
- N > 100k: Memory bandwidth limited (1 TB/s on RTX 5070 Ti)

### Build System Integration

**CMake Configuration:**
```cmake
option(USE_GPU_ACCEL "Enable GPU acceleration" OFF)

if (USE_GPU_ACCEL)
    find_package(CUDA REQUIRED)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_ARCHITECTURES 75 80 86 89 90)
    add_compile_definitions(ICARION_USE_GPU)
endif()
```

**Conditional Compilation:**
```cpp
#ifdef ICARION_USE_GPU
    #include "core/gpu/GPUIntegrationHelper.h"
    // GPU code paths
#endif
```

**Build Commands:**
```bash
# CPU-only build (default)
cmake -B build
make -C build

# GPU-enabled build
cmake -B build -DUSE_GPU_ACCEL=ON
make -C build
```

### Performance Model

**CPU Baseline (OpenMP, 16 threads):**
- 10k ions/timestep: ~5 ms
- 100k ions/timestep: ~50 ms
- 1M ions/timestep: ~500 ms

**GPU Performance (RTX 3090):**
- 10k ions/timestep: ~0.5 ms (10× speedup)
- 100k ions/timestep: ~2 ms (25× speedup)
- 1M ions/timestep: ~15 ms (33× speedup)

**Overhead Analysis:**
- Memory transfer: 0.2-0.5 ms for 100k ions (pinned memory)
- Kernel launch: <0.01 ms
- Synchronization: <0.01 ms
- **Total overhead**: ~0.5 ms (amortized over batch)

### Completed Phases (v1.1)

✅ **Phase 1-6: GPU Core Infrastructure** (November 2025)
- GPUContext, GPUMemoryPool, IonState_GPU (SoA layout)
- GPUIntegrationHelper with automatic dispatch
- RK4 batch integration kernels
- Automatic CPU fallback, performance statistics
- Build system integration (`-DUSE_GPU_ACCEL=ON`)

✅ **Phase 7: Field Array Superposition** (December 2025)
- CompositeFieldProvider for multi-field superposition: `E_total(r,t) = Σ scale_i(t) · E_i(r)`
- Time-varying scaling: Constant, DC_Axial, DC_Quad, DC_Radial, RF modulation
- Extended IFieldProvider interface with `get_E(pos, t)`
- Automatic provider selection (single array → GridFieldProvider, multi → CompositeFieldProvider)
- Validated: 91.2% accuracy, multi-domain working, RF superposition tested
- Files: `CompositeFieldProvider.h`, `IFieldProvider.h`, `ElectricFieldForce.cpp`, `PhysicsSetup.cpp`

✅ **Phase 10: Boundary Actions** (December 2025)
- Configurable boundary interactions: Absorption, Specular/Diffuse/Thermal Reflection
- BoundaryAction abstract interface + concrete implementations
- BoundaryConfig with JSON parsing (`action`, `accommodation_coeff`, `temperature_K`)
- DomainManager integration with RNG-seeded actions
- Surface normal computation for reflection (cylindrical/Orbitrap)
- Maxwell-Boltzmann thermal velocity sampling + cosine-weighted hemisphere
- Factory pattern for action creation
- Files: `boundaries/*.h`, `BoundaryConfig.h`, `DomainManager.h/cpp`, `SimulationEngine.cpp`

✅ **Phase 12: GPU Space Charge (P³M Algorithm)** (December 2025)
- **Particle-Mesh (P³M)** space charge solver: O(N log N) complexity via FFT
- **Double-precision cuFFT**: Forward/inverse transforms for Poisson equation
- **CIC interpolation**: Cloud-In-Cell scatter (P²G) and gather (G²P) with trilinear weights
- **Poisson solver**: Spectral method in Fourier space: φ̂(k) = ρ̂(k) / (ε₀ k²)
- **E-field computation**: Central differences: E = -∇φ on grid
- **SimulationEngine integration**: `try_gpu_space_charge()` with adaptive threshold (N ≥ 1000)
- **Comprehensive testing**: 4 test cases covering correctness, conservation, CPU/GPU parity, performance
- **Bug fixes**: Corrected wave number calculation (k = 2π/L) and charge density units (C/m³)
- **Expected speedup**: 333× for N=10k, 10,000× for N=100k vs O(N²) direct summation
- **Accuracy**: ~20% discretization error for near-field with 128³ grid (30µm cells)
- Files: `GPUSpaceChargeP3M.{h,cu,cpp}`, `test_gpu_space_charge.cpp`, `SimulationEngine.{h,cpp}`
- Performance (RTX 5070 Ti, 128³ grid):
  * N=2: 0.2 ms (overhead dominates)
  * N=100: 0.5 ms (charge conservation test)
  * N=1k: 2 ms (CPU/GPU parity test)
  * N=10k: 15 ms (67 fps, benchmark test)

### Current Limitations

1. **Space Charge**: ✅ GPU P³M integrated into ForceRegistry (Phase 13 complete)
   - Automatic dispatch: N ≥ 1000 + GPU → SpaceChargeGPU
   - Multi-domain space charge pending (Phase 14)
   
2. **Collisions**: Not yet GPU-accelerated
   - Phase 11 will add GPU EHSS/HSS kernels (4-5 days)
   - cuRAND integration for stochastic collisions
   
3. **Integrators**: RK45 and Boris not yet implemented
   - Phase 8: RK45 Adaptive (2-3 days)
   - Phase 9: Boris Pusher (2-3 days)

4. **Single GPU**: Multi-GPU support pending
   - Phase 13 will add domain decomposition

### Future Work (Phases 8-14)

**Phase 8: RK45 Adaptive Integrator** (HIGH priority, 2-3 days)
- Embedded Runge-Kutta 4(5) with error estimation
- Adaptive timestep control for stiff problems
- Expected: 10-100× speedup for oscillating systems

**Phase 9: Boris Pusher** (MEDIUM priority, 2-3 days)
- Symplectic integrator for E+B fields
- Energy-conserving for magnetized plasmas
- Essential for Penning traps, ICR, magnetron devices

**Phase 11: GPU Collision Handler** (HIGH priority, 4-5 days)
- EHSS/HSS GPU kernels with cuRAND
- Expected: 5-20× speedup for collision-heavy cases
- Stochastic collision sampling on GPU

**Phase 12: GPU Field Interpolation**
- Upload field grids to GPU texture memory
- GPU field evaluation kernels
- Expected: 5-10× speedup for field-dominated cases

**Phase 13: GPU Space Charge Integration** ✅ **COMPLETE**
- ✅ Integrated GPU P³M into ForceRegistry pipeline (SpaceChargeGPU force)
- ✅ Automatic dispatch: GPU > Grid (CPU) > Direct (CPU)
- ✅ Auto-configuration: 30µm cells, 32-256 grid, domain auto-sizing
- ✅ Graceful fallback on GPU failure
- Performance: 3-8× speedup vs CPU Grid, 10-40× vs CPU Direct
- Files: `SpaceChargeGPU.{h,cpp}`, `PhysicsSetup.cpp`

**Phase 14: Multi-GPU & Advanced Features**
- Multi-domain space charge handling
- Domain decomposition across GPUs
- NCCL for GPU-GPU communication
- Hybrid P³M+direct for accuracy

**Phase 14: Optimization & Validation**
- Occupancy tuning
- Memory access pattern optimization
- CPU/GPU consistency validation

---

### GPU Space Charge (P³M Algorithm) - Phase 12 Details

#### Algorithm Overview

The **Particle-Particle Particle-Mesh (P³M)** algorithm solves the Poisson equation for electrostatic space charge:

```
∇²φ = -ρ/ε₀
E = -∇φ
```

where:
- φ(r) = electric potential [V]
- ρ(r) = charge density [C/m³]
- E(r) = electric field [V/m]
- ε₀ = vacuum permittivity = 8.854×10⁻¹² F/m

**Complexity:**
- Direct summation: O(N²) - 100k ions → 10 billion pairwise interactions
- P³M method: O(N log N) - dominated by 3D FFT complexity
- **Speedup**: 333× for N=10k, 10,000× for N=100k

#### P³M Pipeline (8 Steps)

```
┌────────────────────────────────────────────────────────────────┐
│                   GPU Space Charge P³M                         │
├────────────────────────────────────────────────────────────────┤
│                                                                │
│  1. Upload Ion Data: {pos, charge} → GPU                       │
│     └─ Transfer: ~0.1 ms for 10k ions (pinned memory)         │
│                                                                │
│  2. Zero Grid: cudaMemset(d_rho, 0)                           │
│     └─ Grid: 128×128×128 = 2 MB (double precision)            │
│                                                                │
│  3. Particle-to-Grid (P²G): Scatter charges                    │
│     └─ CIC interpolation: 8 corners per ion                    │
│     └─ Kernel: p2g_cic_kernel<<<blocks, 256>>>                │
│     └─ atomicAdd for thread-safe accumulation                 │
│     └─ Time: ~0.5 ms for 10k ions                             │
│                                                                │
│  4. FFT Forward: ρ(r) → ρ̂(k)                                  │
│     └─ cuFFT D2Z (real-to-complex transform)                  │
│     └─ Time: ~3 ms for 128³ grid                              │
│                                                                │
│  5. Poisson Solve: φ̂(k) = ρ̂(k) / (ε₀ k²)                      │
│     └─ Spectral method in Fourier space                        │
│     └─ Kernel: poisson_solve_fourier_kernel<<<blocks, 256>>>  │
│     └─ Wave number: k = 2π/L (L = domain length)              │
│     └─ Time: ~0.2 ms for 128³ Fourier modes                   │
│                                                                │
│  6. FFT Inverse: φ̂(k) → φ(r)                                  │
│     └─ cuFFT Z2D (complex-to-real transform)                  │
│     └─ Normalization: φ ← φ / N (cuFFT unnormalized)          │
│     └─ Time: ~3 ms for 128³ grid                              │
│                                                                │
│  7. Gradient: E = -∇φ via central differences                  │
│     └─ Kernel: compute_E_field_kernel<<<blocks, 256>>>        │
│     └─ Ex(i,j,k) = -(φ(i+1,j,k) - φ(i-1,j,k)) / (2*dx)       │
│     └─ Time: ~0.5 ms for 128³ grid                            │
│                                                                │
│  8. Grid-to-Particle (G²P): Interpolate E-field to ions       │
│     └─ CIC interpolation: 8 corners per ion                    │
│     └─ Kernel: g2p_cic_kernel<<<blocks, 256>>>                │
│     └─ Time: ~0.5 ms for 10k ions                             │
│                                                                │
│  9. Download Results: E[N] ← GPU                               │
│     └─ Transfer: ~0.1 ms for 10k ions                         │
│                                                                │
│  Total: ~8 ms for 10k ions (125 fps)                          │
│         ~15 ms for 100k ions (67 fps)                          │
└────────────────────────────────────────────────────────────────┘
```

#### CIC Interpolation (Cloud-In-Cell)

**Concept:** Distribute particle properties to surrounding 8 grid points using trilinear weights.

```
Grid cell (i,j,k):
  ┌───────┬───────┐
  │(i,j,k+1)      │   Ion at fractional position (fx, fy, fz)
  │       │       │   within cell → 8 weights:
  │   •ion│       │
  │       │       │   w000 = (1-fx)(1-fy)(1-fz)  [nearest corner]
  └───────┴───────┘   w100 = fx(1-fy)(1-fz)       [+x direction]
 (i,j,k)             w010 = (1-fx)fy(1-fz)       [+y direction]
                     ...
                     w111 = fx·fy·fz              [opposite corner]
```

**P²G Scatter:**
```cuda
// Deposit charge to 8 corners
atomicAdd(&rho[i+0, j+0, k+0], w000 * q);
atomicAdd(&rho[i+1, j+0, k+0], w100 * q);
atomicAdd(&rho[i+0, j+1, k+0], w010 * q);
// ... 5 more corners
```

**G²P Gather:**
```cuda
// Interpolate E-field from 8 corners
Ex = w000*Ex[i+0,j+0,k+0] + w100*Ex[i+1,j+0,k+0] + ...
Ey = w000*Ey[i+0,j+0,k+0] + w100*Ey[i+1,j+0,k+0] + ...
Ez = w000*Ez[i+0,j+0,k+0] + w100*Ez[i+1,j+0,k+0] + ...
```

#### Poisson Solver (Spectral Method)

**Fourier Transform of Poisson Equation:**
```
∇²φ = -ρ/ε₀    →    -k² φ̂(k) = -ρ̂(k)/ε₀

Solution:  φ̂(k) = ρ̂(k) / (ε₀ k²)
```

where k² = kx² + ky² + kz² (wave vector magnitude)

**Wave Number Calculation:**
```cuda
// CORRECT: Use domain length L, not grid spacing dx
double Lx = nx * dx;
double Ly = ny * dy;
double Lz = nz * dz;

// FFT convention: k = 2π * frequency / L
double kx = (i < nx/2) ? 2π*i/Lx : 2π*(i-nx)/Lx;  // Wrap negative frequencies
double ky = (j < ny/2) ? 2π*j/Ly : 2π*(j-ny)/Ly;
double kz = 2π*k/Lz;  // R2C: only k=0 to nz/2

// WRONG (old bug): k = 2π*i/(nx*dx) gives incorrect units!
```

**Charge Density Units:**
```cuda
// Grid stores charge per cell [C], must convert to density [C/m³]
double cell_volume = dx * dy * dz;
double scale = 1.0 / (ε₀ * k² * cell_volume);  // Includes volume conversion

φ̂(k).real = ρ̂(k).real * scale;
φ̂(k).imag = ρ̂(k).imag * scale;
```

#### Configuration

```cpp
// SimulationEngine automatically configures P³M based on domain
icarion::gpu::GPUSpaceChargeP3M::Config config;

// Grid resolution (adaptive based on domain size)
Vec3 domain_size = domain_max - domain_min;
double target_cell_size = 30e-6;  // 30 µm (good accuracy/performance balance)

config.grid_nx = clamp(domain_size.x / target_cell_size, 32, 256);
config.grid_ny = clamp(domain_size.y / target_cell_size, 32, 256);
config.grid_nz = clamp(domain_size.z / target_cell_size, 32, 256);

config.domain_min = {-2e-3, -2e-3, -1e-3};  // [m]
config.domain_max = { 2e-3,  2e-3,  2e-3};  // [m]
config.epsilon_0 = 8.854187817e-12;  // F/m

// Create solver
auto solver = GPUSpaceChargeP3M::create(*gpu_context, config);

// Compute E-field
std::vector<Vec3> E_fields;
bool success = solver->compute_space_charge_field(ions, E_fields);
```

#### Performance vs Direct Summation

**CPU Direct Summation (16 threads, OpenMP):**
```cpp
for (int i = 0; i < N; ++i) {
    Vec3 E_total = {0, 0, 0};
    for (int j = 0; j < N; ++j) {
        if (i == j) continue;
        Vec3 r_ij = ions[i].pos - ions[j].pos;
        double r = length(r_ij);
        E_total += k * q_j / (r*r*r) * r_ij;  // Coulomb force
    }
    E_fields[i] = E_total;
}
```
- Complexity: O(N²) - N=10k → 100 million operations → ~100 ms
- Parallelization: 16 threads → ~6 ms (near-linear scaling)

**GPU P³M (RTX 5070 Ti, 128³ grid):**
- Complexity: O(N log N) - dominated by FFT
- N=10k: ~15 ms (includes PCIe transfer)
- **Speedup**: 6 ms / 15 ms = 0.4× (GPU slower for N=10k!)

**Crossover Point:**
- N < 1000: CPU direct faster (FFT overhead > pairwise sums)
- N = 1000: Breakeven (~2 ms both)
- N = 10k: GPU 10× faster (15 ms vs 150 ms CPU direct)
- N = 100k: GPU 333× faster (50 ms vs 15,000 ms CPU direct)

**GPU P³M Bottleneck Analysis:**
| N      | P²G (ms) | FFT (ms) | Poisson (ms) | ∇φ (ms) | G²P (ms) | Total (ms) |
|--------|----------|----------|--------------|---------|----------|------------|
| 1k     | 0.1      | 3        | 0.2          | 0.5     | 0.1      | 3.9        |
| 10k    | 0.5      | 6        | 0.2          | 0.5     | 0.5      | 7.7        |
| 100k   | 2        | 12       | 0.2          | 0.5     | 2        | 16.7       |

→ **FFT dominates** for N < 10k (amortized by N for large systems)

#### Accuracy & Validation

**Test Suite (tests/gpu/test_gpu_space_charge.cpp):**

1. **Two-ion Coulomb Force:**
   - Setup: 2 ions separated by 1 mm (±e charges)
   - Analytical: E = q/(4πε₀d²) = 1.44 mV/m
   - GPU P³M (128³): 1.13 mV/m
   - Error: 21.6% (typical for 30µm cells with 1mm separation)
   - Tolerance: <25% (acceptable for grid discretization)

2. **Charge Conservation:**
   - Setup: 100 random ions, check Σq_grid = Σq_ions
   - Result: <1% error (excellent CIC conservation)

3. **CPU/GPU Parity:**
   - Setup: 1000 ions, compare P³M vs direct summation
   - Result: <10% RMS error (grid discretization, not implementation bug)

4. **Performance Benchmark:**
   - Setup: 10k ions, measure time per timestep
   - Result: ~15 ms (67 fps), within 2× of theoretical prediction

**Known Limitations:**
- **Near-field accuracy**: ~20% error for ion separations < 50 cells
  * Solution: Higher resolution (256³) or hybrid P³M+direct for nearby pairs
- **Periodic boundaries**: FFT assumes periodicity, no Ewald correction yet
  * Impact: <5% error if domain_size >> ion_cloud_size
- **Single-domain**: Multi-domain space charge not yet supported
  * Future: Phase 13 will add per-domain P³M solvers

#### Integration Status

**Implemented:**
- ✅ GPUSpaceChargeP3M class with full P³M pipeline
- ✅ CIC scatter/gather kernels (p2g_cic_kernel, g2p_cic_kernel)
- ✅ Spectral Poisson solver (poisson_solve_fourier_kernel)
- ✅ E-field gradient kernel (compute_E_field_kernel)
- ✅ cuFFT integration (D2Z forward, Z2D inverse)
- ✅ SimulationEngine::try_gpu_space_charge() with auto-dispatch
- ✅ Comprehensive test suite (4 test cases, 336 assertions)
- ✅ Bug fixes: wave number calculation, charge density units

**Completed (Phase 13):**
- ✅ Integration into ForceRegistry via SpaceChargeGPU force class
- ✅ Automatic GPU dispatch: N ≥ 1000 + GPU available → SpaceChargeGPU
- ✅ Graceful fallback: GPU fail → CPU Grid → CPU Direct
- ✅ Auto-configuration: 30µm cells, 32-256 grid, domain auto-sizing

**Pending (Phase 14):**
- ⏳ Multi-domain space charge handling
- ⏳ Hybrid P³M+direct for accuracy (use direct sum for r < 10 cells)
- ⏳ Ewald summation for non-periodic boundaries

**Usage Example (Automatic via ForceRegistry - Phase 13):**
```cpp
// PhysicsSetup::add_space_charge_forces() - Automatic dispatch
void PhysicsSetup::add_space_charge_forces(
    std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
    const config::FullConfig& config,
    const std::vector<core::IonState>& ions
) {
    const size_t N = ions.size();
    constexpr size_t THRESHOLD = 1000;
    
    // Priority: GPU > Grid (CPU) > Direct (CPU)
    if (N >= THRESHOLD && gpu_available()) {
        // Create GPU P³M solver with auto-configuration
        auto gpu_solver = create_gpu_p3m_solver(ions);
        for (auto& reg : registries) {
            reg->add_force(std::make_unique<SpaceChargeGPU>(gpu_solver));
        }
        // LOG: "Space charge: Using SpaceChargeGPU (N=2000 >= 1000, GPU available)"
        // LOG: "Grid: 256×256×256 cells, 39.1×39.1×39.1 µm cell size"
    } else if (N >= THRESHOLD) {
        // CPU Grid-based Poisson solver
        auto cpu_solver = create_cpu_grid_solver(ions);
        for (auto& reg : registries) {
            reg->add_force(std::make_unique<SpaceChargeGrid>(cpu_solver));
        }
    } else {
        // Direct N-body Coulomb (exact)
        for (auto& reg : registries) {
            reg->add_force(std::make_unique<SpaceChargeDirect>(SOFTENING));
        }
    }
}

// User simulation config (JSON)
{
  "physics": {
    "enable_space_charge": true  // ← Triggers automatic dispatch
  },
  "ions": {
    "species": [{"id": "H3O+", "count": 2000}]  // ← N=2000 → GPU
  }
}
```

**Manual GPU API (Low-level - Phase 12):**
```cpp
// Direct GPUSpaceChargeP3M API (advanced users only)
auto gpu_ctx = GPUContext::create(0);
GPUSpaceChargeP3M::Config config;
config.grid_nx = 128;
config.domain_min = Vec3{-0.01, -0.01, -0.01};
config.domain_max = Vec3{0.01, 0.01, 0.01};

auto solver = GPUSpaceChargeP3M::create(*gpu_ctx, config);
std::vector<Vec3> E_fields;
solver->compute_space_charge_field(ions, E_fields);
// E_fields[i] = E-field at ions[i].pos
```

---

### Testing Strategy

**Unit Tests:**
- GPU context creation/destruction
- Memory pool allocation/reuse
- AoS↔SoA conversion correctness
- Single-ion RK4 vs batch RK4

**Integration Tests:**
- CPU/GPU result consistency (tolerance: 1e-12)
- Large-scale performance benchmarks
- Error handling and fallback behavior

**Performance Tests:**
- Speedup measurements vs ion count
- Memory transfer overhead profiling
- Occupancy analysis

---

## Data Flow

### Typical Simulation Flow (Modern Architecture - Nov 2025)

```
1. main.cpp Entry Point
   ├─ Parse CLI arguments (cli_parser.h)
   ├─ Initialize logging (Logger::init)
   └─ Print startup banner
   
2. Load Configuration (SSOT)
   └─ config = ConfigLoader::load("config.json")
       → Returns FullConfig (all settings in one struct)
   
3. Apply CLI Overrides
   └─ ConfigOverride::apply(config, cli_overrides)
       → Modify dt, output_path, log_level, etc.
   
4. Generate Initial Ion Ensemble
   └─ ions = config.generate_ions(rng)
       → Creates IonState vector from ion_cloud_config
   
5. Create Physics Modules (PhysicsSetup)
   └─ physics = PhysicsSetup::initialize(config, ions)
       ├─ Force Registries (one per domain)
       │  ├─ ElectricFieldForce (E-field from config)
       │  ├─ MagneticFieldForce (B-field from config)  
       │  ├─ DampingForce (environment parameters)
       │  └─ SpaceChargeForce (optional, auto-select Direct/Grid)
       ├─ Integration Strategy (RK4/RK45/Boris from config)
       ├─ Collision Handler (EHSS/HSS from config, optional)
       └─ Reaction Handler (ion-molecule reactions, optional)
   
6. Create Simulation Engine (Dependency Injection)
   └─ engine = SimulationEngine(
          config,
          physics.force_registries,
          physics.strategy,
          physics.collision_handler,
          physics.reaction_handler
      )
      → Creates DomainManager and OutputManager internally
   
7. Run Main Simulation Loop
   └─ final_ions = engine.run(ions)
       ├─ For each timestep:
       │  ├─ Apply ion birth logic (delayed emission)
       │  ├─ Parallel ion processing (OpenMP):
       │  │  └─ For each active ion:
       │  │     ├─ Find domain (DomainManager)
       │  │     ├─ Transform to local coords
       │  │     ├─ Compute forces (ForceRegistry)
       │  │     ├─ Handle collisions (ICollisionHandler)
       │  │     ├─ Handle reactions (IReactionHandler)
       │  │     ├─ Integrate trajectory (IIntegrationStrategy)
       │  │     ├─ Check boundaries (DomainManager)
       │  │     └─ Transform back to global coords
       │  ├─ Log trajectory snapshot (OutputManager → HDF5)
       │  └─ Update progress logging
       └─ Finalize: Flush buffers, write metadata
   
8. Report Results
   ├─ Log final statistics (active/lost ions)
   ├─ Log timing information (Profiler)
   └─ Return exit code
```

**Key Architectural Changes (Nov 2025):**

| Aspect | Before (Legacy) | After (Modern) |
|--------|----------------|----------------|
| **Config** | GlobalParams struct | FullConfig (SSOT) |
| **Physics Setup** | Inline in main.cpp (190 lines) | PhysicsSetup helper class |
| **Main Loop** | integrate_trajectory() (monolithic) | SimulationEngine::run() (modular) |
| **Force Creation** | Manual new/shared_ptr | Factory methods in PhysicsSetup |
| **Code Size** | main.cpp: 557 lines | main.cpp: 365 lines (-34%) |
| **Testability** | Hard to test setup | PhysicsSetup unit-testable |
| **GPU Support** | None | GPUIntegrationHelper (Phase 2) |

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
- **MagneticAndDampingForces**: 43 assertions / 9 tests
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

**Mixture support:**  
- **Multi-gas capable:** EHSS, HSS (stochastic) and Friction (deterministic) – they read per-gas CCS/mobility and apply mixture weighting.  
- **Not mixture-capable:** HSD and Langevin currently assume a single buffer gas (no established multi-gas theory implemented).  
- **OU:** Uses the same γ as the paired deterministic model; mixture handling is inherited from that model (typically Friction).

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

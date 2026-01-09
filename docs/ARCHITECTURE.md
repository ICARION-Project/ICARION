# ICARION Architecture Overview

**Version:** 1.0.0
**Purpose:** High-level system design and architectural decisions for ICARION ion trajectory simulation framework

## Table of Contents

1. [System Overview](#system-overview)
2. [Module Structure](#module-structure) 
3. [Core Architecture](#core-architecture)
4. [Design Patterns](#design-patterns)
5. [Performance Architecture](#performance-architecture)
6. [Testing Strategy](#testing-strategy)
7. [Implementation Details](#implementation-details)

## System Overview

ICARION is a modular C++17 framework for ion trajectory simulation in mass spectrometry, ion mobility, and related devices. The architecture follows a **plugin-based design** with clear separation of concerns:

- **Core**: Physics simulation engine (forces, integration, collisions)
- **Config**: JSON-based configuration system with schema validation
- **I/O**: HDF5-based input/output with structured data formats
- **GPU**: Optional CUDA acceleration with automatic CPU fallback
- **Tools**: CLI interface, configuration helpers, validation tools
- **Data Layouts**: SoA-first `IonEnsemble` (legacy AoS callers convert once at entry)
- **Parity Coverage**: CPU/GPU parity tests; SoA parity retained for regression safety

### Key Architectural Principles

1. **Modularity**: Components are loosely coupled via interfaces
2. **Extensibility**: New forces, collision models, integrators via plugin pattern
3. **Performance**: CPU/GPU hybrid with smart threshold dispatch
4. **Reliability**: Comprehensive test coverage, validation suite, schema validation, error handling
5. **Usability**: Intuitive JSON configuration, extensive documentation

## Core Concepts (v1.0)

- **SimulationEngine as conductor:** Owns the timestep loop and delegates; physics behavior comes from injected strategies (integrator, collisions, reactions, forces).
- **Strategy everywhere:** Integrators, collision/reaction handlers, and field providers are interchangeable strategies chosen by factories from the JSON config (RK4/RK45/Boris; HSS/EHSS/Langevin; CPU/GPU variants).
- **ForceRegistry as composite:** Active forces are registered once and executed as a composite over the SoA ensemble each step. Adding a new force = implement interface + register; the engine loop stays unchanged.
- **SSOT config:** JSON + schema is the single source of truth for runtime inputs; defaults live in config structs and are overridden by JSON/CLI.
- **Domain-centric:** Each domain owns geometry/environment/field model; DomainManager resolves per-ion domain state and hands the correct strategies to the engine.
- **CPU/GPU split with transparent fallback:** CPU paths are canonical; GPU helpers wrap the same strategies and take over only when enabled/thresholded, otherwise falling back without diverging user logic.

## Module Structure
```
src/
├── core/
│   ├── config/        # FullConfig types, loader, overrides, schema hooks
│   ├── integrator/    # SimulationEngine, DomainManager, OutputManager, strategies/
│   ├── physics/       # Forces, collisions, reactions, contexts
│   ├── gpu/           # GPUContext, helpers, CUDA kernels (conditional)
│   ├── types/         # IonState, IonEnsemble (SoA), Vec3, collision types
│   ├── log/           # Logger wrappers and sinks
│   ├── io/            # Output helpers (HDF5, writers)
│   ├── param/         # Legacy/bridge parameter helpers
│   └── utils/         # Math, safety guards, profiling hooks used by core
├── fieldsolver/       # Field computation (prototype/future)
├── main/              # CLI entrypoint and setup wiring
├── optimizer/         # Parameter optimization (future-facing)
└── utils/             # Shared utilities (logging, math, profiling, CLI, common helpers)

tests/                 # Unit/integration/physics/gpu suites
validation/            # Reproducible validation configs, scripts, reports
schema/                # JSON Schemas for configs
tools/                 # Developer tools and scripts
```

### Core Dependencies

- **Required**: C++17, CMake 3.16+, HDF5, jsoncpp, OpenSSL, Threads, BLAS
- **Optional**: CUDA 11.0+ (GPU acceleration), OpenMP 4.5+ (threading)
- **External (fetched if missing)**: cxxopts (CLI), spdlog (logging)
- **Header-only**: nlohmann/json (CLI snapshot + fieldsolver helpers)

## Core Architecture

### Simulation Engine

`SimulationEngine` is the central orchestrator. The engine runs on the SoA `IonEnsemble`; any legacy AoS callers are converted once at the boundary before entering the timestep loop.

**Separation of concerns via strategies/factories:** SimulationEngine orchestrates, but domain logic is delegated to pluggable components:
- Geometry/fields: `IDomainGeometry` + `IFieldModel` (wired by DomainManager)
- Forces: `IForce` implementations aggregated by `ForceRegistry`
- Collisions: `ICollisionHandler` selected by factory
- Reactions: `IReactionHandler` selected by factory (CPU stochastic handler or GPU wrapper, depending on config)
- Integrators: `IIntegrationStrategy` selected by factory
Swapping a strategy/factory changes behavior without touching the loop.

1. Initialize `DomainManager` and `OutputManager`, log metadata (AoS init uses a temporary conversion in the SoA path).
2. GPU-capable strategies/handlers are constructed during setup when `simulation.enable_gpu` is true; boundary helpers remain unused by the main loop.
3. Per-timestep loop:
   - Apply ion birth timing (SoA birth flags).
   - Per-ion processing (OpenMP-capable): domain lookup, domain property updates, collisions (if handler), reactions (if enabled), integrator step, boundary checks/domain transitions, time update, safety checks.
   - Output write every `write_interval` steps.
4. Finalize output and optional numerical safety report; GPU stats if used.

**Note:** GPU boundary checks have helpers but are not wired into the main loop yet; space charge uses `SpaceChargeGPUModel` when enabled and available.

**Key Files:**
- `src/core/integrator/SimulationEngine.{h,cpp}` - Main simulation loop
- `src/core/integrator/DomainManager.{h,cpp}` - Multi-domain coordination via geometry/field strategies

### Force System

**Architecture Pattern:** Registry + Strategy (SoA-only)

Forces implement `IForce` and register with `ForceRegistry`. There is a single SoA interface (`compute(const IonEnsemble&, size_t, t, ForceContext)`). AoS hooks were removed from the hot path; tests wrap single ions into a scratch SoA when needed.

**Built-in Forces:**
- `ElectricFieldForce` - E-field via `IFieldModel` (analytical or grid-backed); SSOT: PhysicsSetup injects `FieldProviderModel` when grid data exist, otherwise `AnalyticalFieldModel` fallback
- `MagneticFieldForce` - B-field (Boris integrator compatible)
- `SpaceCharge models` - `ForceRegistry` feeds `ISpaceChargeModel` implementations (Direct/Grid/GPU) produced by `SpaceChargeModelFactory`
- `DampingForce` - Drag depending on chosen deterministic collision model

**Space-Charge Pipeline**
- `ISpaceChargeModel` – interface implemented by:
  - `SpaceChargeDirectModel` (per-domain, O(N²))
  - `SpaceChargeGridModel` (geometry-aware Poisson solver)
  - `SpaceChargeGPUModel` (experimental GPU P³M; enabled via `physics.enable_space_charge_gpu`)
- `SpaceChargeModelFactory` selects the model per-domain (GPU if available, otherwise Direct for small N, else Grid) and `ForceRegistry` injects the Coulomb field into the SoA loop without AoS fallbacks.

**Key Files:**
- `src/core/physics/forces/ForceRegistry.{h,cpp}` - Force management
- `src/core/physics/forces/*Force.{h,cpp}` - Force implementations

### Integrator System

**Architecture Pattern:** Strategy + Factory

Integration strategies implement `IIntegrationStrategy::step(...)` on the SoA container and may optionally expose a `step_batch(...)` hook. GPU-enabled strategies wrap the CPU implementation and dispatch batches when the entire ion set for a timestep satisfies domain/E-field constraints.

**Available Strategies:**
- `RK4Strategy` – 4th order Runge-Kutta (fixed step)
- `RK45Strategy` - Dormand-Prince (adaptive, per-ion dt/state; OpenMP disabled for adaptive RK45 in SimulationEngine)
- `BorisStrategy` – Boris pusher
- `GPUIntegrationStrategy` – wrapper that delegates to RK4/RK45/Boris and hands off batches to the CUDA helper when possible (auto CPU fallback)

**Parallelism note:** RK45 maintains per-ion adaptive state; SimulationEngine disables OpenMP for adaptive RK45. Batch paths (CPU/GPU) are only used when all active ions share the same `dt`.

**Key Files:**
- `src/core/integrator/strategies/IntegrationStrategyFactory.h`
- `src/core/integrator/strategies/GPUIntegrationStrategy.{h,cpp}`
- `src/core/integrator/strategies/*Strategy.{h,cpp}` (CPU implementations)

### GPU Acceleration

**Architecture Pattern:** Hybrid CPU/GPU with automatic dispatch

GPU acceleration uses threshold-based dispatch (default integration/collisions: 5000 ions). Batch GPU integration requires uniform `dt`, a single active domain, and a supported force mix (no space charge, magnetic, or damping unless explicitly enabled in the GPU strategy; GPU damping is not wired from config yet). It expects a grid-backed E-field; non-grid or snapshot-based providers yield zero-field GPU integration, so avoid GPU dispatch in those cases. Collisions/reactions use batch hooks only when `dt` is uniform (per-domain for collisions, global for reactions) and fall back on failure or low counts.

**GPU Features (current state):**
- Integration: `GPUIntegrationStrategy` uses `GPUIntegrationHelper` to run RK4/RK45/Boris batches when one domain + `ElectricFieldForce` is active and no unsupported forces are present; `GridFieldProvider` + `FieldArray` is required for non-zero fields (non-grid or snapshot providers yield zero-field integration)
- Collisions: `GPUCollisionHandler` wraps EHSS/HSS CPU models, advertises `supports_batch()`, and `SimulationEngine::perform_collisions()` groups ions per domain before invoking the GPU helper with automatic CPU fallback. Multi-gas mixtures are supported (up to 8 components per domain; additional entries are truncated with a warning)
- Reactions: `GPUReactionHandler` wraps the stochastic CPU handler and invokes `GPUReactionBackend` batches when enabled and `dt` is uniform across active ions; otherwise it falls back to CPU
- Space charge: P³M helper exposed through `SpaceChargeGPUModel`; enabled when `physics.enable_space_charge_gpu=true` and `ICARION_USE_GPU` is defined (falls back to Grid/Direct otherwise)
- Boundary checks: Helper exists for absorption/cylindrical only, not wired into the loop
- Automatic CPU fallback on errors or below-threshold counts

- **Key Files:**
  - `src/core/integrator/strategies/GPUIntegrationStrategy.{h,cpp}` - Integration wrapper
  - `src/core/gpu/core/GPUIntegrationHelper.{h,cpp}` - CUDA kernels/dispatch
  - `src/core/gpu/core/GPUContext.{h,cpp}` - CUDA context management
  - `src/core/gpu/collisions/GPUCollisionHelper.{h,cpp}`
  - `src/core/gpu/spacecharge/GPUSpaceChargeP3M.{h,cpp}`
  - `src/core/gpu/*.cu`

### Configuration System

**Architecture Pattern:** Schema-driven JSON with validation

All configurations use JSON Schema for validation:

```bash
# Validate any config (from repo root)
python3 schema/validator.py my_config.json
```

**Schema Files:**
- `schema/icarion-config.schema.json` - Master schema
- `schema/simulation.schema.json` - Simulation parameters
- `schema/physics.schema.json` - Physics models
- `schema/domain.schema.json` - Geometry and fields

Find all allowed schemas in the schema/ folder!

**Key Files:**
- `src/core/config/loader/ConfigLoader.{h,cpp}` - JSON loading
- `src/core/config/types/*Config.h` - Configuration structures

## Design Patterns

### Plugin Architecture

**Forces, Collisions, Integrators** use plugin pattern for extensibility:

1. **Define Interface** (`IForce`, `ICollisionHandler`, etc.)
2. **Implement Plugin** (inherit from interface)
3. **Register with Factory** (explicit registration in factory switch/registry)
4. **Configure via JSON** (string-based selection)

### RAII and Smart Pointers

- **Memory Management**: Automatic cleanup via RAII
- **GPU Memory**: CUDA memory wrapped in RAII classes
- **File I/O**: HDF5 files auto-close via destructors

### Single Source of Truth (SSOT)

- **Configuration**: JSON files are authoritative for runtime settings; defaults live in config structs
- **Physics**: Formula implementations in single location
- **Documentation**: Schema files define validation rules and are referenced by docs

## Performance Architecture

### CPU Optimization

- **OpenMP**: Parallel force computation for large ion counts
- **SIMD**: Vectorized operations where applicable
- **Memory Layout**: Structure-of-Arrays for cache efficiency

### GPU Optimization

- **Threshold Dispatch**: Automatic CPU/GPU selection based on problem size
- **Coalesced Memory**: GPU memory access patterns optimized
- **Texture Memory**: Field array interpolation via texture cache

### I/O Optimization

- **HDF5 Compression**: Automatic compression for trajectory data
- **Chunked Writing**: Efficient large dataset handling

## Testing Strategy

### Test Organization
```
tests/
├── unit/           # Isolated component tests (forces, integrators, registries)
├── integration/    # Multi-component workflows (main.cpp path, config, I/O)
├── integrator/     # SimulationEngine-focused cases
├── physics/        # Physics accuracy/edge cases
├── collision/      # Collision-model tests
├── gpu/            # GPU/CPU parity and performance probes
├── instruments/    # Instrument-specific behaviors
├── io/, utils/, helpers/  # I/O, shared helpers, fixtures
└── config/         # Config parsing/override/schema tests
```

Validation suite (scientific reproducibility) lives in `validation/`:
- `configs/`, `scripts/`, `results/`, `figures/`, `logs/`, plus reports like `VALIDATION_REPORT_v1.0.md`.

### Validation Levels

1. **Unit Tests**: Individual functions and classes
2. **Integration Tests**: Module interactions
3. **Physics Tests**: Scientific accuracy (energy conservation, etc.)
4. **Parity Tests**: CPU/GPU result matching
5. **End-to-End Tests**: Complete simulation workflows

## Implementation Details

For detailed implementation specifics, see:

- **Force Implementation**: Source files in `src/core/physics/forces/`
- **GPU Kernels**: [GPU Architecture Guide](GPU_ARCHITECTURE.md) and `src/core/gpu/`
- **Configuration Schema**: Schema files in `schema/` directory
- **Testing Guide**: Test files in `tests/` directory
- **API Reference**: [Developer Guide](DEVELOPERS_GUIDE.md)

---

**Last Updated:** December 2025 (v1.0)

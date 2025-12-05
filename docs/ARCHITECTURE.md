# ICARION Architecture Overview

**Version:** v1.0
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

- **Required**: C++17, CMake 3.16+, HDF5, spdlog
- **Optional**: CUDA 11.0+ (GPU acceleration), OpenMP (threading)
- **External**: cxxopts (CLI), nlohmann/json (config parsing)

## Core Architecture

### Simulation Engine

`SimulationEngine` is the central orchestrator. The engine runs on the SoA `IonEnsemble`; any legacy AoS callers are converted once at the boundary before entering the timestep loop.

**Separation of concerns via strategies/factories:** SimulationEngine orchestrates, but domain logic is delegated to pluggable components:
- Geometry/fields: `IDomainGeometry` + `IFieldModel` (wired by DomainManager)
- Forces: `IForce` implementations aggregated by `ForceRegistry`
- Collisions: `ICollisionHandler` selected by factory
- Reactions: `IReactionHandler` selected by factory
- Integrators: `IIntegrationStrategy` selected by factory
Swapping a strategy/factory changes behavior without touching the loop.

1. Initialize `DomainManager` and `OutputManager`, log metadata (AoS init uses a temporary conversion in the SoA path).
2. Optionally initialize GPU helpers (integration/collisions; space charge/boundary helpers exist but are not yet dispatched).
3. Per-timestep loop:
   - Apply ion birth timing (SoA birth flags).
   - Per-ion processing (OpenMP-capable): domain lookup, domain property updates, collisions (if handler), reactions (if enabled), integrator step, boundary checks/domain transitions, time update, safety checks.
   - Output write every `write_interval` steps.
4. Finalize output and optional numerical safety report; GPU stats if used.

**Note:** GPU space-charge and GPU boundary checks have helpers but are not wired into the main loop yet; CPU paths run instead.

**Key Files:**
- `src/core/integrator/SimulationEngine.{h,cpp}` - Main simulation loop
- `src/core/integrator/DomainManager.{h,cpp}` - Multi-domain coordination via geometry/field strategies

### Force System

**Architecture Pattern:** Registry + Strategy (SoA-first)

Forces implement `IForce` and register with `ForceRegistry`. The SoA path is primary (`compute_soa(...)`), with `compute(...)` retained for tests and legacy hooks.

**Built-in Forces:**
- `ElectricFieldForce` - E-field via `IFieldModel` (analytical or grid-backed); SSOT: PhysicsSetup injects `FieldProviderModel` when grid data exist, otherwise `AnalyticalFieldModel` fallback
- `MagneticFieldForce` - B-field (Boris integrator compatible)
- `SpaceChargeDirect` - Ion-ion Coulomb interactions
- `DampingForce` - Drag depending on chosen deterministic collision model

**Key Files:**
- `src/core/physics/forces/ForceRegistry.{h,cpp}` - Force management
- `src/core/physics/forces/*Force.{h,cpp}` - Force implementations

### Integrator System

**Architecture Pattern:** Strategy + Factory

Integration strategies implement `IIntegrationStrategy::step(IonEnsemble&, size_t ion_idx, ...)` on the SoA container.

**Available Strategies:**
- `RK4Strategy` - 4th order Runge-Kutta (default, stable, fixed step size)
- `RK45Strategy` - Runge-Kutta-Fehlberg (adaptive step size) with Dormand-Prince coefficients
- `BorisStrategy` - Boris pusher (for strong magnetic fields)

**Key Files:**
- `src/core/integrator/strategies/IntegrationStrategyFactory.h` - Factory
- `src/core/integrator/strategies/*Strategy.{h,cpp}` - Implementations

### GPU Acceleration

**Architecture Pattern:** Hybrid CPU/GPU with automatic dispatch

GPU acceleration uses threshold-based dispatch (default integration/collisions: 5000 ions):
- N < threshold → CPU path (OpenMP-capable)
- N ≥ threshold → GPU integration helper (RK4/RK45/Boris) and optional GPU collision helper (HSS/EHSS).

**GPU Features (current state):**
- Integration: RK4/RK45/Boris batch kernels
- Collisions: HSS/EHSS batch helper with CPU fallback
- Space charge: P³M helper exists but not called from the main loop yet
- Boundary checks: Helper exists for absorption/cylindrical only, not wired into the loop
- Automatic CPU fallback on errors or below-threshold counts

**Key Files:**
- `src/core/gpu/core/GPUContext.{h,cpp}` - CUDA context management
- `src/core/gpu/core/GPUIntegrationHelper.{h,cpp}` - GPU integration dispatch
- `src/core/gpu/collisions/GPUCollisionHelper.{h,cpp}` - GPU collision processing
- `src/core/gpu/spacecharge/GPUSpaceChargeP3M.{h,cpp}` - GPU space charge helper
- `src/core/gpu/*.cu` - CUDA kernels

### Configuration System

**Architecture Pattern:** Schema-driven JSON with validation

All configurations use JSON Schema for validation:

```bash
# Validate any config
python3 schema/validator.py schema/icarion-config.schema.json my_config.json
```

**Schema Files:**
- `schema/icarion-config.schema.json` - Master schema
- `schema/simulation.schema.json` - Simulation parameters
- `schema/physics.schema.json` - Physics models
- `schema/domain.schema.json` - Geometry and fields

Find all allowed schemes in the schema/ folder!

**Key Files:**
- `src/core/config/loader/ConfigLoader.{h,cpp}` - JSON loading
- `src/core/config/types/*Config.h` - Configuration structures

## Design Patterns

### Plugin Architecture

**Forces, Collisions, Integrators** use plugin pattern for extensibility:

1. **Define Interface** (`IForce`, `ICollisionHandler`, etc.)
2. **Implement Plugin** (inherit from interface)
3. **Register with Factory** (automatic discovery)
4. **Configure via JSON** (string-based selection)

### RAII and Smart Pointers

- **Memory Management**: Automatic cleanup via RAII
- **GPU Memory**: CUDA memory wrapped in RAII classes
- **File I/O**: HDF5 files auto-close via destructors

### Single Source of Truth (SSOT)

- **Configuration**: JSON files are authoritative (not hardcoded defaults)
- **Physics**: Formula implementations in single location
- **Documentation**: Schema files generate validation and docs

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
- **Metadata Caching**: Reduce repeated schema validation

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

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
7. [Future Architecture](#future-architecture)

## System Overview

ICARION is a modular C++17 framework for ion trajectory simulation in mass spectrometry and ion mobility. The architecture follows a **plugin-based design** with clear separation of concerns:

- **Core**: Physics simulation engine (forces, integration, collisions)
- **Config**: JSON-based configuration system with schema validation
- **I/O**: HDF5-based input/output with structured data formats
- **GPU**: Optional CUDA acceleration with automatic CPU fallback
- **Tools**: CLI interface, configuration helpers, validation tools

### Key Architectural Principles

1. **Modularity**: Components are loosely coupled via interfaces
2. **Extensibility**: New forces, collision models, integrators via plugin pattern
3. **Performance**: CPU/GPU hybrid with smart threshold dispatch
4. **Reliability**: Comprehensive test coverage, schema validation, error handling
5. **Usability**: Intuitive JSON configuration, extensive documentation

## Module Structure

```
src/
├── core/              # Core simulation engine
│   ├── physics/       # Forces, collisions, reactions
│   ├── integrator/    # Time integration strategies
│   ├── config/        # Configuration loading and validation
│   └── gpu/          # CUDA acceleration
├── fieldsolver/       # Electric/magnetic field computation
├── utils/             # Shared utilities (logging, math, I/O)
├── optimizer/         # Parameter optimization (future)
└── main/             # CLI application and setup
```

### Core Dependencies

- **Required**: C++17, CMake 3.16+, HDF5, spdlog
- **Optional**: CUDA 11.0+ (GPU acceleration), OpenMP (threading)
- **External**: cxxopts (CLI), nlohmann/json (config parsing)

## Core Architecture

### Simulation Engine

`SimulationEngine` is the central orchestrator:

```cpp
// Main simulation loop
for (double t = 0; t < total_time; t += dt) {
    // 1. Apply forces and integrate (CPU or GPU)
    integrator->integrate_timestep(ions, dt, t);
    
    // 2. Handle collisions with buffer gas
    collision_handler->process_collisions(ions, dt);
    
    // 3. Handle chemical reactions
    reaction_handler->process_reactions(ions, dt);
    
    // 4. Apply boundary conditions
    boundary_manager->check_boundaries(ions);
    
    // 5. Write output (every N steps)
    if (step % write_interval == 0) {
        output_manager->write_trajectories(ions, t);
    }
}
```

**Key Files:**
- `src/core/integrator/SimulationEngine.{h,cpp}` - Main simulation loop
- `src/core/integrator/DomainManager.{h,cpp}` - Multi-domain coordination

### Force System

**Architecture Pattern:** Registry + Strategy

Forces implement `IForce` interface and register with `ForceRegistry`:

```cpp
class IForce {
public:
    virtual Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) = 0;
};
```

**Built-in Forces:**
- `ElectricFieldForce` - E-field from analytical or HDF5 sources
- `MagneticFieldForce` - B-field (Boris integrator compatible)
- `SpaceChargeForce` - Ion-ion Coulomb interactions
- `DampingForce` - Velocity-dependent drag

**Key Files:**
- `src/core/physics/forces/ForceRegistry.{h,cpp}` - Force management
- `src/core/physics/forces/*Force.{h,cpp}` - Force implementations

### Integrator System

**Architecture Pattern:** Strategy + Factory

Integration strategies implement `IIntegrationStrategy`:

```cpp
class IIntegrationStrategy {
public:
    virtual void integrate_single_step(IonState& ion, double dt, double t, 
                                     const ForceContext& ctx) = 0;
};
```

**Available Strategies:**
- `RK4Strategy` - 4th order Runge-Kutta (default, stable)
- `RK45Strategy` - Runge-Kutta-Fehlberg (adaptive step size)
- `BorisStrategy` - Boris push (for strong magnetic fields)

**Key Files:**
- `src/core/integrator/strategies/IntegrationStrategyFactory.h` - Factory
- `src/core/integrator/strategies/*Strategy.{h,cpp}` - Implementations

### GPU Acceleration

**Architecture Pattern:** Hybrid CPU/GPU with automatic dispatch

GPU acceleration uses threshold-based dispatch (default: 5000 ions):
- N < threshold → CPU integration (lower overhead)
- N ≥ threshold → GPU integration (higher throughput)

**GPU Features:**
- All integrators (RK4, RK45, Boris)
- Collision models (HSS, EHSS)
- Space charge forces (planned)
- Automatic fallback on errors

**Key Files:**
- `src/core/gpu/GPUContext.{h,cpp}` - CUDA context management
- `src/core/gpu/GPUIntegrationHelper.{h,cpp}` - GPU integration dispatch
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
├── unit/           # Isolated component tests
├── integration/    # Multi-component interaction tests
├── gpu/           # GPU-specific tests and CPU/GPU parity
├── physics/       # Physics accuracy validation
└── helpers/       # Test utilities and fixtures
```

### Validation Levels

1. **Unit Tests**: Individual functions and classes
2. **Integration Tests**: Module interactions
3. **Physics Tests**: Scientific accuracy (energy conservation, etc.)
4. **Parity Tests**: CPU/GPU result matching
5. **End-to-End Tests**: Complete simulation workflows

### Key Test Files

- `tests/integrator/test_rk45_boris_parity.cpp` - Integration accuracy
- `tests/gpu/test_gpu_*.cpp` - GPU functionality and parity
- `tests/physics/test_*_physics.cpp` - Physics model validation

## Future Architecture

### Planned Extensions

1. **Advanced GPU Features**: Space charge on GPU, field solver acceleration
2. **Distributed Computing**: MPI support for large-scale simulations
3. **Python Interface**: Python bindings for scripting and analysis
4. **Real-time Visualization**: Live trajectory plotting during simulation
5. **Machine Learning**: ML-assisted parameter optimization

### Architectural Evolution

- **Plugin System**: More runtime-discoverable plugins
- **Microservice Architecture**: Separate field solver, reaction engine services
- **Cloud Integration**: Kubernetes deployment, cloud storage backends

## Implementation Details

For detailed implementation specifics, see:

- **Force Implementation**: [Force System Documentation](src/core/physics/forces/README.md)
- **GPU Kernels**: [GPU Implementation Guide](src/core/gpu/README.md) 
- **Configuration Schema**: [Schema Documentation](schema/README.md)
- **Testing Guide**: [Testing Documentation](tests/README.md)
- **API Reference**: [Developer Guide](DEVELOPERS_GUIDE.md)

---

**Last Updated:** December 2025 (v1.0)
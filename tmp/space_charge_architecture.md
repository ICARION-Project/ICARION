# Space-Charge Architecture Proposal (vNext)

## Goals
- Replace ad-hoc direct/grid/GPU branching in `PhysicsSetup` with strategy-based `ISpaceChargeModel`.
- Leverage `IDomainGeometry`/`IFieldModel` for domain-aware grids/boundary conditions instead of global bounding boxes.
- Enable future geometry-specific solvers (e.g., Orbitrap, cylindrical) and easier GPU integration.

## Proposed Components

### 1. ISpaceChargeModel Interface
```cpp
class ISpaceChargeModel {
public:
    virtual ~ISpaceChargeModel() = default;
    virtual void update_fields(const core::IonEnsemble& ions, double time_s) = 0;
    virtual core::Vec3 sample_electric_field(size_t ion_idx) const = 0;
    virtual std::string name() const = 0;
};
```
- `update_fields` computes the field for current ion distribution (deposit + solve).
- `sample_electric_field` returns `E` at ion `ion_idx` (SoA index).
- Optional helpers: `void on_step_begin(double dt)` for caching, `double energy()` for diagnostics.

### 2. Concrete Models
- `SpaceChargeDirectModel`: current `SpaceChargeDirect` O(N^2) interaction.
- `SpaceChargeGridModel`: wraps `Grid3D + PoissonSolver + depositCharge`. Gets grid/BC info from `IDomainGeometry`.
- `SpaceChargeGPUModel`: wraps `gpu::GPUSpaceChargeP3M`.
- Future: geometry-specific (e.g., `OrbitrapSpaceChargeModel`) using analytical BCs.

### 3. Domain/Geometry Hooks
Extend `IDomainGeometry` (or companion) with:
- `BoundingBox bounding_box(double margin)` → axis-aligned (global coords) sized to domain, not ions.
- `void apply_spacecharge_dirichlet(Grid3D&, std::vector<char>& mask, std::vector<double>& values)` to specify electrode potentials.
- `bool supports_spacecharge_grid()`, `SpaceChargeGridPreference preference()` to guide factory.

### 4. SpaceChargeModelFactory
```cpp
class SpaceChargeModelFactory {
public:
    static std::unique_ptr<ISpaceChargeModel> create(
        const config::DomainConfig& domain,
        size_t ion_count,
        const config::FullConfig& config,
        gpu::GPUContext* gpu_ctx);
};
```
Decision tree:
1. Domain disables space charge → return nullptr.
2. Ion count < threshold → `SpaceChargeDirectModel`.
3. Ion count >= threshold:
   - If GPU enabled and domain preference allows → `SpaceChargeGPUModel`.
   - Else if geometry supports grid → `SpaceChargeGridModel` with geometry-defined bbox & BCs.
   - Else fallback (e.g., `DirectModel`).

### 5. ForceRegistry / SimulationEngine integration
- `ForceRegistry` holds a `std::unique_ptr<ISpaceChargeModel>` and exposes a `SpaceChargeForce` wrapper implementing `IForce`.
- `SpaceChargeForce::compute_soa` fetches `E = model.sample_electric_field(i)` and returns `q * E`.
- SimulationEngine lifecycle per timestep:
  1. Before per-ion loop: `space_charge_model->update_fields(ensemble, t)` (if model present).
  2. Integrator calls ForceRegistry; space-charge force uses latest field.

### 6. Testing Plan
- Reuse existing tests (`test_space_charge_force`, `test_space_charge_integration`, `test_poisson_solver`).
- Add `test_space_charge_model_factory` (decision tree), `test_space_charge_model_grid` (geometry-based bbox), GPU parity test stub.

### 7. Migration Steps
1. Introduce `ISpaceChargeModel`, implement adapters for existing Direct/Grid/GPU logic.
2. Update `PhysicsSetup` to instantiate per-domain models via factory.
3. Modify `ForceRegistry` to accept `set_space_charge_model(...)` and add `SpaceChargeForce` automatically when non-null.
4. Update SimulationEngine to call `update_fields` each timestep.
5. Extend `IDomainGeometry` with bbox/BC hooks (Orbitrap, Cylindrical first).
6. Add docs/tests.

## Benefits
- Clear separation (strategy pattern) of space-charge implementations.
- Geometry-aware grids/BCs → improved physical fidelity.
- Easier to add future models (e.g., instrument-specific solvers, GPU improvements).

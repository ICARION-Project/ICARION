# Space-Charge Solver Refactor Plan

## Goals
1. Introduce `ISpaceChargeModel` abstraction with Direct/Grid/GPU implementations.
2. Use `IDomainGeometry` (and `IFieldModel`) to provide bounding boxes & boundary conditions for grid-based solvers.
3. Simplify solver selection logic, document fallback strategy, and provide regression tests.

## Phase 1 – Interfaces & Direct Model (Day 1)
- [x] Define `ISpaceChargeModel` interface (update_fields/sample_E/name).
- [x] Implement `SpaceChargeDirectModel` by adapting existing `SpaceChargeDirect` logic.
- [ ] Add `SpaceChargeForce` adapter implementing `IForce` using an `ISpaceChargeModel`. *(Re-evaluated: direct accumulation inside `ForceRegistry` – optional helper still pending.)*
- [x] Update `ForceRegistry` to own optional space-charge model + integrate contribution in compute path.
- [x] SimulationEngine change: call `space_charge_model->update_fields(ensemble, t)` once per timestep before force evaluations.
- [x] Tests: `test_space_charge_model_direct` exercises model + ForceRegistry integration. (Legacy force tests retained for now.)

## Phase 2 – Grid Model & Geometry Hooks (Day 2-3)
- [ ] Extend `IDomainGeometry` with:
  - [x] `BoundingBox global_bounding_box(double margin) const`.
  - [ ] `void apply_spacecharge_dirichlet(Grid3D&, std::vector<char>& mask, std::vector<double>& values) const` (default no-op).
- [ ] Implement `SpaceChargeGridModel`:
  - Construct grid from geometry bbox (configurable padding/cell size).
  - Deposit charge (`depositCharge.cpp`), call `PoissonSolver` using geometry-supplied BCs.
  - Cache `E` field and serve through `sample_E`.
- [ ] Update `PoissonSolver` doc/comments to clarify solver selection (keep GS + Multigrid, deprecate extra modes for now).
- [x] Modify `PhysicsSetup::add_space_charge_forces` to instantiate per-domain models via new `SpaceChargeModelFactory` (direct path active; legacy grid/GPU branch remains as fallback until GridModel lands).
- [ ] Tests: expand `test_space_charge_integration` to cover geometry-based bbox; add `test_cylindrical_spacecharge_bbox` (unit-level).

## Phase 3 – GPU Model Integration (Day 4)
- [ ] Wrap `gpu::GPUSpaceChargeP3M` in `SpaceChargeGPUModel` (implement interface).
- [ ] Extend factory to try GPU when `config.physics.enable_space_charge_gpu` and thresholds satisfied.
- [ ] Ensure CPU fallback path (if GPU fails) returns Grid/Direct.
- [ ] Tests: basic GPU stub test or compile-time skip; keep existing GPU helper tests for regression.

## Phase 4 – Docs & Validation (Day 5)
- [ ] Update `docs/ARCHITECTURE.md` + `docs/DEVELOPERS_GUIDE.md` to describe `ISpaceChargeModel` design.
- [ ] Document solver selection + geometry-driven grids in `docs/CONFIG_GUIDE.md` (space charge section).
- [ ] Run validation scripts (IMS drift, mixture thermalization) to confirm behavior.

## Risks & Mitigations
- Lifecycle: ensure `update_fields` invoked before `compute_soa` reads (unit test with instrumentation).
- Performance: geometry bbox may be large → allow overrides or heuristics (e.g., clamp to ion extents if geometry is huge).
- Boundary conditions: start with simple BC (Dirichlet=0) for geometries lacking data; flag TODO for electrode potentials.
- GPU: keep optional; failure must degrade gracefully.

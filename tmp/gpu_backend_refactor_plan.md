# GPU Backend Refactor Plan

## Motivation
Current GPU features (integrators, collisions, space charge) are partly wired through ad-hoc helpers and threshold checks inside `SimulationEngine`. Space charge now uses a clean factory/strategy pattern (`ISpaceChargeModel`). To make GPU support first-class and consistent we need factories/adapters for the remaining modules.

## Goals
1. Introduce polymorphic GPU-aware backends for integration strategies and collision handlers (and optionally reactions/output as follow-up).
2. Centralize GPU selection logic in factories (SSOT: configuration + thresholds + device availability).
3. Maintain automatic CPU fallbacks and deterministic logging of chosen backends.
4. Keep GPU kernels untouched initially; focus on architecture + wiring + tests.

## Phase 0 – Inventory & Requirements
- [ ] Document current GPU entry points:
  - `GPUIntegrationHelper` (RK4/RK45/Boris)
  - `GPUCollisionHelper` (HSS/EHSS)
  - `GPUContext`, device detection, thresholds in `SimulationEngine`.
- [ ] Confirm configuration knobs (`simulation.enable_gpu`, potential per-module flags).
- [ ] Decide default thresholds and logging expectations (reuse current numbers).

## Phase 1 – Integration Strategy Factory Upgrade
- [ ] Extend `IIntegrationStrategy` to allow GPU-backed implementations (option A: derive `GPUEnabledStrategy` that wraps helper; option B: inject `IIntegratorBackend`).
- [ ] Implement GPU wrappers:
  - [ ] `RK4StrategyGPU`
  - [ ] `RK45StrategyGPU`
  - [ ] `BorisStrategyGPU`
  (Each should internally reference `GPUIntegrationHelper`, handle allocation, and fall back to CPU if helper fails mid-run.)
- [ ] Update `IntegrationStrategyFactory`:
  - [ ] Evaluate config + ion count + GPU availability.
  - [ ] Return GPU strategy when enabled and thresholds met; otherwise CPU strategy.
  - [ ] Emit log lines summarizing the choice (mirror SpaceCharge factory style).
- [ ] SimulationEngine changes:
  - [ ] Remove direct calls to `GPUIntegrationHelper`; rely solely on selected strategy.
  - [ ] Ensure existing SoA tests still pass (parity between CPU/GPU strategies).
- [ ] Tests:
  - [ ] Compile-time test to ensure GPU strategies build (stub when CUDA disabled).
  - [ ] Parity test reusing existing `test_rk45_boris_parity` but parameterized by backend.
  - [ ] Update `tests/README.md` to reflect strategy selection.

## Phase 2 – Collision Handler Factory Upgrade
- [ ] Introduce `ICollisionBackend` or GPU-aware handler base class.
- [ ] Create `GPUCollisionHandler` adapter that uses `GPUCollisionHelper`.
- [ ] Extend `CollisionHandlerFactory::create(...)`:
  - [ ] Accept `FullConfig` / GPU preference.
  - [ ] If GPU enabled + threshold met + handler supports GPU (HSS/EHSS), return GPU variant.
  - [ ] Provide logging describing chosen backend / fallback reason.
- [ ] SimulationEngine should only interact with `ICollisionHandler`; remove bespoke GPU threshold checks.
- [ ] Tests:
  - [ ] CPU/GPU parity tests for collision statistics (reuse existing GPU thermalization tests but run through new handler).
  - [ ] CPU-only build stub test ensuring factory returns CPU handler gracefully.
  - [ ] Documentation updates (Architecture + Developers Guide).

## Phase 3 – Field Interpolation (Stretch)
- [ ] Evaluate extending `IFieldModel` / `IFieldProvider` with GPU-aware backends so grid-based field maps can be uploaded once and interpolated on-device.
- [ ] Update `PhysicsSetup` field loading to optionally create `GPUFieldProvider` objects when GPU integration is active (detect via config + `enable_gpu`).
- [ ] Ensure `ElectricFieldForce` remains agnostic (receives a field model that can answer `E(pos, t)` even when running inside GPU strategies).
- [ ] Tests: parity between CPU and GPU interpolation for a sample field array; stub test for CPU-only builds.

## Phase 4 – Optional: Reaction Handler & Output
- [ ] Evaluate if reactions need GPU analog (currently CPU only; might remain so until demand arises).
- [ ] Consider output manager GPU path (not required for v1.0; keep on roadmap).

## Phase 4 – Documentation & Tooling
- [ ] Update docs/ARCHITECTURE.md and docs/DEVELOPERS_GUIDE.md with new factory diagrams (Integration, Collision).
- [ ] Update CLI/CONFIG guide if new flags introduced (e.g., per-module enablement).
- [ ] Ensure tests/README.md references new GPU-oriented CTests.
- [ ] Add release notes entry describing unified GPU backend architecture.

## Risks & Mitigations
- **Device availability mid-run:** GPU strategy should detect failures and degrade to CPU with logging; unit tests must cover fallback path.
- **State synchronization:** Strategies need to ensure IonEnsemble SoA stays in sync when GPU buffer updates occur (reuse existing helper logic).
- **Build configurations:** CPU-only builds must compile (provide stub classes) and tests should confirm no GPU code is invoked.
- **Performance regression:** After refactor, benchmark to confirm throughput matches prior helper-based wiring.

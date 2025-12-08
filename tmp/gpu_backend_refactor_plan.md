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

## Phase 4 – Reaction Handler Parity (new)

### Requirements *(status: scaffolded in code)*
- Preserve SSOT: reaction logic reads `ReactionDatabase`, `SpeciesDatabase`, and `EnvironmentConfig` directly, just like the CPU `StochasticReactionHandler`.
- Reuse the existing `IReactionHandler` interface so `SimulationEngine` remains agnostic; we simply plug in a GPU-backed implementation.
- Deterministic logging and fallback: if GPU is unavailable, disabled, or the handler encounters an unsupported feature, control must fall back to the CPU handler with a single warning.

### Architecture
1. **GPUReactionBackend (new helper)**
   - Mirrors `GPUCollisionHelper`: owns device buffers for ion properties relevant to reactions (species ids, charge state, cached coefficients) and manages cuRAND streams.
   - Provides a batch API `process_reactions(core::IonEnsemble&, const ReactionDatabase&, const SpeciesDatabase&, const EnvironmentConfig&, double dt)`.
   - Threshold-aware to avoid GPU launches for tiny ensembles.

2. **GPUReactionHandler (Strategy)**
   - Implements `IReactionHandler`.
   - Internally owns a `GPUReactionBackend` and a CPU fallback pointer.
   - `handle_reaction(...)` simply marks ions for GPU processing; actual batch execution happens in the same integration loop phase as CPU reactions to keep behavior identical.

3. **ReactionHandlerFactory Changes**
   - Accepts `FullConfig`, GPU enable flags, and optional thresholds.
   - Logs the chosen backend (CPU vs GPU) at creation time, mirroring Collision/Integration factories.
   - Returns `GPUReactionHandler` when:
     * `physics.enable_reactions` is true,
     * GPU support is enabled,
     * the selected reaction model is supported on GPU (initially stochastic MC),
     * GPUContext creation succeeds.
   - Otherwise returns the existing CPU handler.

4. **SimulationEngine Flow** *(implemented)* 
   - Unchanged from the caller perspective; engine still holds an `IReactionHandler`.
   - Batch execution order remains: collisions → reactions → forces. GPU handler will upload/download as needed but must expose identical semantics.

### Testing & Docs
- Extend `tests/physics/reactions` with a GPU parity suite (e.g., thermalization or product branching comparisons) gated by `ICARION_USE_GPU`. *(pending real backend)*
- Updated `docs/ARCHITECTURE.md` and `docs/DEVELOPERS_GUIDE.md` with the new wrapper details. `tests/README.md` entry to follow once parity tests exist.
- Add release note entry once the GPU backend does real work.

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

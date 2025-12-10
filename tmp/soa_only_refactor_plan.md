# SoA-only Force/Integrator Refactor Plan

## Scope
- Remove AoS force path entirely: single SoA force interface and registry.
- Update integrators (RK4/RK45/Boris) and SimulationEngine to use only `IonEnsemble` views.
- Fix tests/mocks and documentation references.

## Steps
1) **Interface/Registry**
   - `IForce`: replace `compute` + `compute_batch` with a single SoA method (e.g., `compute(const IonEnsemble&, size_t, double, const ForceContext&)`).
   - `ForceRegistry`: drop AoS overload; SoA-only `compute_total_force` using the new method.
   - ForceContext: ensure fields used are compatible with SoA-only (remove AoS expectations).

2) **Integrators/Engine**
   - RK4/RK45/Boris: rework stage evaluations to operate on `IonEnsemble` (no `IonState` copies). Use SoA helper for accelerations; remove AoS calls.
   - SimulationEngine: ensure force evaluations go through SoA registry; remove AoS branches.
   - GPU parity paths: align to new API.

3) **Tests/Mocks**
   - Update force mocks in integrator/force tests to implement the new single method.
   - Adjust documentation snippets that mention `compute_soa`/AoS.
   - Run/adjust integrator and force-related CTests.

4) **Validation**
   - Build and run targeted tests: `test_force_registry`, `test_rk4_strategy`, `test_rk45_strategy`, `test_boris_strategy`, `test_simulation_engine`, GPU parity tests if applicable.

## Notes
- Expect intermediate breakage while integrators are migrated; keep changes localized per step.
- Consider committing in stages: interface/registry, integrators/engine, tests/docs.

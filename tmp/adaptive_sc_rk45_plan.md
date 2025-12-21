Adaptive RK45 with Stage-Synchronous Space Charge
=================================================

Goal
- Support adaptive Dormand–Prince (RK45) with space-charge fields refreshed at each stage.
- Avoid the current hard prohibition of SC+adaptive.

Current State
- Uniform (fixed-dt) RK4 and a custom fixed-dt DP5 block rebuild SC per stage in SimulationEngine.
- Adaptive RK45Strategy runs its own substeps with a single SC update per macro-step; SC+adaptive throws.

Plan

1) Expose staged computation
- Add an interface to RK45Strategy (or a helper) to run one DP step with a given dt and return k1–k6 (pos/vel deltas) and y4/y5 estimates, without internal error control.
- Alternatively, lift the DP coefficients into a shared helper to avoid code duplication between the existing fixed-dt DP block and the adaptive path.

2) Stage-synchronous SC rebuild
- For each stage (c2..c6), materialize stage positions/velocities for all ions, set current_time_ to t + c_i*dt, call update_space_charge_models(ensemble), then compute accelerations using the rebuilt SC field.
- On step reject (error > tol), rerun the stages with the reduced dt (SC rebuilds must be repeated).

3) Adaptive loop integration
- Implement an adaptive stepper in SimulationEngine (or a SC-aware wrapper around RK45Strategy):
  - Propose dt, run DP stages with SC rebuilds, compute error norm, accept/reject, update dt.
  - Initially use uniform dt across ions (per current RK45 usage); per-ion adaptive dt can remain out of scope.
- Keep existing boundary handling (checked after accepted step) and RNG determinism.

4) Gating
- Enable SC+adaptive only for CPU, single-domain electric/damping setups initially; keep GPU off.
- Maintain the current throw for unsupported cases until this path is stable; add a feature flag to revert to the old behavior if needed.

5) Validation
- Unit test: single ion in a simple SC field (e.g., uniform charge sphere) – compare against analytic or high-accuracy reference.
- Regression: compare adaptive SC vs fixed-step SC on a small cloud; ensure no crashes on step rejects and reasonable performance.

Risks
- Performance hit: SC rebuilds per stage and per retry; adaptive steps may be slower than fixed dt for SC-heavy runs.
- Complexity: touching RK45 internals and SimulationEngine loop; needs thorough testing.

Notes (current implementation)
- Adaptive SC path can be disabled via `ICARION_ADAPTIVE_SC=0` (reverts to legacy throw); default is enabled.
- SC fields are rebuilt at every RK stage on accept/reject attempts; expect substantial cost vs fixed-step RK4.
- Initial validation: `test_simulation_engine_adaptive_sc` covers enabled path and env guard. Further physics validation still required.

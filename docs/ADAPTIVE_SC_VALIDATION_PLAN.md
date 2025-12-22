Adaptive RK45 + Space-Charge Validation/Bench Plan
==================================================

Scope
- Validate adaptive RK45 with stage-synchronous space charge (CPU direct/grid) against a high-accuracy fixed-step RK4 reference.
- Measure performance/scaling of SC rebuilds and rejection behavior.
- Provide guidance for tolerances; keep GPU out (v1.0 GPU is disabled).

1) Physics parity (small cloud)
- Setup: single cylindrical domain, direct SC, N=100 ions, random sphere placement, zero E-field, no collisions/reactions.
- Reference: RK4 with very small dt (e.g., 1e-11 s) for 1–5 ns total.
- Test path: RK45 adaptive with stage-refresh SC, default tolerances.
- Metrics: RMS position/velocity error vs. reference at end, energy drift, acceptance/rejection counts.
- Implement as a reusable validation script (under validation/scripts/physics/) that dumps CSV of RK45 vs RK4 and a short summary.

2) Performance/scaling microbench
- Setup: cylindrical or stacked domains, grid SC (not GPU), N≈1e4 ions.
- Run RK45 adaptive for a few steps with SC enabled; log:
  - time per step, per-stage SC rebuild wall-clock
  - number of accepted/rejected steps
  - min/avg/max dt
- Output: CSV + log snippet to spot scaling issues; abort if per-step wall time exceeds a threshold (e.g., >1s for N=1e4).

3) Rejection-path coverage
- Introduce a solver knob (e.g., accept_at_dt_min=false) to force failure after MAX_REJECT_ATTEMPTS for test-only configs.
- Add a unit test that sets absurd tolerances and a strong SC field, asserts the detailed exception message triggers.

4) Tolerance guidance
- Sweep rtol/atol on the small-cloud case; report error vs. runtime to recommend defaults for SC-heavy runs.
- Document in CONFIG_GUIDE once results are in.

Out of scope
- GPU space charge (runtime-disabled in v1.0).
- Large-scale P3M; stay on CPU Direct/Grid for validation.

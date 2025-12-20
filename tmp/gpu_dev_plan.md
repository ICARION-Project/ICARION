GPU Development Plan (v1.0 safety + enablement)
===============================================

Goals
- v1.0: prevent invalid GPU runs, document limits, and deliver feature parity (E/B + damping + space charge) with caching and tests.

1) Consolidate current state
- Hard-gate GPU when unsupported forces are present: space charge, damping, magnetic, composite/multiple electric forces; log “experimental” with reason. Runtime check in GPUIntegrationStrategy + CLI note.
- Surface GPU disablement in output metadata.

2) Minimal v1.0 safety
- Ensure reproducibility: per-ion RNG on GPU collisions or disable GPU collisions entirely.
- Disable GPU when space charge is enabled; warn clearly.

3) Architecture clean-up (must land for v1.0 if we enable GPU)
- Decouple GPU strategy from “exactly one ElectricFieldForce”: accept a vector of force functors; bail only on unsupported types.
- Persist device field uploads across domain switches; cache by (provider ptr, domain id).
- Unify integrator dispatch: no duplicate selection logic in GPUIntegrationStrategy; reuse SimulationEngine integrator choice.

4) Feature enablement sequence (target v1.0) with concrete deliverables
- Multi-domain E/B: accept multiple ElectricFieldForces and domain offsets; extend GPUIntegrationHelper to take per-ion domain offsets. Tests: CPU vs GPU drift for 2-domain E-only.
- Damping on GPU: implement friction kernel (Teff scaling) on device; share CCS lookup (device table) with CPU data; add kernel path in GPUIntegrationHelper. Tests: heating curve vs CPU, mobility match within tolerance.
- Space charge on GPU: Phase 1 hybrid (CPU SC field compute per stage, GPU integrates E+B+D with injected SC field arrays); Phase 2 GPU P³M kernel. Tests: single-cell SC force vs CPU, stage-synchronous drift with SC.
- Magnetic: allow static B provider (no map uploads) in GPU helper; reject maps until upload path exists. Tests: cyclotron freq vs CPU.
- Boundary actions: wire absorption/cylindrical helper into GPU timestep; fall back to CPU for others.

4a) Kernel/interface sketch (Damping)
- Inputs per ion: pos (x,y,z), vel (vx,vy,vz), mass, charge, CCS lookup (per-species table or per-ion CCS), Teff flag, gas params (T, density, gas mass), dt.
- Force: F_damp = - (m / tau) * v with tau from mobility/CCS; scale mobility by sqrt(T0/Teff) with Teff = T + m_gas/(3*kB)*|v_drift|^2 (reuse CPU formula).
- Output: updated velocities (and optionally positions if fused with E/B), optional energy stats.
- Integration: either fuse into existing RK4/RK45 kernels (E + D) or run a separate damping step per stage.

4b) Kernel/interface sketch (Hybrid space charge)
- CPU computes SC field grids stage-synchronously. GPU ingests SC grid (Ex,Ey,Ez) as extra provider and sums E_total = E_field + E_SC.
- GPUIntegrationHelper: add hooks to upload SC grids, cache uploads keyed by pointer/hash to avoid re-upload if unchanged.
- Tests: single-particle SC force vs CPU; stage-synchronous drift with SC matches CPU tolerance.

4c) Data/memory considerations
- Device species table: mass, charge, CCS_HSS/EHSS per gas; gas properties per domain.
- Align SoA buffers for coalesced access; reuse ion_state buffers for v/p.
- Per-stage SC uploads are heavy: prefer pinned host buffers and reuse device allocations.

5) Performance & validation
- Benchmark CPU vs GPU across ion counts; measure PCIe overhead with/without cached fields.
- Add GPU unit/acceptance tests: E-only drift, damping heating, simple SC force spot-check (if implemented).
- Update validation scripts to flag `enable_gpu` + unsupported forces as invalid.

6) Docs & release notes
- Mark GPU as experimental until 3–4 land; document hard gates and supported force matrix.
- Record fallback/refusal behavior in CONFIG_GUIDE and RELEASE_NOTES.

GPU Development Plan (v1.0 safety → v1.1 enablement)
====================================================

Goals
- v1.0: prevent invalid GPU runs, document limits.
- v1.1: enable GPU with feature parity (E/B + damping + space charge) and better caching.

1) Consolidate current state
- Hard-gate GPU when unsupported forces are present: space charge, damping, magnetic, composite/multiple electric forces; log “experimental” with reason. Runtime check in GPUIntegrationStrategy + CLI note.
- Surface GPU disablement in output metadata.

2) Minimal v1.0 safety
- Ensure reproducibility: per-ion RNG on GPU collisions or disable GPU collisions entirely.
- Disable GPU when space charge is enabled; warn clearly.

3) Architecture clean-up for v1.1
- Decouple GPU strategy from “exactly one ElectricFieldForce”: accept a vector of force functors; bail only on unsupported types.
- Persist device field uploads across domain switches; cache by (provider ptr, domain id).
- Unify integrator dispatch: no duplicate selection logic in GPUIntegrationStrategy; reuse SimulationEngine integrator choice.

4) Feature enablement sequence
- Multi-domain E/B: support multiple electric field forces with domain offsets; optional static B-field via provider.
- Damping on GPU: port friction kernel with Teff scaling; share CCS lookup path with CPU.
- Space charge on GPU: either GPU P³M or CPU SC + GPU E/B hybrid; must be stage-synchronous (RK4 first).
- Boundary actions on GPU: wire absorption/cylindrical helper into timestep.

5) Performance & validation
- Benchmark CPU vs GPU across ion counts; measure PCIe overhead with/without cached fields.
- Add GPU unit/acceptance tests: E-only drift, damping heating, simple SC force spot-check (if implemented).
- Update validation scripts to flag `enable_gpu` + unsupported forces as invalid.

6) Docs & release notes
- Mark GPU as experimental until 3–4 land; document hard gates and supported force matrix.
- Record fallback/refusal behavior in CONFIG_GUIDE and RELEASE_NOTES.

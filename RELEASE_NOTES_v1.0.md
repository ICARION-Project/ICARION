# ICARION v1.0 Release Notes

## Highlights
- Per-ion adaptive RK45 with deterministic RNG handling and per-ion timesteps (no global min-dt clamping); OpenMP determinism is covered by tests.
- Reproducibility improvements: inline resolved config snapshot in HDF5, runtime RK45 settings, RNG scope, physics/GPU thresholds, hashes for config/species/reaction DBs and field array files, derived summaries; output buffer cap option.
- GPU EHSS plumbing and parity/thermalization sanity checks retained for developers, but the primary runtime GPU path is **disabled for v1.0** (CPU-only in `SimulationEngine`; other GPU helpers remain buildable for dev/testing).
- New tests: RK45 per-ion dt/OpenMP determinism, SimulationEngine per-ion dt, batch dt fallback, RNG determinism across compaction (NoCollisions), GPU EHSS mapping/parity.
- Forces/Integrator: RK45 uses the SoA force path end-to-end; forces implement `compute_soa` and space-charge indexing is corrected (no AoS staging).

## Validation (CPU v1.0)
- Full physics validation is available in `validation/VALIDATION_REPORT_v1.0.md` (thermalization, IMS, Orbitrap, LQIT, quadrupole, FT-ICR, TOF, drift/gas flow, combined drift, reactions, space charge, CPU performance).
- GPU performance suite is skipped because the runtime GPU path is disabled in v1.0; configs remain for future releases.

## Known Limitations
- **RK45 performance/architecture:** RK45 now uses the SoA force path (no AoS staging) and per-ion state caches; adaptive remains per-ion on CPU (no batch path) by design. SimulationEngine remains monolithic. No performance regressions are measured in CI; see `validation/` for performance/physics sweeps.
- **External inputs embedding:** Config, species/reaction databases, and field arrays are now embedded verbatim in HDF5 alongside SHA256 hashes (larger files); per-ion seeds are not stored (only the seed scheme is recorded).
- **GPU status:** GPU collision (EHSS/HSS) and GPU integration are experimental; the primary runtime GPU path in `SimulationEngine` is disabled for v1.0. CUDA builds still compile helpers for dev/testing, with CPU fallback when batches are unavailable.
- **Magnetic field providers:** Analytical/uniform B-fields work; external B-field map providers are stubbed out and not wired in v1.0.
- **Batch constraints:** Batch integration/collision runs only when active ions share a uniform `dt`; mixed-dt runs fallback to per-ion paths.
- **Friction collision model:** Requires reduced mobility input; when absent, it falls back to CCS-based damping (approximate) and emits debug logs. Validated for low-to-moderate E/N (1–10 Td at 100–5000 Pa); use stochastic models outside this envelope or when diffusion matters.

## Reproducibility
- HDF5 metadata includes: resolved config JSON (runtime), RK45 runtime tolerances, RNG scheme/seed, physics handler names and GPU thresholds, domain/species/reaction counts, SHA256 hashes for config/species/reaction DBs and field array files.
- Output buffer cap configurable via `output.buffer_byte_cap` or `--buffer-byte-cap`; hard enforcement to avoid OOM.
- Per-ion RNG streams resize on compaction/addition; deterministic compaction tests run with NoCollisions.

## Notes
- GPU features require `-DUSE_GPU_ACCEL=ON` and a CUDA-capable device; GPU tests are only built when `USE_GPU_ACCEL` is enabled.
- Validation/benchmarks live under `validation/` and are not part of CI.

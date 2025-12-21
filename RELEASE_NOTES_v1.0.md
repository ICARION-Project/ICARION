# ICARION v1.0 Release Notes

## Highlights
- Per-ion adaptive RK45 (OpenMP-safe) with deterministic RNG handling and per-ion timesteps (no global min-dt clamping).
- Reproducibility improvements: inline resolved config snapshot in HDF5, runtime RK45 settings, RNG scope, physics/GPU thresholds, hashes for config/species/reaction DBs and field array files, derived summaries; output buffer cap option.
- GPU EHSS plumbing and parity/thermalization sanity checks retained for developers, but the runtime GPU path is **hard-disabled for v1.0** (any `enable_gpu=true` falls back to CPU).
- New tests: RK45 per-ion dt/OpenMP determinism, SimulationEngine per-ion dt, batch dt fallback, RNG determinism across compaction, GPU EHSS mapping/parity.
- Forces/Integrator: RK45 uses the SoA force path end-to-end; forces implement `compute_soa` and space-charge indexing is corrected (no AoS staging).

## Known Limitations
- **RK45 performance/architecture:** RK45 now uses the SoA force path (no AoS staging) and per-ion state caches; adaptive remains serial-only by design. SimulationEngine remains monolithic. No performance regressions are measured in CI; see `validation/` for performance/physics sweeps.
- **External inputs embedding:** Config, species/reaction databases, and field arrays are now embedded verbatim in HDF5 alongside SHA256 hashes (larger files); per-ion seeds are not stored (only the seed scheme is recorded).
- **GPU status:** GPU collision (EHSS/HSS) and GPU integration are experimental and runtime-disabled for v1.0; any GPU request is ignored and falls back to CPU.
- **Magnetic field providers:** Analytical/uniform B-fields work; external B-field map providers are stubbed out and not wired in v1.0.
- **Batch constraints:** Batch integration/collision runs only when active ions share a uniform `dt`; mixed-dt runs fallback to per-ion paths.

## Reproducibility
- HDF5 metadata includes: resolved config JSON (runtime), RK45 runtime tolerances, RNG scheme/seed, physics handler names and GPU thresholds, domain/species/reaction counts, SHA256 hashes for config/species/reaction DBs and field array files, buffer cap setting.
- Output buffer cap configurable via `output.buffer_byte_cap` or `--buffer-byte-cap`; hard enforcement to avoid OOM.
- Per-ion RNG streams resize on compaction/addition; deterministic tests cover compaction.

## Notes
- GPU features require `-DUSE_GPU_ACCEL=ON` and a CUDA-capable device; GPU tests skip otherwise.
- Validation/benchmarks live under `validation/` and are not part of CI.

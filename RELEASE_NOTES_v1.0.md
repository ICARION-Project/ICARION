# ICARION v1.0 Release Notes

## Highlights
- Per-ion adaptive RK45 (OpenMP-safe) with deterministic RNG handling and per-ion timesteps (no global min-dt clamping).
- Reproducibility improvements: inline resolved config snapshot in HDF5, runtime RK45 settings, RNG scope, physics/GPU thresholds, hashes for config/species/reaction DBs and field array files, derived summaries; output buffer cap option.
- GPU EHSS: geometry upload and species mapping enabled; experimental parity/thermalization sanity checks added.
- New tests: RK45 per-ion dt/OpenMP determinism, SimulationEngine per-ion dt, batch dt fallback, RNG determinism across compaction, GPU EHSS mapping/parity.

## Known Limitations
- **RK45 performance/architecture:** RK45 still uses AoS copies and carries state internally; SimulationEngine remains monolithic. No performance regressions are measured in CI; see `validation/` for performance/physics sweeps.
- **External inputs embedding:** Species/reaction databases and field arrays are hashed (SHA256) but not embedded; per-ion seeds are not stored (only the seed scheme is recorded).
- **GPU status:** GPU collision (EHSS/HSS) and GPU integration are experimental; GPU EHSS physics beyond parity checks is not validated.
- **Batch constraints:** Batch integration/collision runs only when active ions share a uniform `dt`; mixed-dt runs fallback to per-ion paths.

## Reproducibility
- HDF5 metadata includes: resolved config JSON (runtime), RK45 runtime tolerances, RNG scheme/seed, physics handler names and GPU thresholds, domain/species/reaction counts, SHA256 hashes for config/species/reaction DBs and field array files, buffer cap setting.
- Output buffer cap configurable via `output.buffer_byte_cap` or `--buffer-byte-cap`; hard enforcement to avoid OOM.
- Per-ion RNG streams resize on compaction/addition; deterministic tests cover compaction.

## Notes
- GPU features require `-DUSE_GPU_ACCEL=ON` and a CUDA-capable device; GPU tests skip otherwise.
- Validation/benchmarks live under `validation/` and are not part of CI.

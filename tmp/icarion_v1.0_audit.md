**Physical Correctness**
- ✅ Collisions/reactions now use the integrator’s accepted per-ion `dt`; pre-adaptation bias removed.
- ⚠ GPU EHSS uploads geometry with deterministic species→index map and warns on missing species; parity/thermalization sanity added, but physics still experimental.
- ✅ Missing mixture concentrations default to zero (no buffer-gas fallback).
- ✅ Single-gas HSS/OU unchanged and consistent.

**Numerical Soundness**
- ✅ RK45 parallel races avoided by forcing serial mode with warning.
- ✅ RNG pool resizes with ensemble to prevent desync after compaction/addition.
- ⚠ Adaptive timestep feedback still uses min suggested `dt` across ions (conservative, not tuned).
- ✅ Fixed-step paths remain deterministic with per-ion RNG.

**Architecture & Maintainability**
- ⚠ RK45 remains non-reentrant; mitigation is serial enforcement, not refactor.
- ⚠ SimulationEngine still large orchestrator; untouched.
- ✅ ForceRegistry separation retained.

**Performance & Parallelization**
- ⚠ RK45 scaling absent (serial-only); RK45 AoS copies/allocation overhead unchanged.
- ⚠ GPU collision helper thresholded; adaptive/mixed-`dt` runs may fall back to CPU.
- ✅ Fixed-step + SoA path unchanged.

**I/O & Reproducibility**
- ✅ HDF5 embeds resolved config snapshot (in-memory), integrator/RK45 params (runtime), RNG scope, physics handler/GPU thresholds, derived summaries, input hashes (config/species/reaction DB, field arrays).
- ✅ Output buffer cap (`output.buffer_byte_cap` / `--buffer-byte-cap`) with hard enforcement.
- ⚠ DB contents not embedded (only hashes); field arrays hashed, not stored; per-ion seeds not recorded (scheme only).

**Release Readiness (v1.0)**
- ✅ Prior blockers mitigated: RK45 races/dt misalignment, GPU EHSS misrepresentation (geometry upload), missing metadata.
- ⚠ Remaining: RK45 serial-only, RK45 performance, GPU EHSS still experimental, DB/field data not embedded.

Verdict: Blockers addressed; residual risks are performance/parallel limits and metadata completeness. Proceed with caution for v1.0; parallel RK45 and GPU EHSS not production-grade.

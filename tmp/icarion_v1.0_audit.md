**Physical Correctness**
- ⚠ Space-charge fields refresh at each RK4/RK45 stage on CPU; adaptive RK45 now runs stage-synchronously when `ICARION_ADAPTIVE_SC` is enabled (default). GPU paths deliberately drop SC because GPU is disabled for v1.0.
- ⚠ GPU integration is intentionally marked experimental/disabled for v1.0; only a single electric field force is supported and other forces (damping/B/SC) are ignored (`GPUIntegrationStrategy.cpp:119-133`). Acceptable only because GPU is not shipped for production.
- ✅ Boundary actions are invoked on exit/entry using geometry intersections (`SimulationEngine.cpp:507-612`); collisions/reactions still use pre-step environments by design.
- ⚠ Collisions/reactions use pre-integration domain/environment (`SimulationEngine.cpp:443-497`); mid-step domain crossings bias gas properties. Accepted as macro-step approximation and documented.
- ❌ Magnetic field providers are stubbed out (`MagneticFieldForce.cpp:73-79`); map-driven B-fields are ignored and replaced by analytical defaults.
- ⚠ Damping “Friction” mode falls back to ion CCS without species DB, emits debug to stderr (`DampingForce.cpp:146-255`); mixture CCS and logging are uncontrolled.
- ✅ Core collision kernels (HSS/EHSS/OU) match documented formulas and guard zero-velocity cases.

**Numerical Soundness**
- ✅ Trajectory output now stores per-ion times and species indices (flattened buffers); snapshots no longer mix adaptive-step states without timestamps.
- ⚠ GPU collision batching discards per-ion RNG streams and uses helper RNG (`GPUCollisionHandler.cpp:84-121`); tolerable only because GPU is disabled in v1.0.
- ⚠ Collisions/reactions assume constant `dt_used_per_ion` within the pre-step domain; mid-step domain changes or adaptive `dt` may bias rates because environments aren’t recomputed for substeps. Space-charge substep syncing covers RK4 and adaptive RK45 on CPU; GPU still drops SC.
- ✅ CPU RNG seeding is deterministic per ion (`SimulationEngine.cpp:430-433`); RK45 error control clamps `dt` sensibly (`RK45Strategy.cpp`).

**Architecture & Maintainability**
- ❌ `SimulationEngine` is a monolith entangling domains, forces, collisions, reactions, output, RNG, and OpenMP; boundary actions and field models are interleaved, increasing change risk.
- ⚠ GPU strategy is tightly coupled to registry shape (exactly one electric force) and duplicates integrator selection; acceptable only because GPU is flagged experimental/disabled for v1.0.
- ❌ Field models are instantiated both in `PhysicsSetup` and `DomainManager`; engine backfills them, violating SSOT and complicating lifetime guarantees.
- ⚠ Space-charge models are shared via registries and deduped by pointer comparison (`SimulationEngine.cpp:400-415`), which is brittle if models move.
- ✅ ForceRegistry/IForce layering is clear and SoA-aware; collision/reaction handlers follow Strategy with batch hooks.

**Performance & Parallelization**
- ❌ Domain lookup is O(domains) per ion twice per step (`SimulationEngine.cpp:443-525`) with no spatial index; dominates OpenMP runs at scale.
- ⚠ GPU integration caches field uploads when the provider pointer is stable but still refuses multi-domain or multi-force batches, so GPU often falls back and wastes PCIe bandwidth (tolerated because GPU is disabled for v1.0).
- ❌ Space-charge grid model recomputes and caches per-ion fields every step (`SpaceChargeGridModel.cpp:87-101`), even when unchanged, adding O(N) memory traffic before every substep.
- ✅ OpenMP loops use guided scheduling with chunking (~128) to reduce imbalance when many ions are inactive or born mid-run.
- ✅ SoA layouts and batch hooks exist, offering a path to vectorization once blockers are fixed.

**I/O & Reproducibility**
- ✅ Environment cache in `IonEnsemble` is refreshed from domain config each step and on domain switch and is documented as a macro-step approximation (CONFIG_GUIDE.md); mid-step domain changes for collisions/reactions remain approximate but accepted.
- ✅ HDF5 embeds config JSON plus species/reaction DBs and field arrays as blobs (plus hashes); external inputs are preserved for reruns (file size cost).
- ✅ Validation/analysis scripts now consume embedded inputs and new `species_id_indices` (schema-safe helpers added); legacy `species_ids` assumptions pruned.
- ✅ Species IDs now stored as indices in trajectory output; flattened buffers avoid per-step varlen strings (performance win).
- ✅ Metadata records git hash, build info, RNG scheme, config JSON embedding (`hdf5Writer.cpp:60-210`).

- **Release Readiness**
- ✅ v1.0 CPU-only release is acceptable: GPU is explicitly disabled/experimental; stage-synchronous SC works for RK4/RK45 on CPU (adaptive enabled by default via `ICARION_ADAPTIVE_SC`); boundary actions applied; per-ion times/species indices written; external inputs embedded; env cache approximation documented.
- ⚠ v1.1: spatial index for domains, quiet DampingForce logs, reintroduce deterministic GPU path with full force coverage (SC/damping/B/stochastic), and persist GPU field uploads across domain switches; consider compression/streaming for large embedded field grids.

Verdict: I would recommend publication at this stage because the CPU path is physically consistent and reproducible with documented macro-step environment assumptions, GPU acceleration is intentionally disabled and labeled experimental for v1.0, and all external inputs are embedded; remaining issues are performance (domain lookup, SC rebuilds) and future GPU coverage deferred to v1.1.

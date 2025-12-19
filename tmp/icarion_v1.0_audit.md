**Physical Correctness**
- ❌ Space-charge fields are frozen for the whole macro-step: `SimulationEngine.update_space_charge_models` recomputes once per step (`SimulationEngine.cpp:400-415`), but RK4/RK45 sub-stages reuse stale fields (`ForceRegistry.cpp:38-55`, `SpaceChargeGridModel.cpp:87-106`), so Coulomb forces lag ion motion.
- ❌ GPU integration ignores everything except a single electric field provider: `GPUIntegrationStrategy` requires exactly one `ElectricFieldForce` with a provider and applies only E/B; damping, space charge, magnetic forces, and composites are dropped (`GPUIntegrationStrategy.cpp:119-133`). Magnetic instruments and space-charge runs on GPU are invalid.
- ❌ Boundary physics are dead: `DomainManager` builds `BoundaryAction`s but the engine never invokes them; boundary exits simply deactivate ions (`SimulationEngine.cpp:507-547`), so configured reflections/absorptions are ignored.
- ❌ Collisions/reactions use pre-integration domain/environment (`SimulationEngine.cpp:443-497`); if an ion crosses a boundary mid-step, gas properties for stochastic events are wrong.
- ❌ Magnetic field providers are stubbed out (`MagneticFieldForce.cpp:73-79`); map-driven B-fields are ignored and replaced by analytical defaults.
- ⚠ Damping “Friction” mode falls back to ion CCS without species DB, emits debug to stderr (`DampingForce.cpp:146-255`); mixture CCS and logging are uncontrolled.
- ✅ Core collision kernels (HSS/EHSS/OU) match documented formulas and guard zero-velocity cases.

**Numerical Soundness**
- ❌ Output stores one global time stamp per snapshot while ions advance with per-ion/adaptive `dt`; per-ion times are not written (`SimulationEngine.cpp:321-329`, `hdf5Writer.cpp:91-195`), so trajectories mix states at different times and are not replayable.
- ❌ GPU collision batching discards per-ion RNG streams and uses helper RNG (`GPUCollisionHandler.cpp:84-121`), breaking CPU/GPU parity and reproducibility.
- ⚠ Collisions/reactions assume constant `dt_used_per_ion` within the pre-step domain; mid-step domain changes or adaptive `dt` may bias rates because environments aren’t recomputed for substeps.
- ✅ CPU RNG seeding is deterministic per ion (`SimulationEngine.cpp:430-433`); RK45 error control clamps `dt` sensibly (`RK45Strategy.cpp`).

**Architecture & Maintainability**
- ❌ `SimulationEngine` is a monolith entangling domains, forces, collisions, reactions, output, RNG, and OpenMP; boundary actions and field models are interleaved, increasing change risk.
- ❌ GPU strategy is coupled to registry shape (exactly one electric force) and duplicates integrator selection; adding any force disables GPU silently (`GPUIntegrationStrategy.cpp:119-133`).
- ❌ Field models are instantiated both in `PhysicsSetup` and `DomainManager`; engine backfills them, violating SSOT and complicating lifetime guarantees.
- ⚠ Space-charge models are shared via registries and deduped by pointer comparison (`SimulationEngine.cpp:400-415`), which is brittle if models move.
- ✅ ForceRegistry/IForce layering is clear and SoA-aware; collision/reaction handlers follow Strategy with batch hooks.

**Performance & Parallelization**
- ❌ Domain lookup is O(domains) per ion twice per step (`SimulationEngine.cpp:443-525`) with no spatial index; dominates OpenMP runs at scale.
- ❌ GPU integration reuploads field grids every batch (`GPUIntegrationHelper.cpp:70-150`) and refuses multi-domain or multi-force batches, so GPU often falls back and wastes PCIe bandwidth.
- ❌ Space-charge grid model recomputes and caches per-ion fields every step (`SpaceChargeGridModel.cpp:87-101`), even when unchanged, adding O(N) memory traffic before every substep.
- ⚠ OpenMP uses fixed `schedule(static,256)` regardless of active ion density, hurting balance with many inactive/born flags.
- ✅ SoA layouts and batch hooks exist, offering a path to vectorization once blockers are fixed.

**I/O & Reproducibility**
- ❌ SSOT is claimed but environment is cached in `IonEnsemble` and manually updated (`SimulationEngine.cpp:475-483`); missed updates desync from config.
- ❌ HDF5 stores hashes of external inputs but not the field arrays or DB contents; combined with missing per-ion times/seeds, long-term reproducibility is weak.
- ✅ Species IDs now stored as indices in trajectory output; flattened buffers avoid per-step varlen strings (performance win).
- ✅ Metadata records git hash, build info, RNG scheme, config JSON embedding (`hdf5Writer.cpp:60-210`).

**Release Readiness**
- ❌ Blockers for v1.0: fix GPU path to honor all active forces (at least space charge, damping/magnetic) or disable GPU; store per-ion times (or enforce uniform time) in outputs; refresh space-charge fields within substeps or document the approximation explicitly.
- ⚠ v1.1: spatial index for domains, cache/persist GPU field uploads, quiet DampingForce logs, embed external inputs, improve OpenMP scheduling.
- ✅ Boundary actions are now applied using geometry intersections; CPU collision/reaction kernels and electric-only force paths are usable for single-domain, fixed-`dt` studies.

Verdict: I would not recommend publication at this stage, because GPU support drops key forces, outputs misrepresent adaptive timelines, and space-charge handling is physically incomplete. Boundary handling is improved but does not remove the main blockers.

# Space Charge: Current State and Next Steps

## Current State
- **CPU**
  - Direct Coulomb force (O(N²)) is exact but slow; valid for all geometries.
  - Grid-based Poisson solver (`SpaceChargeSolver` + `SpaceChargeGrid`) uses a rectangular box grid with simple boundary assumptions (no cylindrical/Orbitrap BC). Accuracy depends on embedding the domain in a box; not geometry-accurate for cylinders/Orbitrap.
- **GPU**
  - P³M helper (`GPUSpaceChargeP3M`) implemented and unit-tested (parity/benchmarks), but **not wired into SimulationEngine**.
  - Limitations: rectangular domain only, P2G/G2P effectively CIC only, ignores `active` flags, no bounds checks, buffers not resized if N grows.
  - GPU boundary helper and GPU space-charge helper both not called from main loop.

## Suitability
- For cylindrical/Orbitrap domains: the box-grid CPU solver is a poor match; direct O(N²) is correct but slow. GPU P³M is not integrated and also box-only.

## Next Steps (Plan)
1) **Documentation hygiene:** audit src/ comments for performance claims/limitations to avoid misleading users.
2) **Wire GPU P³M:** call `try_gpu_space_charge` with active-ion count, add bounds/active handling, and resize buffers on demand.
3) **Geometry-aware solvers:** add cylindrical boundary handling on CPU (grid in r–z or cylindrical Poisson) and GPU (cylindrical discretization or Orbitrap-specific potential).
4) **Config/control:** add explicit switches/thresholds for space charge backend (direct vs grid vs GPU), and log warnings for geometry mismatch.
5) **Testing/validation:** extend parity/physics tests for cylindrical/Orbitrap cases; add performance benchmarks.

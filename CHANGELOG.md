# ICARION Changelog

All notable changes to this project follow [Semantic Versioning](https://semver.org/) and are documented in this file.

## [1.0.1] - 2026-06-02

### Fixed

- **EHSS mixture-path collision rate** (`EHSSCollisionHandler`): multi-component gas
  configs incorrectly used the bulk relative velocity `|v_ion − v_drift|` (≈ 0 for cold
  ions) instead of a thermally-sampled per-component velocity, severely underestimating
  the EHSS collision rate.  Each neutral velocity is now sampled from the
  Maxwell-Boltzmann distribution and reused consistently for rate selection and momentum
  transfer.

- **EHSS first-hit atom selection** (`CollisionKernels::detect_atom_collision`): when
  multiple overlapping atoms intersect the collision trajectory the first atom in array
  order was returned instead of the geometrically earliest contact (smallest `s_hit`).
  This biased the contact normal and momentum-transfer cosine distribution.  The fix
  iterates all atoms and selects the one with the smallest `s_hit`.

## [1.0.0] - 2025-12-04

- Initial v1.0.0 release tag.
- Stable JSON configuration schema (SSOT) with validation.
- Reproducible HDF5 output (version, git hash, build flags, RNG seed).
- Instruments: IMS, RF quadrupole, Orbitrap, TOF, LQIT, FT-ICR.
- Physics: HSS/EHSS/Langevin/HardSphere collisions; Arrhenius reactions; space-charge models (direct/grid/GPU experimental).
- Integrators: RK4, RK45, Boris (GPU backend experimental).
- Validation suite shipped in-repo with scripts/results/figures.
- Documentation: architecture, CLI usage, config schema, validation report.

[1.0.1]: https://github.com/ICARION-Project/ICARION/releases/tag/v1.0.1
[1.0.0]: https://github.com/ICARION-Project/ICARION/releases/tag/v1.0.0

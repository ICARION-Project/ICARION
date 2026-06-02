# ICARION Changelog

All notable changes to this project follow [Semantic Versioning](https://semver.org/) and are documented in this file.

## [1.0.1] - 2026-06-02

### Fixed

- **EHSS mixture-path collision rate** (`EHSSCollisionHandler`): the rate calculation for
  multi-component (or single-component list) gas configurations incorrectly used the bulk
  relative velocity `|v_ion − v_drift|` instead of the thermally-sampled per-component
  velocity `|v_ion − v_neutral_i|`.  For cold ions (`T_ion ≪ T_gas`), the bulk relative
  velocity is near zero, which caused the EHSS collision rate — and therefore the
  thermalization rate — to be severely underestimated.  This manifested as EHSS appearing
  *slower* than HSS in standard thermalization benchmarks, contrary to the expected
  behaviour for molecules with `σ_EHSS > σ_HSS`.  The fix samples each neutral velocity
  from the Maxwell-Boltzmann distribution at the gas temperature, consistent with the HSS
  handler which was not affected by this bug.  The sampled velocity is reused for the
  actual scattering step so that rate selection and momentum transfer are fully consistent.

## [1.0.0] - 2025-12-04

- Initial v1.0.0 release tag.
- Stable JSON configuration schema (SSOT) with validation.
- Reproducible HDF5 output (version, git hash, build flags, RNG seed).
- Instruments: IMS, RF quadrupole, Orbitrap, TOF, LQIT, FT-ICR.
- Physics: HSS/EHSS/Langevin/HardSphere collisions; Arrhenius reactions; space-charge models (direct/grid/GPU experimental).
- Integrators: RK4, RK45, Boris (GPU backend experimental).
- Validation suite shipped in-repo with scripts/results/figures.
- Documentation: architecture, CLI usage, config schema, validation report.

[1.0.0]: https://github.com/ICARION-Project/ICARION/releases/tag/v1.0.0

# ICARION Changelog

All notable changes to this project follow [Semantic Versioning](https://semver.org/) and are documented in this file.

## [1.0.2] - 2026-01-07

- Performance: batch RK4 CPU path and thread-safe profiling updates.
- Benchmarking: thread/ion scaling runner and plotting helper.
- Cleanup: suppress unused parameter warnings in boundary/waveform helpers.

## [1.0.0] - 2025-12-04

- Initial v1.0.0 release tag.
- Stable JSON configuration schema (SSOT) with validation.
- Reproducible HDF5 output (version, git hash, build flags, RNG seed).
- Instruments: IMS, RF quadrupole, Orbitrap, TOF, LQIT, FT-ICR.
- Physics: HSS/EHSS/Langevin/HardSphere collisions; Arrhenius reactions; space-charge models (direct/grid/GPU experimental).
- Integrators: RK4, RK45, Boris (GPU backend experimental).
- Validation suite shipped in-repo with scripts/results/figures.
- Documentation: architecture, CLI usage, config schema, validation report.

[1.0.2]: https://github.com/ICARION-Project/ICARION/releases/tag/v1.0.2
[1.0.0]: https://github.com/ICARION-Project/ICARION/releases/tag/v1.0.0

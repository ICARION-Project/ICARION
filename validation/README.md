# ICARION Validation Suite

**Version:** 1.0.0  
**Last Updated:** 2026-01-09  
**Branch:** `release/v1.0.0-prep`  
**Goal:** Systematic validation of physics and performance for the v1.0.0 release. The suite ships with this repo and does not require external assets.

**Validation vs CTests**

This directory contains high-fidelity validations (runtime up to ~30 minutes). For regression tests, see `tests/` (CTests).

| Aspect | CTests (`tests/`) | Validation Suite (`validation/`) |
|--------|-------------------|----------------------------------|
| Purpose | CI/CD regression testing | Scientific validation for publications |
| Runtime | Seconds to minutes (benchmarks/GPU can be longer) | Up to ~30 minutes |
| Ensemble | Varies (10-50k ions depending on test) | 1k-10k ions typical |
| Output | Pass/Fail | Plots, tables, metrics |
| Example | `test_gas_flow_transport.cpp` | `validate_gas_flow_transport.py` |

---

## Scope and Status

| Suite | Scope | Configs | Entry Script | Status | Notes |
|-------|-------|---------|--------------|--------|-------|
| Thermalization | HSS/EHSS thermalization | 90 (generated) | `scripts/thermalization/run_thermalization_tests.sh` | Complete | Configs generated, not checked in |
| Instruments | IMS, LQIT, Orbitrap, Quadrupole, TOF | `configs/instruments/` | `scripts/run_instrument_suite.sh` | Complete | FT-ICR configs are generated, not checked in |
| Transport | Drift, gas flow, mixtures | `configs/physics/transport/` (generated) | `scripts/run_physics_suite.sh` | Partial | Diffusion validation pending |
| Space Charge | Coulomb expansion, solver parity | `configs/physics/spacecharge/` | `scripts/run_physics_suite.sh spacecharge` | Complete | Samples checked in |
| Reactions | First-order + bimolecular | `configs/physics/reactions/` | `scripts/run_physics_suite.sh reactions` | Complete | Samples checked in |
| Performance (CPU) | Scaling + model overhead | 18 | `scripts/performance/run_performance_suite.sh` | Complete | CPU only in v1.0.0 |
| Performance (GPU) | Speedup vs CPU | 31 (generated) | `scripts/performance/run_performance_suite.sh --gpu-only` | Skipped | Runtime GPU disabled in v1.0.0 |

---

## Quick Start

Primary runners (cover the full suite areas):
```bash
cd /home/chsch95/ICARION/validation

# Instruments
./scripts/run_instrument_suite.sh

# Physics (includes thermalization, transport, reactions, space charge)
./scripts/run_physics_suite.sh

# Performance (CPU only in v1.0.0)
./scripts/performance/run_performance_suite.sh
```

Example single-runner (covered by the physics runner, but useful for focused runs):
```bash
# Generate thermalization configs (if missing)
python3 scripts/thermalization/generate_thermalization_configs.py

# Run only thermalization
./scripts/thermalization/run_thermalization_tests.sh full
```

---

## Suites

### Primary runners

#### Instruments
- Configs: `configs/instruments/`
- Runner: `scripts/run_instrument_suite.sh`
- Analysis: `scripts/run_instrument_analysis.sh`
- Instrument-specific generators and analyzers:
  - IMS: `scripts/instrumentation/configs/generate_ims_configs.py`, `scripts/instrumentation/ims/analyze_ims_drift.py`
  - LQIT: `scripts/instrumentation/configs/generate_lqit_configs.py`, `scripts/instrumentation/configs/generate_lqit_mass_scan_suite.py`
  - Orbitrap: `scripts/instrumentation/configs/generate_orbitrap_configs.py`, `scripts/instrumentation/orbitrap/analyze_orbitrap_frequencies.py`
  - Quadrupole: `scripts/instrumentation/configs/generate_quadrupole_stability_map.py`, `scripts/instrumentation/quadrupole/analyze_quadrupole_stability_map.py`
  - TOF: `scripts/instrumentation/tof/analyze_tof_flight_time.py`
  - FT-ICR: `scripts/instrumentation/configs/generate_fticr_configs.py`, `scripts/instrumentation/fticr/analyze_fticr_frequencies.py`
- Results (baseline): `results/v1.0.0_test/instruments/`
- Results (per-run): `runs/<run-id>/results/instruments/`

#### Physics (includes thermalization, transport, reactions, space charge)
- Runner: `scripts/run_physics_suite.sh`
- Analysis: `scripts/run_physics_analysis.sh`
- Configs: `configs/physics/`
- Results: `results/` (suite-specific folders)

#### Performance
- CPU: `scripts/performance/run_performance_suite.sh`, `scripts/performance/run_performance_analysis.sh`
- GPU: `scripts/performance/generate_gpu_performance_configs.py` (generated configs; GPU runtime disabled in v1.0.0)
- Results (baseline): `results/v1.0.0_test/performance/logs/`
- Results (per-run): `runs/<run-id>/results/performance/logs/`

### Physics components (covered by the physics runner)

#### Thermalization
- Config generator: `scripts/thermalization/generate_thermalization_configs.py`
- Runner: `scripts/thermalization/run_thermalization_tests.sh`
- Analysis: `scripts/thermalization/analyze_thermalization_complete.py`
- Configs: `configs/physics/thermalization/` (generated)
- Results: `results/v1.0.0_test/physics/thermalization/`

#### Transport and Mixtures
- Config generator: `scripts/instrumentation/configs/generate_transport_drift_configs.py`
- Drift analysis: `scripts/physics/analyze_transport_drift.py`
- Gas flow: `scripts/physics/validate_gas_flow_transport.py`
- Combined drift: `scripts/physics/validate_combined_drift.py`
- Mixture mobility: `scripts/physics/validate_gas_mixture_mobility.py`
- Mixture thermalization: `scripts/physics/validate_mixture_thermalization.py`
- Results: `results/combined_drift/`, `results/gas_flow_transport/`, `results/gas_mixture_mobility/`, `results/mixture_thermalization/`
- Logs: `logs/`

#### Space Charge
- Configs: `configs/physics/spacecharge/`
- Scripts: `scripts/physics/validate_space_charge_adaptive_parity.py`, `scripts/physics/bench_space_charge_adaptive.py`, `scripts/physics/analyze_spacecharge.py`
- Results: `results/v1.0.0_test/physics/spacecharge/`

#### Reactions
- Configs: `configs/physics/reactions/`
- Scripts: `scripts/physics/validate_reaction_kinetics.py`, `scripts/physics/analyze_reactions.py`, `scripts/physics/plot_reaction_validation.py`
- Results: `results/physics/reactions/`

---

## Directory Structure

```
validation/
|-- configs/
|   |-- physics/
|   |   |-- thermalization/     # generated by scripts/thermalization/generate_thermalization_configs.py
|   |   |-- transport/          # generated by scripts/instrumentation/configs/generate_transport_drift_configs.py
|   |   |-- spacecharge/        # committed samples
|   |   `-- reactions/          # committed samples
|   |-- instruments/
|   |   |-- ims/
|   |   |-- lqit/
|   |   |-- orbitrap/
|   |   |-- quadrupole/
|   |   |-- tof/
|   |   `-- quadrupole_old_partial/  # legacy configs
|   `-- performance/            # generated benchmarks
|-- figures/                    # plots (generated)
|-- logs/                       # analysis summaries
|-- results/                    # test outputs (HDF5, logs)
|-- runs/                       # per-run structured artifacts (generated, gitignored)
`-- README.md                   # this file
```

**Config policy:** Generators emit JSON configs under `validation/configs/<category>/`. If a folder is empty, re-run the relevant generator. Legacy pre-generated configs may still live under `validation/scripts/configs/` for reference.

**Results layout:**
- `validation/results/v1.0.0_test/instruments/<instrument>/...`
- `validation/results/<suite>/` for gas flow, combined drift, mixture mobility/thermalization, reactions
- `validation/results/v1.0.0_test/performance/logs/`
- `validation/results/v1.0.0_test/` for frozen baselines
- Legacy runs may exist under `validation/results/physics/<suite>/`

---

## Per-run output folders (recommended)

Historically, many scripts wrote directly into `validation/logs/`, `validation/figures/`, and `validation/results/`, which makes it hard to compare runs over time.

The suite runners can create a per-run folder and export `ICARION_VALIDATION_RUN_DIR` so validators write into a clean, timestamped structure:

```
validation/runs/<run-id>/
|-- logs/
|-- figures/physics/
`-- results/
```

You can override the identifier and output location:

```bash
./scripts/run_physics_suite.sh --run-id my_tag
./scripts/run_physics_suite.sh --run-dir /tmp/icarion_physics_run

./scripts/run_instrument_suite.sh --run-id my_tag
./scripts/run_instrument_suite.sh --run-dir /tmp/icarion_instrument_run

./scripts/performance/run_performance_suite.sh --run-id my_tag
./scripts/performance/run_performance_suite.sh --run-dir /tmp/icarion_performance_run
```

To force legacy/baseline output locations, pass `--baseline-output` to the instrument or performance runners.

---

## Thermalization Test Modes

- Quick mode: representative subset
- Subset mode: one ion species across temperatures
- Full mode: all 90 configs

---

## Notes

- GPU runtime is disabled in v1.0.0; any `enable_gpu=true` falls back to CPU.
- FT-ICR configs are provided under `validation/configs/instruments/fticr/` (you can regenerate them via `scripts/instrumentation/configs/generate_fticr_configs.py`).
- Diffusion validation is not implemented yet.

## Publishing results (recommended)

If you want to publish validation results without checking in large `*.h5` files, see [validation/PUBLISHING.md](PUBLISHING.md).
The flow is: run into `validation/runs/<run-id>/`, then export a git-friendly bundle via `scripts/export_run_bundle.sh`.

# Validation scripts

This folder contains the orchestration scripts for the ICARION validation suite.

## Where to start (entrypoints)

Use these top-level runners in normal workflows:

- `run_all_suites.sh` for the combined physics, instrument, and performance run
- `run_physics_suite.sh` (+ `run_physics_analysis.sh`)
- `run_instrument_suite.sh` (+ `run_instrument_analysis.sh`)
- `performance/run_performance_suite.sh` (+ `performance/run_performance_analysis.sh`)

They are designed to be runnable from **any** current working directory.

## Output locations

Many configs use `output.folder: "validation/results/..."`.
To avoid accidentally creating nested paths like `validation/validation/results/...`, the runners invoke `icarion_main` with the **repo root** as working directory.

Recommended: per-run output folders (do not overwrite baselines):

```
validation/runs/<run-id>/
|-- logs/
|-- figures/
`-- results/
```

Legacy/baseline output is still supported via `--baseline-output` (where implemented).

## Folder map

- `physics/`: physics validators + analyzers (Python)
- `instrumentation/`: instrument-specific generators/runners/analyzers
- `thermalization/`: thermalization generators + runners
- `performance/`: performance benchmark runners

## Notes on legacy helpers

There are a few standalone helper scripts (e.g. single-config runners, ad-hoc benchmarks). These write under `validation/results/...` and should not create `validation/scripts/validation/...` or `validation/validation/...` output trees.

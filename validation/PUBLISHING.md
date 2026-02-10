# Publishing validation results (without large `.h5` files)

Goal: publish scientifically useful validation evidence (report + plots + logs) while keeping the repository lightweight and results reproducible.

## Recommended policy

Keep in git:
- Validation report: `validation/VALIDATION_REPORT_v1.0.md`
- Curated plots (PNG/SVG) referenced by the report
- Curated summary logs/tables (TXT/CSV/JSON)
- Reproducibility metadata and configs **used for the run** (small JSON)

Do **not** keep in git:
- Large simulation artifacts (especially `*.h5`)

## Why nested folders happen

Many configs use a relative output like:

```json
{"output": {"folder": "validation/results/..."}}
```

If `icarion_main` is launched with a working directory inside `validation/`, that relative path becomes `validation/validation/results/...`.
The suite runners fix this by invoking `icarion_main` from the **repo root**.

## Reproducible runs

Run suites with a stable run id:

```bash
./validation/scripts/run_physics_suite.sh --run-id v1.0_physics
./validation/scripts/run_instrument_suite.sh --run-id v1.0_instruments
./validation/scripts/performance/run_performance_suite.sh --run-id v1.0_performance
```

Each run writes a `manifest.json` into the run directory with:
- git commit + dirty flag
- command line
- run id + run directory

## Export a git-friendly bundle

To publish a run without `.h5`, export it:

```bash
./validation/scripts/export_run_bundle.sh \
  --run-dir validation/runs/v1.0_physics \
  --out-dir validation/published/v1.0_physics
```

This copies:
- `manifest.json`
- `logs/`
- `figures/`
- JSON config snapshots found under `results/`

but excludes `*.h5`.

## Suggested published layout (in git)

Create a curated, versioned bundle under `validation/published/` and reference *only* those artifacts from the report.

Recommended naming:

```
validation/published/
`-- v1.0/
    |-- physics/
    |-- instruments/
    `-- performance/
```

Each bundle should contain (small files only):
- `manifest.json` (git commit, command line, run id)
- `README.md` (human summary + how to reproduce)
- `logs/` (key text outputs)
- `figures/` (PNGs/SVGs referenced by the report)
- `configs/` (snapshotted JSON configs used for the run, when available)

Example export:

```bash
./validation/scripts/export_run_bundle.sh \
  --run-dir validation/runs/v1.0_physics \
  --out-dir validation/published/v1.0/physics
```

Then update `validation/VALIDATION_REPORT_v1.0.md` to link to:
- `validation/published/v1.0/physics/figures/...`
- `validation/published/v1.0/physics/logs/...`

This keeps the report stable and avoids linking to transient `validation/runs/...` paths.

## Tips

- Prefer generators + parameter grids over checking in thousands of configs.
- For publication-grade runs, snapshot the exact configs used (the suite runners and some validators already do this for sensitive cases like space charge / reactions).

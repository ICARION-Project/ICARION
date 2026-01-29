# ICARION Tools

## ccs_precompute

Derive gas-specific CCS maps from a reference CCS value and store them in a
species database JSON.

### Build

Use the CMake target:
```bash
cmake --build build --target ccs_precompute
```

The binary is written to your build directory (for example `build/ccs_precompute`
with single-config generators).

### Usage

```bash
build/ccs_precompute \
  --input species.json \
  --output species_with_ccs.json \
  --species H3O+,NH4+ \
  --ref-gas He \
  --ref-ccs-A2 110.0 \
  [--model HSS|EHSS] \
  [--override] \
  [--n-orientations 300]
```

### Options

- `--input`: Species JSON containing a top-level `species` object.
- `--output`: Output JSON path.
- `--species`: Comma-separated species IDs to update.
- `--ref-gas`: Reference gas for `--ref-ccs-A2`.
- `--ref-ccs-A2`: Reference CCS in A^2.
- `--model`: `HSS` (default) or `EHSS`.
- `--override`: Overwrite existing `CCS_<model>` entries.
- `--n-orientations`: Number of orientations for the EHSS OAPA projection
  (default: 300).

### Supported gases

`He`, `Ar`, `CO2`, `Ne`, `N2`, `O2`, `H2O`

### Output fields

For each selected ion species (charge != 0), the tool writes:

- `CCS_reference_gas`: The reference gas name.
- `CCS_model`: `HSS` or `EHSS`.
- `CCS_HSS` or `CCS_EHSS`: Gas-to-CCS map in A^2.

When `--model EHSS` is selected, the tool attempts a geometry-based OAPA
projection using `geometry_file` (resolved relative to the input JSON). If
geometry is missing or fails to load, a warning is emitted and the HSS-derived
map is kept under `CCS_EHSS`.

### References

See `docs/CONFIG_GUIDE.md` for how the CCS maps are referenced in configurations.

---

## ehss_samples_precompute

Precompute orientation-sampled projection areas for EHSS and store them in a
per-species JSON file that can be referenced from `EHSS_samples_file`.

### Build

```bash
cmake --build build --target ehss_samples_precompute
```

### Usage

```bash
build/ehss_samples_precompute \
  --input species.json \
  --output samples.json \
  --species H3O+ \
  [--n-orientations 300] \
  [--n-samples 8000] \
  [--seed 12345]
```

### Options

- `--input`: Species JSON containing a top-level `species` object.
- `--output`: Output samples JSON path.
- `--species`: Target species ID (must include `geometry_file`).
- `--n-orientations`: Number of random orientations to sample.
- `--n-samples`: Monte Carlo samples per orientation.
- `--seed`: RNG seed (default: 12345).

### Notes

- Uses Monte Carlo estimation of the projected union area per orientation.
- Samples are computed for each supported gas radius (same list as above).
- Output fields include `orientations_quat` and `areas_by_gas_m2`.

---

## ehss_samples_sanity

Sanity checks for EHSS samples: summary stats, optional OAPA comparison, and
histogram plots (Agg backend).

### Usage

```bash
tools/ehss_samples_sanity.py --input samples.json --gas N2 --bins 30
tools/ehss_samples_sanity.py --input samples.json --geometry-file molecule.json --plot-dir /tmp/ehss_hists
```

### Options

- `--input`: Samples JSON path.
- `--gas`: Gas key to inspect (default: all).
- `--bins`: Histogram bins (default: 20).
- `--geometry-file`: Geometry JSON for OAPA mean comparison.
- `--plot`: Output PNG for a single gas histogram (requires `--gas`).
- `--plot-dir`: Output directory for per-gas histograms.

### References

See `docs/CONFIG_GUIDE.md` for how `EHSS_samples_file` is used at runtime.

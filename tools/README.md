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

---

## ehss_offline_precompute

Precompute single-gas EHSS collision samples for `EHSS_offline_samples_file`.
The runtime uses these files instead of atom-sphere ray casting in single-gas
EHSS runs.

### Build

```bash
cmake --build build --target ehss_offline_precompute
```

### Usage

```bash
build/src/ehss_offline_precompute \
  --input species.json \
  --output h3o_ehss_he.h5 \
  --species H3O+ \
  --gas He \
  [--format hdf5|json] \
  [--n-orientations 1000] \
  [--n-sigma-samples 10000] \
  [--n-mu-samples 2000]
```

### Runtime Format

JSON and HDF5 use the same logical fields:

- `format`: `ehss_offline_samples`
- `version`: `1`
- `species_id`: species key used to generate the file
- `gas`: gas species for which the file is valid
- `units`: `sigma_eff_m2,mu`
- `orientations_quat`: shape `(N, 4)`, sampled molecular orientations
- `sigma_eff_m2`: shape `(N)`, orientation-specific EHSS cross section
- `mu_samples`: shape `(N, M)`, conditional samples of `cos(theta)` for accepted first-contact collisions

The precompute step finds the first atom contact for each accepted projected ray
and stores only the resulting scattering distribution. The runtime therefore
does not need molecular geometry on the offline path.

---

## interaction_potential_precompute

Precompute interaction-potential collision lookup data for `ipm_samples_file`.
The tool integrates classical trajectories for each orientation and relative
velocity bin, then writes a runtime HDF5 lookup file.

### Build

```bash
cmake --build build --target interaction_potential_precompute
```

### Usage

```bash
build/src/interaction_potential_precompute \
  --input species.json \
  --species H3O+ \
  --output h3o_ipm_he.h5 \
  --gas He \
  [--potential lj1264] \
  [--n-orientations 50] \
  [--n-trials 20000] \
  [--v-bins 32] \
  [--store-full-cdf]
```

### Runtime Format

HDF5 attributes:

- `format`: `ipm_offline_samples`
- `version`: `1`
- `species_id`, `gas`, `gas_model`
- `mixing_rule`, `polarization`, `potential`, `param_model`
- `units`: `logv, sigma_mt_m2, b_max_m, dp_SI`
- `n_orientations`, `v_bins`, `seed`

Required datasets:

- `logv_bins`: shape `(K)`, natural-log relative-speed bins
- `orientations_quat`: shape `(N, 4)`
- `sigma_mt_m2`: shape `(N, K)`, momentum-transfer cross section
- `b_max_m`: shape `(N, K)`, impact-parameter integration cutoff

Optional full-CDF datasets:

- `cdf_offsets`: shape `(N, K)`
- `cdf_counts`: shape `(N, K)`
- `cdf_values`: shape `(S)`, cumulative weights for sampled momentum kicks
- `dp_samples`: shape `(S, 3)`, sampled momentum kicks in SI units

Fallback statistics:

- `dp_stats`: shape `(N, K, 4)`, weighted mean parallel/perpendicular momentum kick and variances

If the full CDF is absent, the runtime samples from `dp_stats`; if present, it
uses `cdf_*` and `dp_samples` for the higher-fidelity momentum-kick distribution.

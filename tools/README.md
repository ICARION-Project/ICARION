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

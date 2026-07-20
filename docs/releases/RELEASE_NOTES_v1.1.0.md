# ICARION v1.1.0 Release Notes

This release covers the audited v1.1 runtime, collision, reaction, TIMS,
analysis, documentation, schema, validation, and packaging surface. FieldSolver
work from the development branch is intentionally deferred to v1.2.

## Highlights

- Reactions now support equilibrium-linked dynamic reverse channels. A forward
  reaction marked with `equilibrium=true` and SI thermochemistry metadata can
  generate its reverse channel automatically, with the reverse effective rate
  computed from the linked forward rate and the temperature-dependent
  equilibrium constant.
- `InteractionPotentialModel` can be used to calculate ion-neutral collisions under the
  influence of long-range interactions. Species use `ipm_samples_file`, and runtime HDF5
  sample tables identify themselves with `format=ipm_offline_samples`.
- EHSS and IPM offline precompute tools are documented as part of the public
  v1.1 surface. The tools generate runtime sample tables instead of requiring
  geometry/contact work during the simulation step.
- Stochastic collision subcycling supports multiple collision opportunities per
  simulation timestep. This is useful for high collision-rate regimes, but for
  highest accuracy the physical timestep should still keep unresolved multiple
  events rare.
- Deep collision diagnostics can write per-ion collision summaries and optional
  sampled/full event rows under `/analysis/deep_collision` for momentum-transfer
  audits.
- TIMS is available as a distinct instrument type with an analytical axial ramp
  field model, axial gas-flow parsing, and a basic TIMS elution example.
- Analysis utilities now include maintained scripts for run summaries, spectra,
  trap stability, transport diagnostics, collision diagnostics, reaction
  kinetics, and space charge cloud broadening.

## Validation

- Dynamic equilibrium reactions have a dedicated validation scenario that
  checks the thermochemistry-derived reverse rate against the expected
  equilibrium population split.
- IPM validation generates a fresh mini sample table, verifies the
  `ipm_offline_samples` HDF5 format, checks positive momentum-transfer cross
  sections, runs an IMS-style simulation, and confirms finite trajectory and
  momentum diagnostics.
- TIMS elution validation checks mobility-sorted release order, peak separation,
  and approximate agreement with the analytical ramp-fraction estimate.

## Data And Examples

- `examples/ims/ims_ipm_basic.json` runs against the global species database.
- `examples/ims/ims_tims_basic.json` demonstrates the TIMS axial ramp field with
  axial gas flow and RF confinement.
- Small bundled IPM example tables are available under
  `data/molecules/precomputed_ipm/` for H3O+ and PentanalH+ in He.
- Larger IPM tables should be regenerated per release or study using
  `interaction_potential_precompute`; the local helper
  `tmp/run_large_ipm_precompute.sh` is intentionally not part of the release
  tree.
- JSON Schema coverage was aligned with v1.1 runtime controls, including
  `space_charge_model`, InteractionPotentialModel controls, compact/minimal
  output, and deep collision diagnostics.

## Packaging

- CMake/CPack package version is `1.1.0`.
- Release packaging keeps the v1.0-style assets: Linux `.deb`, portable Linux
  `.tar.gz`, Windows `.zip`, SHA256 files, and the minimal launchers.
- Installed command-line tools include `icarion`, `icarion_main`,
  `ccs_precompute`, `ehss_samples_precompute`, `ehss_offline_precompute`, and
  `interaction_potential_precompute`.
- Release packages exclude generated example `results/`, legacy
  `precomputed_lr1264/` tables, and analysis result folders.

## Compatibility Notes

- `physics.enable_reactions=false` no longer loads the global reaction database
  fallback, so disabled-reaction examples no longer emit unrelated reaction
  validation warnings.

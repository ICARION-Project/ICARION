# ICARION v1.1.0 Release Notes

Status: draft. This file currently covers the audited reactions, collision, and
TIMS updates only; later runtime, output, fieldsolver, validation, and
documentation blocks should extend it before the final v1.1.0 release.

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

## Compatibility Notes

- `physics.enable_reactions=false` no longer loads the global reaction database
  fallback, so disabled-reaction examples no longer emit unrelated reaction
  validation warnings.
- Old `LR1264_*` sample names and old `lj1264_offline_samples` assumptions are
  not part of the audited v1.1 public surface. Use `InteractionPotentialModel`
  and `ipm_samples_file`.

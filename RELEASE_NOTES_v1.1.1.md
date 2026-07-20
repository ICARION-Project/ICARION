# ICARION v1.1.1 Release Notes

This patch release fixes trajectory output scheduling and HDF5 batching. It
does not change the public configuration API, HDF5 schema, integrators, or
collision models.

## Fixed

- Removed overlapping step-based and time-based trajectory triggers that could
  write adjacent duplicate snapshots.
- Separated trajectory sampling from HDF5 buffer flushing. Reaching buffer
  capacity no longer creates an additional sample.
- Buffered scheduled snapshots are appended in batches instead of forcing
  pathological one-snapshot writes at every output interval.

## Scientific Impact

The issue affects trajectory sampling density, file size, and analyses that
use the stored time axis, including frequency-domain reconstruction. It does
not alter the simulated trajectory integration or collision physics between
output points.

Runs whose HDF5 trajectories were produced with v1.1.0 should be checked for
adjacent time samples before reuse in signal or frequency analysis. Re-running
affected simulations with v1.1.1 is recommended.

## Compatibility

- Configuration compatibility: unchanged from v1.1.0.
- HDF5 schema: unchanged from v1.1.0.
- Package and runtime version: 1.1.1.

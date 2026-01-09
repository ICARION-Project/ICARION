# Analysis Examples

Small, runnable examples for ICARION HDF5 trajectories. Focused on quick visualizations with switches for subsampling and time striding so large ion counts remain manageable.

## Requirements
- Python 3
- `h5py`, `numpy`, `matplotlib`
- Optional: `pillow` (GIF writer) or `ffmpeg` (MP4 writer) for animations

## Scripts
- `analysis/plot_trajectories.py`: Static plot (XY and XZ projection plus optional 3D) for a selected subset of ions.
- `analysis/animate_trajectories.py`: Simple projection animation (GIF/MP4) with striding and sampling to keep size down.
- `analysis/elimination_histograms.py`: Elimination-time histograms and channel counts (radial/axial, Orbitrap inner/outer); also produces per-species breakdowns.
- `analysis/mean_positions.py`: Mean radial and axial positions over time (per species) with ion downsampling and time striding.
- `analysis/thermalization_check.py`: Temperature vs time plus speed distribution vs Maxwell-Boltzmann over the last timesteps.
- `analysis/ims_mobility.py`: Drift mobility / reduced mobility estimation from IMS trajectories (arrival times, K, K0).

## Example invocations
```bash
# Static plot from an LQIT simulation
python analysis/plot_trajectories.py \
  --traj validation/results/v1.0_test/instruments/lqit/lqit_vacuum_q0.700_a0.000.h5 \
  --out analysis/output/lqit_trajectory.png \
  --species H3O+ \
  --max-ions 80 \
  --time-stride 5

# XY animation as GIF (Pillow writer) with automatic downsampling
python analysis/animate_trajectories.py \
  --traj validation/results/v1.0_test/instruments/lqit/lqit_vacuum_q0.700_a0.000.h5 \
  --out analysis/output/lqit_xy.gif \
  --projection xy \
  --max-ions 60 \
  --time-stride 4 \
  --frame-step 2

# Elimination histogram (radial vs axial); auto-detects Orbitrap, otherwise radial/axial
python analysis/elimination_histograms.py \
  --traj validation/results/v1.0_test/instruments/orbitrap/orbitrap_ReserpineH+_V3500.00.h5 \
  --out analysis/output/orbitrap_elimination.png \
  --species ReserpineH+ \
  --max-ions 300 \
  --log-bins \
  --tol-frac 0.05
# Per-species plot is also created as *_per_species.png (disable with --no-per-species).

# Mean radial/axial positions (per species)
python analysis/mean_positions.py \
  --traj validation/results/v1.0_test/instruments/lqit/lqit_vacuum_q0.700_a0.000.h5 \
  --out analysis/output/lqit_mean_positions.png \
  --species H3O+ \
  --time-stride 5 \
  --max-ions 300

# Thermalization check: temperature curve + MB overlay on last speeds
python analysis/thermalization_check.py \
  --traj validation/results/v1.0_test/physics/thermalization/therm_hss_H3Op_300K_20.0Pa.h5 \
  --out analysis/output/thermalization_check.png \
  --time-stride 10 \
  --window 10

# IMS mobility / reduced mobility from arrival times
python analysis/ims_mobility.py \
  --traj validation/results/v1.0_test/instruments/ims/ims_hss_48Vm_200Pa.h5 \
  --out analysis/output/ims_mobility.png \
  --out-csv analysis/output/ims_mobility.csv \
  --species H3O+ \
  --bins 60 \
  --time-stride 2
```

See `--help` on the scripts for more options (e.g., ion caps per species, RNG seed, projection).

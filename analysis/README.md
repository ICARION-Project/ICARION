# Analysis Examples

Small, runnable examples for ICARION HDF5 trajectories. Focused on quick visualizations with switches for subsampling and time striding so large ion counts remain manageable.

## Requirements
- Python 3
- `h5py`, `numpy`, `matplotlib`

## Scripts
- `analysis/plot_trajectories.py`: Static plot (XY and XZ projection plus optional 3D) for a selected subset of ions.
- `analysis/animate_trajectories.py`: Simple projection animation (GIF/MP4) with striding and sampling to keep size down.
- `analysis/elimination_histograms.py`: Elimination-time histograms and channel counts (radial/axial, Orbitrap inner/outer); also produces per-species breakdowns.
- `analysis/mean_positions.py`: Mean radial and axial positions over time (per species) with ion downsampling and time striding.

## Example invocations
```bash
# Static plot from an LQIT simulation
python analysis/plot_trajectories.py \
  --traj results/lqit/lqit_trajectories.h5 \
  --out analysis/output/lqit_trajectory.png \
  --species Mass80_CCS25 Mass88_CCS50 \
  --max-ions 80 \
  --time-stride 5

# XY animation as GIF (Pillow writer) with automatic downsampling
python analysis/animate_trajectories.py \
  --traj results/lqit/lqit_trajectories.h5 \
  --out analysis/output/lqit_xy.gif \
  --projection xy \
  --max-ions 60 \
  --time-stride 4 \
  --frame-step 2

# Elimination histogram (radial vs axial); auto-detects Orbitrap, otherwise radial/axial
python analysis/elimination_histograms.py \
  --traj results/orbitrap/orbitrap_trajectories.h5 \
  --out analysis/output/orbitrap_elimination.png \
  --species ReserpineH \
  --max-ions 300 \
  --log-bins \
  --tol-frac 0.05
# Per-species plot is also created as *_per_species.png (disable with --no-per-species).

# Mean radial/axial positions (per species)
python analysis/mean_positions.py \
  --traj results/lqit/lqit_trajectories.h5 \
  --out analysis/output/lqit_mean_positions.png \
  --species Mass80_CCS25 Mass88_CCS50 \
  --time-stride 5 \
  --max-ions 300
```

See `--help` on the scripts for more options (e.g., ion caps per species, RNG seed, projection).

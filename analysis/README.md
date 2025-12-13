# Analysis Examples

Kleine, lauffähige Auswertungen für die ICARION HDF5-Trajektorien. Der Fokus liegt auf schnellen Visualisierungen mit ein paar Schaltern, damit auch Runs mit vielen Ionen handhabbar bleiben (Subsampling, Zeit-Striding).

## Voraussetzungen
- Python 3
- `h5py`, `numpy`, `matplotlib`

## Skripte
- `analysis/plot_trajectories.py`: Statischer Plot (XY- und XZ-Projektion plus optional 3D) für eine Auswahl an Ionen.
- `analysis/animate_trajectories.py`: Einfache Animation (GIF/MP4) einer Projektion. Nutzt Striding und Sampling, damit große Dateien nicht explodieren.
- `analysis/elimination_histograms.py`: Histogramme der Eliminationszeiten und Kanal-Zählungen (radial/axial, Orbitrap inner/outer).
- `analysis/mean_positions.py`: Mittlere radiale und axiale Positionen über die Zeit (pro Species), mit Ion-Downsampling und Time-Striding.

## Beispielaufrufe
```bash
# Statischer Plot aus einer LQIT-Simulation
python analysis/plot_trajectories.py \
  --traj results/lqit/lqit_trajectories.h5 \
  --out analysis/output/lqit_trajectory.png \
  --species Mass80_CCS25 Mass88_CCS50 \
  --max-ions 80 \
  --time-stride 5

# XY-Animation als GIF (Pillow-Writer) mit automatischem Downsampling
python analysis/animate_trajectories.py \
  --traj results/lqit/lqit_trajectories.h5 \
  --out analysis/output/lqit_xy.gif \
  --projection xy \
  --max-ions 60 \
  --time-stride 4 \
  --frame-step 2

# Eliminationshistogramm (radial vs axial); erkennt Orbitrap automatisch, sonst radial/axial
python analysis/elimination_histograms.py \
  --traj results/orbitrap/orbitrap_trajectories.h5 \
  --out analysis/output/orbitrap_elimination.png \
  --species ReserpineH \
  --max-ions 300 \
  --log-bins \
  --tol-frac 0.05

# Mittlere radiale/axiale Positionen (pro Species)
python analysis/mean_positions.py \
  --traj results/lqit/lqit_trajectories.h5 \
  --out analysis/output/lqit_mean_positions.png \
  --species Mass80_CCS25 Mass88_CCS50 \
  --time-stride 5 \
  --max-ions 300
```

Siehe `--help` bei den Skripten für weitere Optionen (z. B. max. Ionen pro Species, Zufallssaat, Projektion).

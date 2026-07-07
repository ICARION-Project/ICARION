# ICARION Analysis

Maintained Python analysis utilities for ICARION trajectory, spectra, transport, collision, reaction, and space charge outputs.

## Requirements

- Python 3
- `numpy`, `h5py`, `matplotlib`
- `scipy` for Gaussian fitting helpers
- Optional: `pillow` or `ffmpeg` for animations

Install the analysis dependencies with:

```bash
python3 -m pip install -r requirements-analysis.txt
```

Install the package and console entry points in editable mode with:

```bash
python3 -m pip install -e .
```

## Shared Modules

- `analysis.trajectory`: trajectory HDF5 opening, species selection, position loading, domain/environment metadata, arrival detection, death-time helpers.
- `analysis.physics`: constants, temperature/mobility conversions, histogram edges, Gaussian fitting.
- `analysis.spectra`: histogram spectra, FFT frequency spectra, peak picking.
- `analysis.stability`: trap/filter transmission, survival, and radial-envelope helpers.

`analysis.common` remains as a compatibility facade for older imports.

## Release Scope

Tracked and maintained scripts:

- `plot_trajectories.py`
- `animate_trajectories.py`
- `elimination_histograms.py`
- `mean_positions.py`
- `thermalization_check.py`
- `ims_mobility.py`
- `arrival_time_distribution.py`
- `run_report.py`
- `spectrum_analysis.py`
- `trap_stability.py`
- `transport_diagnostics.py`
- `collision_diagnostics.py`
- `reaction_kinetics.py`
- `spacecharge_diagnostics.py`

Local or experimental scripts should live under `tmp/analysis/`, which is ignored by git.

## Examples

Static trajectory plot:

```bash
python3 analysis/plot_trajectories.py \
  --traj validation/results/v1.1.0_test/instruments/lqit/lqit_vacuum_q0.700_a0.000.h5 \
  --out analysis/output/lqit_trajectory.png \
  --species H3O+ \
  --max-ions 80 \
  --time-stride 5
```

Installed equivalent:

```bash
icarion-plot-trajectories \
  --traj validation/results/v1.1.0_test/instruments/lqit/lqit_vacuum_q0.700_a0.000.h5 \
  --out analysis/output/lqit_trajectory.png \
  --species H3O+ \
  --max-ions 80 \
  --time-stride 5
```

Trajectory animation:

```bash
python3 analysis/animate_trajectories.py \
  --traj validation/results/v1.1.0_test/instruments/lqit/lqit_vacuum_q0.700_a0.000.h5 \
  --out analysis/output/lqit_xy.gif \
  --projection xy \
  --max-ions 60 \
  --time-stride 4 \
  --frame-step 2
```

Arrival-time distribution:

```bash
python3 analysis/arrival_time_distribution.py \
  --traj validation/results/v1.1.0_test/instruments/ims/ims_hss_48Vm_200Pa.h5 \
  --out analysis/output/arrival_time_dist.png \
  --out-csv analysis/output/arrival_time_dist.csv \
  --species H3O+ \
  --bins 50
```

IMS mobility:

```bash
python3 analysis/ims_mobility.py \
  --traj validation/results/v1.1.0_test/instruments/ims/ims_hss_48Vm_200Pa.h5 \
  --out analysis/output/ims_mobility.png \
  --out-csv analysis/output/ims_mobility.csv \
  --species H3O+ \
  --bins 60
```

`ims_mobility.py` derives field strength from HDF5 metadata by default:
`EN_Td` plus gas pressure/temperature is preferred, then axial/drift voltage
over domain length. `--field-VM` and `--field-V` remain explicit manual
overrides and are marked as such in the CSV `field_source` column.

Run summary:

```bash
python3 analysis/run_report.py \
  --traj validation/results/v1.1.0_test/instruments/ims/ims_hss_48Vm_200Pa.h5 \
  --out analysis/output/run_report.csv
```

TOF spectrum:

```bash
python3 analysis/spectrum_analysis.py \
  --traj validation/results/v1.1.0_test/instruments/tof/tof_H3O+_V2000.h5 \
  --mode tof \
  --acceleration-voltage 2000 \
  --flight-distance 0.1 \
  --out analysis/output/tof_spectrum.png \
  --out-csv analysis/output/tof_spectrum.csv
```

Orbitrap/FTICR-style frequency spectrum:

```bash
python3 analysis/spectrum_analysis.py \
  --traj validation/results/v1.1.0_test/instruments/orbitrap/orbitrap_basic.h5 \
  --mode frequency \
  --coordinate z \
  --resample uniform \
  --out analysis/output/frequency_spectrum.png \
  --out-csv analysis/output/frequency_spectrum.csv
```

Frequency mode rejects nonuniform time grids unless `--resample uniform` is
specified. The generated sidecar `<out-csv>.meta.json` records input time-step
statistics, resampling `dt`, frequency resolution, and Nyquist frequency.

Trap/filter stability:

```bash
python3 analysis/trap_stability.py \
  --traj validation/results/v1.1.0_test/instruments/quadrupole/quad_a+0.0000_q0.7000.h5 \
  --radius 0.01 \
  --out analysis/output/trap_stability.png \
  --out-csv analysis/output/trap_stability.csv
```

Transport diagnostics:

```bash
python3 analysis/transport_diagnostics.py \
  --traj validation/results/v1.1.0_test/instruments/ims/ims_hss_48Vm_200Pa.h5 \
  --out analysis/output/transport_diagnostics.png \
  --out-csv analysis/output/transport_diagnostics.csv
```

Collision diagnostics:

```bash
python3 analysis/collision_diagnostics.py \
  --traj validation/results/v1.1.0_test/physics/deep_collision.h5 \
  --out-csv analysis/output/collision_diagnostics.csv
```

Reaction kinetics:

```bash
python3 analysis/reaction_kinetics.py \
  --traj validation/results/v1.1.0_test/physics/reactions/reaction_demo.h5 \
  --out analysis/output/reaction_kinetics.png \
  --out-csv analysis/output/reaction_kinetics.csv
```

Space charge cloud broadening:

```bash
python3 analysis/spacecharge_diagnostics.py \
  --traj validation/results/v1.1.0_test/physics/spacecharge/spacecharge_demo.h5 \
  --out analysis/output/spacecharge_diagnostics.png \
  --out-csv analysis/output/spacecharge_diagnostics.csv
```

Real data sweep:

```bash
python3 analysis/run_realdata_sweep.py \
  --results-dir results \
  --out-dir analysis/output/realdata_sweep_fresh
```

The sweep runner writes a fresh `summary.json`, marks expected negative path
checks explicitly, and avoids carrying stale failures forward.

## Tests

Run the Python analysis checks with stdlib `unittest`:

```bash
python3 -m unittest discover -s tests/analysis
```

Release-oriented analysis gates:

```bash
python3 -m unittest discover -s tests/analysis -p test_analysis_release.py
python3 -m unittest discover -s tests/analysis -p test_analysis_scientific.py
ctest --test-dir build -R 'PythonAnalysis' --output-on-failure
```

Cheap syntax/import gate:

```bash
python3 -m py_compile $(git ls-files 'analysis/*.py')
python3 - <<'PY'
import importlib, subprocess
mods = [p[:-3].replace('/', '.') for p in subprocess.check_output(['git', 'ls-files', 'analysis/*.py'], text=True).splitlines()]
for mod in mods:
    importlib.import_module(mod)
print(f"imported {len(mods)} analysis modules")
PY
```

# ICARION Launcher Guide

ICARION release packages include minimal launchers for users who do not want to
start every simulation and quick analysis from a terminal. The launchers are
small convenience tools. They do not replace the command-line interface or the
Python analysis scripts.

Use the launcher when you want to:

- select a JSON configuration file
- start one ICARION run
- inspect the run log while it runs
- create basic plots from an existing trajectory HDF5 file

Use the terminal or Python scripts directly when you need:

- batch runs
- custom command-line flags
- publication-quality figures
- instrument-specific analysis parameters
- reproducible scripted post-processing

## Windows Zip

Download and extract `icarion-<version>-windows-x64.zip`.

If the file came from the GitHub Actions artifact page, extract the outer
artifact zip first. The inner `icarion-<version>-windows-x64.zip` is the actual
package.

Start the launcher by double-clicking:

```text
ICARION-Launcher.cmd
```

The Windows package layout should look like this:

```text
bin/
analysis/
data/
examples/
schema/
docs/
requirements-analysis.txt
ICARION-Launcher.cmd
ICARION-Launcher.ps1
```

## Linux Tarball

Download and extract `icarion-<version>-Linux.tar.gz`, then run:

```bash
chmod +x ICARION-Launcher.sh
./ICARION-Launcher.sh
```

The Linux tarball layout should look like this:

```text
bin/
share/icarion/analysis/
share/icarion/data/
share/icarion/examples/
share/icarion/schema/
share/icarion/requirements-analysis.txt
docs/
ICARION-Launcher.sh
```

The Linux launcher uses `zenity` for file pickers when it is available. Without
`zenity`, it falls back to terminal prompts.

## Debian Package

Install the package:

```bash
sudo apt install ./icarion_<version>_amd64.deb
```

Start the launcher:

```bash
icarion-launcher
```

The installed examples and analysis scripts are under `/usr/share/icarion`.

## Running A Simulation

Choose `Run simulation`, then select a JSON config file. A good first test is:

Windows:

```text
examples\ims\ims_basic.json
```

Linux tarball:

```text
share/icarion/examples/ims/ims_basic.json
```

Debian package:

```text
/usr/share/icarion/examples/ims/ims_basic.json
```

The launcher runs ICARION from the package root so relative paths inside the
example configs continue to work. Logs are written to:

```text
launcher-logs/
```

The exact simulation output location depends on the selected config. Example
configs usually write under a `results/` directory.

## Running Basic Analysis

Analysis actions use Python. Install the Python requirements before using the
analysis buttons.

Windows:

```powershell
py -3 -m pip install -r .\requirements-analysis.txt
```

Linux tarball:

```bash
python3 -m pip install -r ./share/icarion/requirements-analysis.txt
```

Debian package:

```bash
python3 -m pip install -r /usr/share/icarion/requirements-analysis.txt
```

Then choose one of the analysis actions and select an ICARION trajectory file
with extension `.h5` or `.hdf5`.

Available actions:

- `IMS mobility...`: arrival-time histogram plus inferred mobility, with a
  per-species CSV summary.
- `Arrival times...`: arrival-time distributions, including a per-species plot.
- `Trajectories...`: quick trajectory plot with conservative ion limits.
- `Mean positions...`: mean position over time plus CSV output.
- `Eliminations...`: elimination histograms, including per-species output when
  available.
- `Animate...`: creates a GIF animation with conservative frame and ion limits.

Analysis outputs are written to:

```text
analysis-output/
```

Analysis logs are written to:

```text
launcher-logs/
```

On Windows, image outputs are shown in a preview window when the analysis
finishes successfully. The preview includes buttons to close the preview, open
the image, or open the output folder. GIF files may need to be opened with the
system viewer.

## Terminal Equivalents

Everything the launcher does can also be done from a terminal.

Run a Windows simulation:

```powershell
.\bin\icarion.exe .\examples\ims\ims_basic.json
```

Run a Linux tarball simulation:

```bash
./bin/icarion ./share/icarion/examples/ims/ims_basic.json
```

Run an installed Debian simulation:

```bash
icarion /usr/share/icarion/examples/ims/ims_basic.json
```

Run Windows IMS analysis directly:

```powershell
py -3 .\analysis\ims_mobility.py --traj .\results\ims\ims_trajectories.h5 --out .\analysis-output\ims.png --out-csv .\analysis-output\ims.csv
```

Run Linux tarball arrival-time analysis directly:

```bash
python3 ./share/icarion/analysis/arrival_time_distribution.py \
  --traj ./results/ims/ims_trajectories.h5 \
  --out ./analysis-output/arrival-times.png \
  --out-csv ./analysis-output/arrival-times.csv
```

Show all options for a script:

```bash
python3 /usr/share/icarion/analysis/ims_mobility.py --help
```

## Limitations

The launchers intentionally expose only safe defaults. They are useful for a
first run, a quick sanity check, or a basic figure. They are not designed as a
complete graphical analysis environment.

Important limitations:

- Python is required for analysis and is not bundled into the ICARION binary
  package.
- Analysis scripts expect ICARION trajectory HDF5 files with the relevant
  trajectory, species, and domain metadata.
- IMS mobility inference needs meaningful drift-field and drift-length metadata,
  or explicit parameters when running the Python script directly.
- TIMS or other zero-field setups can produce invalid or infinite mobility
  estimates if treated as a simple drift-tube IMS analysis.
- Animation can be slow for large files. The launcher limits ions and frames to
  keep the output usable.
- The GUI exposes only common analysis actions. The Python scripts expose more
  options through `--help`.

For careful analysis, keep the launcher output as a first look and then rerun
the relevant Python script from a terminal with explicit options.

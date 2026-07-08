# Release packages and launcher

For new users, the release packages are the suggested starting point. They avoid
local CMake setup and include examples, schemas, documentation, launchers, and
basic analysis scripts.

## Package types

ICARION releases can provide:

- a Linux Debian package: `icarion_<version>_amd64.deb`
- a portable Linux archive: `icarion-<version>-Linux.tar.gz`
- a Windows archive: `icarion-<version>-windows-x64.zip`
- matching `.sha256` checksum files

Use the release asset for your platform. GitHub Actions artifacts may contain an
extra outer zip; extract that first and then extract the inner distributable
archive.

## Debian package

Install:

```bash
sudo apt install ./icarion_<version>_amd64.deb
```

Installed commands include:

- `icarion`
- `icarion-launcher`
- `icarion_main`
- `ccs_precompute`
- `ehss_samples_precompute`
- `ehss_offline_precompute`
- `interaction_potential_precompute`

Installed resources are under `/usr/share/icarion`:

- `/usr/share/icarion/examples`
- `/usr/share/icarion/schema`
- `/usr/share/icarion/data`
- `/usr/share/icarion/analysis`

Run an installed example:

```bash
icarion /usr/share/icarion/examples/ims/ims_basic.json
```

Start the launcher:

```bash
icarion-launcher
```

The installed launcher writes simulation output, analysis output, and logs under
`~/ICARION-runs` by default. Set `ICARION_RUN_DIR` to change this location.

## Portable Linux archive

Extract and run:

```bash
tar -xf icarion-<version>-Linux.tar.gz
cd icarion-<version>-Linux
chmod +x ICARION-Launcher.sh
./ICARION-Launcher.sh
```

Terminal run:

```bash
./bin/icarion ./share/icarion/examples/ims/ims_basic.json
```

The launcher uses `zenity` for file pickers when available and falls back to
terminal prompts otherwise.

## Windows archive

Extract `icarion-<version>-windows-x64.zip` and double-click:

```text
ICARION-Launcher.cmd
```

Terminal run:

```powershell
.\bin\icarion.exe .\examples\ims\ims_basic.json
```

The Windows package contains `bin/`, `examples/`, `schema/`, `data/`,
`analysis/`, and launcher scripts.

Release packages intentionally exclude generated example `results/`, legacy
development sample folders, and analysis result folders. Small curated IPM
example tables under `data/molecules/precomputed_ipm/` are included.

## Launcher scope

The launcher is intentionally minimal. It is useful for:

- selecting a JSON configuration file
- starting one simulation
- watching the run log
- running basic analysis on an existing HDF5 file

Use the CLI or Python scripts directly for batch runs, parameter sweeps,
publication figures, or fully reproducible post-processing.

## Analysis requirements

The analysis actions require Python 3 and the package requirements shipped with
the release. Installing these packages in a virtual environment is recommended,
especially on Linux distributions that protect the system Python installation.

Create and activate a virtual environment if you do not already use one:

```bash
python3 -m venv .venv
. .venv/bin/activate
python -m pip install --upgrade pip
```

Windows:

```powershell
py -3 -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -r .\requirements-analysis.txt
```

Portable Linux:

```bash
python -m pip install -r ./share/icarion/requirements-analysis.txt
```

Debian package:

```bash
python -m pip install -r /usr/share/icarion/requirements-analysis.txt
```

The same scripts can be run directly for more control:

```bash
python3 /usr/share/icarion/analysis/ims_mobility.py --help
```

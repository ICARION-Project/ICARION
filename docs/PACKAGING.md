# Packaging ICARION

ICARION can be distributed as an installable CMake/CPack package. The default
packaging setup builds:

- a Debian package on Linux (`.deb`)
- a portable install archive (`.tar.gz`)

## Build Packages

From the repository root:

```bash
cmake -B build-package -S . \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr \
  -DUSE_GPU_ACCEL=OFF \
  -DICARION_BUILD_TESTS=OFF
cmake --build build-package -j"$(nproc)"
cmake --build build-package --target package
```

The generated packages are written to `build-package/`.

## GitHub Releases

The `Release Packages` GitHub Actions workflow builds release assets
automatically for tags matching `v*`.

```bash
git tag v1.0.1-beta.1
git push origin v1.0.1-beta.1
```

The workflow creates a draft GitHub Release with:

- `icarion_<version>_amd64.deb`
- `icarion-<version>-Linux.tar.gz`
- `icarion-<version>-windows-x64.zip`
- matching `.sha256` files

The same workflow can also be started manually from the GitHub Actions tab. A
manual run uploads the packages as workflow artifacts, but does not create a
GitHub Release.

## Windows Zip

The Windows release artifact is built on GitHub Actions with MSVC and vcpkg.
It is distributed as a `.zip` first, because missing runtime DLLs or package
layout issues are easier to diagnose before adding an installer.

The zip contains a flat distribution layout:

- `bin/`
- `analysis/`
- `data/`
- `examples/`
- `schema/`
- `docs/`
- `requirements-analysis.txt`
- top-level `README.md`, `CHANGELOG.md`, and `LICENSE`
- `ICARION-Launcher.cmd` and `ICARION-Launcher.ps1`

For the simplest Windows use, double-click:

```text
ICARION-Launcher.cmd
```

The launcher opens a file picker for a JSON config, starts `bin\icarion.exe`,
and streams the run log in a small window. It can also run a small set of
common analysis scripts for an existing trajectory file and writes the output to
`analysis-output/`.

After extracting the zip, run from a terminal:

```powershell
.\bin\icarion.exe .\examples\ims\ims_basic.json
```

When downloading from the GitHub Actions **Artifacts** section, GitHub wraps the
release asset in an additional artifact zip. Extract that outer zip first; the
inner `icarion-*-windows-x64.zip` is the distributable package. GitHub Release
assets are uploaded directly without this extra wrapper.

## Portable Linux Archive

The Linux `.tar.gz` archive uses the same portable layout as the Windows zip:

- `bin/`
- `share/icarion/analysis`
- `share/icarion/data`
- `share/icarion/examples`
- `share/icarion/schema`
- `docs/`
- `share/icarion/requirements-analysis.txt`
- `ICARION-Launcher.sh`

After extracting the archive:

```bash
chmod +x ICARION-Launcher.sh
./ICARION-Launcher.sh
```

The launcher starts a simulation or runs one of the packaged analysis scripts.
It uses `zenity` when available and otherwise falls back to terminal prompts.
Direct terminal use is also supported:

```bash
./bin/icarion ./share/icarion/examples/ims/ims_basic.json
```

## Install A Debian Package

```bash
sudo apt install ./build-package/icarion_*.deb
```

Installed command-line tools:

- `icarion` convenience alias for `icarion_main`
- `icarion-launcher`
- `icarion_main`
- `ccs_precompute`
- `ehss_samples_precompute`

Installed resource layout:

- `/usr/share/icarion/analysis`
- `/usr/share/icarion/data`
- `/usr/share/icarion/examples`
- `/usr/share/icarion/schema`
- `/usr/share/doc/ICARION`

Example run after installation:

```bash
icarion /usr/share/icarion/examples/ims/ims_basic.json
```

Minimal launcher after installation:

```bash
icarion-launcher
```

## Minimal Launcher Analysis

The packaged launchers are intentionally small convenience tools. They are meant
for opening a config, starting a run, selecting a trajectory file, and producing
basic plots without typing the full command line.

Available analysis actions:

- IMS mobility
- arrival-time distributions, including per-species output
- trajectory plots
- mean-position plots
- elimination histograms
- trajectory animation

Analysis results are written next to the package or current run context:

- `analysis-output/` for figures, GIFs, and CSV files
- `launcher-logs/` for launcher and analysis logs

The analysis actions require Python 3 and the packages listed in
`requirements-analysis.txt`. They are not bundled into the ICARION binary
package.

On Windows:

```powershell
py -3 -m pip install -r .\requirements-analysis.txt
```

For the portable Linux archive:

```bash
python3 -m pip install -r ./share/icarion/requirements-analysis.txt
```

For the Debian package:

```bash
python3 -m pip install -r /usr/share/icarion/requirements-analysis.txt
```

The same analysis scripts can be run directly from a terminal. This exposes all
script options and is the recommended path for more detailed or reproducible
post-processing:

```powershell
py -3 .\analysis\ims_mobility.py --traj .\run.h5 --out .\analysis-output\ims.png --out-csv .\analysis-output\ims.csv --help
```

```bash
python3 ./share/icarion/analysis/arrival_time_distribution.py --traj ./run.h5 --out ./analysis-output/arrival.png --help
```

For an installed Debian package:

```bash
python3 /usr/share/icarion/analysis/ims_mobility.py --help
```

### Analysis Limitations

The launchers are not a full analysis workbench. They use conservative defaults
so that large trajectory files remain usable on normal laptops, especially for
animation. For publication-quality figures, batch processing, custom filters, or
instrument-specific parameters, run the Python scripts directly and inspect
their `--help` output.

Some analysis scripts require trajectory HDF5 files with the expected ICARION
trajectory, species, and domain metadata. IMS mobility inference additionally
depends on meaningful drift-field and drift-length information. TIMS or other
zero-field configurations can produce invalid or infinite mobility estimates
unless the relevant physical parameters are supplied explicitly through the
Python script options.

The launcher preview only shows generated image files that the desktop toolkit
can load. GIF animations are written to disk and can be opened with the system
viewer when preview support is limited.

## Notes

The Debian package uses `dpkg-shlibdeps` through CPack to derive runtime shared
library dependencies from the linked binaries. Build dependencies such as
CMake, compilers, Eigen, nlohmann-json, and C++ headers are still needed only on
the machine that creates the package.

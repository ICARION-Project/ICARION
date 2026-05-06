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
- `data/`
- `examples/`
- `schema/`
- `docs/`
- top-level `README.md`, `CHANGELOG.md`, and `LICENSE`
- `ICARION-Launcher.cmd` and `ICARION-Launcher.ps1`

For the simplest Windows use, double-click:

```text
ICARION-Launcher.cmd
```

The launcher opens a file picker for a JSON config, starts `bin\icarion.exe`,
and streams the run log in a small window.

After extracting the zip, run from a terminal:

```powershell
.\bin\icarion.exe .\examples\ims\ims_basic.json
```

When downloading from the GitHub Actions **Artifacts** section, GitHub wraps the
release asset in an additional artifact zip. Extract that outer zip first; the
inner `icarion-*-windows-x64.zip` is the distributable package. GitHub Release
assets are uploaded directly without this extra wrapper.

## Install A Debian Package

```bash
sudo apt install ./build-package/icarion_*.deb
```

Installed command-line tools:

- `icarion` convenience alias for `icarion_main`
- `icarion_main`
- `ccs_precompute`
- `ehss_samples_precompute`

Installed resource layout:

- `/usr/share/icarion/data`
- `/usr/share/icarion/examples`
- `/usr/share/icarion/schema`
- `/usr/share/doc/ICARION`

Example run after installation:

```bash
icarion /usr/share/icarion/examples/ims/ims_basic.json
```

## Notes

The Debian package uses `dpkg-shlibdeps` through CPack to derive runtime shared
library dependencies from the linked binaries. Build dependencies such as
CMake, compilers, Eigen, nlohmann-json, and C++ headers are still needed only on
the machine that creates the package.

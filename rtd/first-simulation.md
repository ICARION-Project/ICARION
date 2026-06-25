# First simulation

This page walks through the basic workflow of running an ICARION simulation.

## 1. Choose a configuration file

ICARION simulations are defined by JSON configuration files. Example configurations are available in the repository, for example under `examples/ims/`.

A typical first run is an ion mobility simulation:

```bash
./build/src/icarion_main examples/ims/ims_basic.json
```

Depending on your installation (e.g., when you use the Linux Debian package), the executable may also be available as:

```bash
icarion examples/ims/ims_basic.json
```

The input configuration files are described in more detail in the
[configuration guide](configuration.md). For installed packages and launchers,
see [Release packages and launcher](release-packages.md).

## 2. Validate the configuration

Before running longer simulations, validate the configuration without propagating ions:

```bash
./build/src/icarion_main --dry-run examples/ims/ims_basic.json
```

or:

```bash
./build/src/icarion_main --validate-config examples/ims/ims_basic.json
```

If ICARION is installed as a package, replace `./build/src/icarion_main` with `icarion`.
For the full list of CLI flags, see [CLI reference](cli-reference.md).

## 3. Run with a fixed seed

ICARION uses an explicit random number generator seed for reproducibility. To override the configured `simulation.rng_seed` from the command line:

```bash
./build/src/icarion_main --seed 42 examples/ims/ims_basic.json
```

## 4. Inspect the output

A simulation typically writes an HDF5 trajectory file and a resolved configuration snapshot. For `examples/ims/ims_basic.json`, the default output is `results/ims/ims_trajectories.h5` plus `results/ims/ims_trajectories.config.json`.

Simulation files are written as HDF5 files; see [Output files](output-files.md)
and [Validation](validation.md) for details. If the run fails or produces
unexpected output, use [Troubleshooting](troubleshooting.md).

# CLI reference

ICARION simulations are started from a JSON configuration file. In an installed
package the command is `icarion`; in a local build tree the executable is usually
`./build/src/icarion_main`.

```bash
icarion examples/ims/ims_basic.json
./build/src/icarion_main examples/ims/ims_basic.json
```

Use the build tree path when you compiled ICARION from source and did not install
it. Use `icarion` after installing a release package or after adding the build
output to `PATH`.

## Basic checks

Validate a configuration before spending time on a run:

```bash
icarion --dry-run examples/ims/ims_basic.json
icarion --validate-config examples/ims/ims_basic.json
icarion --validate-schema examples/ims/ims_basic.json
```

The checks cover different layers:

- `--dry-run` loads the configuration and exits before the simulation loop.
- `--validate-config` prints runtime validation warnings and errors.
- `--validate-schema` checks the JSON file against the schema files (to check if all fields are valid).

For first time debugging, use both schema and runtime validation. Note 
that a valid schema file can still be physically questionable, for example 
because the timestep is too large or all ions start outside the useful geometry.

## Reproducibility

Set the random seed explicitly when comparing runs:

```bash
icarion --seed 42 examples/ims/ims_basic.json
```

This overrides `simulation.rng_seed` from the configuration. For stochastic
collisions and reactions, keep the seed fixed while changing one physical
parameter at a time.

## Output control

Override the output location from the command line:

```bash
icarion --output-dir ./results/test1 --output ims_test1.h5 examples/ims/ims_basic.json
```

Useful flags:

- `--output FILE` or `-o FILE`: override the HDF5 trajectory filename.
- `--output-dir DIR`: override `output.folder`.
- `--buffer-byte-cap BYTES`: cap the in-memory trajectory buffer; `0` disables
  the cap.

For anything other than short test runs, prefer an explicit output directory and a descriptive filename.
See [Output files](output-files.md) for the HDF5 layout.

## Logging

```bash
icarion --verbose --log-file debug.log examples/ims/ims_basic.json
icarion --log-level DEBUG --log-format json --log-file run.jsonl config.json
```

Supported log levels are `DEBUG`, `INFO`, `WARN`, and `ERROR`. Supported log
formats are `text` and `json`.

## Runtime overrides

`--set KEY=VALUE` changes selected configuration values after loading the JSON
file. It can be passed multiple times.

```bash
icarion \
  --set simulation.dt_s=1e-10 \
  --set simulation.rng_seed=999 \
  --set physics.collision_model=EHSS \
  --set output.folder=./results/ehss \
  examples/ims/ims_basic.json
```

Supported high-value overrides include:

- `simulation.dt_s`, `simulation.total_time_s`, `simulation.write_interval`
- `simulation.rng_seed`, `simulation.integrator`, `simulation.rk45_min_step_s`
- `simulation.enable_gpu`, `simulation.enable_openmp`
- `physics.collision_model`
- `physics.enable_reactions`
- `physics.enable_space_charge`, `physics.enable_space_charge_gpu`
- `physics.enable_ou_thermalization`
- `output.folder`, `output.trajectory_file`, `output.trajectory_mode`
- `output.print_progress`, `output.buffer_byte_cap`
- `species_database` or `database.species`
- `reaction_database` or `database.reactions`

Unknown keys are ignored with a warning. For reproducible studies, record the
full command line together with the base configuration file or use adapted configuration files without CLI overrides.

## Performance flags

```bash
icarion --threads 8 --benchmark config.json
icarion --profile --profile-output profile.json config.json
```

- `--threads N` sets the OpenMP thread count to N for CPU builds with OpenMP.
- `--benchmark` prints a timing summary after the run.
- `--profile` enables profiling instrumentation.
- `--profile-output FILE` writes profiling data to JSON or CSV.

CPU scaling is simulation-dependent. For many CPU runs, a moderate thread
count is more efficient than using every hardware thread.

## Information flags

These flags print information and exit:

```bash
icarion --dump-build-info
icarion --dump-hdf5-schema
icarion --dump-config-schema
icarion --list-collision-models
icarion --list-integrators
icarion --check-deps
```

Use `--dump-build-info` when reporting a problem. It shows the version, git hash,
compiler/build mode, GPU status, OpenMP status, and dependency information.

## Exit codes

ICARION returns `0` for success and `1` for invalid arguments, validation
failures, or runtime errors. Shell scripts should check the exit code before
starting post-processing.

```bash
icarion config.json
if [ $? -ne 0 ]; then
  echo "ICARION run failed"
  exit 1
fi
```

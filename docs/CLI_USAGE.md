# ICARION Command-Line Interface (CLI) Guide

**Version:** 1.0.0  
**Last Updated:** December 2025

---

## Table of Contents

- [Quick Start](#quick-start)
- [Basic Usage](#basic-usage)
- [Core Options](#core-options)
- [Logging Options](#logging-options)
- [Output Control](#output-control)
- [Advanced Configuration](#advanced-configuration)
- [Information Flags](#information-flags)
- [Examples](#examples)
- [Exit Codes](#exit-codes)

---

## Quick Start

```bash
# Run a simulation
icarion config.json

# Validate configuration without running
icarion --dry-run config.json

# Run with verbose logging
icarion --verbose config.json

# Override configuration parameters
icarion --set simulation.dt_s=1e-10 config.json
```

---

## Basic Usage

```
icarion [OPTIONS] <config.json>
```

The configuration file is **required** for all simulation runs. It can be specified:
- As a positional argument: `icarion config.json`
- Using the flag: `icarion --config config.json`

---

## Core Options

### `--help`, `-h`
Show comprehensive help message with all available options.

```bash
icarion --help
```

### `--version`, `-v`
Display version information, git commit hash, and build type.

```bash
icarion --version
# Output:
# ICARION v1.0.0
# Git commit: 4a2cbee
# Build type: Full (with FieldSolver & Optimizer)
# GPU acceleration: Disabled
```

### `--config <FILE>`
Specify the configuration file (alternative to positional argument).

```bash
icarion --config examples/ims/ims_basic.json
```

### `--seed <N>`
Override the random number generator seed from the configuration file.  
Useful for reproducibility or running multiple independent simulations.

```bash
icarion --seed 42 config.json
```

### `--dry-run`
Validate the configuration file without actually running the simulation.  
Checks for syntax errors, missing fields, and parameter consistency.

```bash
icarion --dry-run config.json
```

### `--validate-config`
Similar to `--dry-run`, but provides more detailed validation output  
including warnings and errors.

```bash
icarion --validate-config config.json
```

### `--no-reactions`
Disable chemical reactions during the simulation.  
Only ion-neutral collisions will be computed.

```bash
icarion --no-reactions config.json
```

---

## Logging Options

### `--log-level <LEVEL>`
Set the logging verbosity level.

**Available levels:**
- `DEBUG` - Detailed debugging information
- `INFO` - General informational messages (default)
- `WARN` - Warning messages only
- `ERROR` - Error messages only

```bash
icarion --log-level DEBUG config.json
```

### `--verbose`
Enable verbose output (equivalent to `--log-level DEBUG`).

```bash
icarion --verbose config.json
```

### `--log-file <FILE>`
Write log output to a file instead of the console.

```bash
icarion --log-file simulation.log config.json
```

### `--log-format <FORMAT>`

Set the log output format. Available formats:

- `text` - Human-readable text output (default)
- `json` - Machine-readable JSON format for automated analysis

```bash
# JSON format for pandas/automated analysis
icarion --log-format json config.json

# Text format (default)
icarion --log-format text config.json
```

See [JSON_LOGGING.md](JSON_LOGGING.md) for analysis examples with pandas.

**Combining logging options:**
```bash
icarion --log-level DEBUG --log-file debug.log config.json
icarion --log-format json --log-file simulation.log config.json
```

---

## Output Control

### `--output <FILE>`, `-o <FILE>`
Override the output HDF5 filename specified in the configuration.

```bash
icarion --output results_run1.h5 config.json
```

### `--output-dir <DIR>`
Override the output directory specified in the configuration.

```bash
icarion --output-dir ./results/experiment1 config.json
```

### `--buffer-byte-cap <BYTES>`
Cap the in-memory trajectory buffer to avoid OOM; `0` (default) disables the cap.

```bash
# Cap at 512 MB
icarion --buffer-byte-cap 536870912 config.json
```

---

## Performance Options

### `--threads <N>`
Set the number of OpenMP threads for CPU parallelization.  
Overrides the `OMP_NUM_THREADS` environment variable.

**Usage:**
```bash
# Use 8 CPU threads
icarion --threads 8 config.json

# Single-threaded execution
icarion --threads 1 config.json

# Optimal for most systems (4-8 threads)
icarion --threads 4 config.json
```

**Performance Notes:**
- For CPU-only simulations, 4-8 threads often give best efficiency; beyond that, bandwidth can limit scaling.
- Ignored if OpenMP is not enabled at compile time.

### `--benchmark`
Enable profiling and print a timing summary after the run (uses the internal profiler).

```bash
icarion --benchmark config.json
```

**Output includes:**
- Timing summary for instrumented phases (ion generation, physics setup, main loop, output)
- Total simulation wall-clock time

### `--profile`
Enable profiling instrumentation for detailed performance analysis (same profiler as `--benchmark`).

```bash
icarion --profile config.json
```

### `--profile-output <FILE>`
Write profiling data to a file (JSON or CSV format).

```bash
# JSON format (default)
icarion --profile --profile-output profile.json config.json

# CSV format for Excel/pandas
icarion --profile --profile-output profile.csv config.json
```

---

## Advanced Configuration

### `--set <KEY=VALUE>`
Override individual configuration parameters from the command line.  
Can be specified multiple times to override multiple values.

**Supported keys:**

**Simulation parameters:**
- `simulation.dt_s` - Timestep [seconds]
- `simulation.total_time_s` - Total simulation time [seconds]
- `simulation.write_interval` - Output write interval
- `simulation.rng_seed` - Random number seed
- `simulation.integrator` - Integration method (RK4, RK45, Boris)
- `simulation.rk45_min_step_s` - Absolute minimum step size for RK45 [seconds]
- `simulation.enable_gpu` - Enable GPU acceleration (true/false; requires CUDA build, otherwise CPU fallback)
- `simulation.enable_openmp` - Enable OpenMP parallelization (true/false; requires OpenMP build)

**Physics parameters:**
- `physics.collision_model` - Collision model (NoCollisions, HardSphere/HSD, Langevin, Friction, HSS, EHSS, InteractionPotentialModel)
- `physics.enable_reactions` - Enable chemical reactions (true/false)
- `physics.enable_space_charge` - Enable space charge effects (true/false)
- `physics.enable_space_charge_gpu` - Prefer GPU-based P³M space charge (if built with CUDA; falls back to CPU if unavailable)
- `physics.space_charge_model` - Space charge backend (`auto`, `direct`, `grid`, `gpu`)
- `physics.enable_ou_thermalization` - Enable Ornstein-Uhlenbeck thermalization for deterministic collision models (requires custom wiring of `gamma_for_ou`; incompatible with HSS/EHSS)
- `physics.collision_subcycles_per_step` - Split stochastic collision updates into this many micro-steps
- `physics.collision_multi_event_mode` - Enable approximate multi-collision micro-subcycling (true/false)
- `physics.collision_max_events_per_step` - Minimum subcycle count used by multi-event mode
- `physics.ipm_orientation_mode` - InteractionPotential orientation sampling (`random` or `fixed`)
- `physics.ipm_fixed_orientation_index` - Fixed IPM orientation index when mode is `fixed`
- `physics.ipm_vrel_log_prefix` - Optional InteractionPotential relative-speed CSV log prefix
- `physics.ipm_momentum_log_prefix` - Optional InteractionPotential momentum CSV log prefix

**Output parameters:**
- `output.folder` - Output directory path
- `output.trajectory_file` - Trajectory filename
- `output.trajectory_mode` - Output mode (`full` or `minimal`)
- `output.print_progress` - Print progress updates (true/false)
- `output.buffer_byte_cap` - Optional trajectory buffer memory cap in bytes (`0` disables)

**Database paths:**
- `species_database` (or `database.species`) - Path to species database
- `reaction_database` (or `database.reactions`) - Path to reaction database

**Examples:**
```bash
# Change timestep
icarion --set simulation.dt_s=1e-10 config.json

# Use different collision model
icarion --set physics.collision_model=Langevin config.json

# Multiple overrides
icarion --set simulation.rng_seed=999 \
        --set physics.collision_model=EHSS \
        --set output.folder=./results \
        config.json
```

---

## Information Flags

These flags provide information and exit immediately. `--validate-schema` requires a configuration file; the others do not.

### `--dump-build-info`
Show detailed build configuration and features.

```bash
icarion --dump-build-info
```

**Output includes:**
- ICARION version and git commit hash
- Build date and compiler version
- C++ standard
- GPU acceleration status
- OpenMP status
- Build mode (Release/Debug)
- Dependencies and their versions
- Schema and example paths

### `--dump-hdf5-schema`
Display the HDF5 output file schema documentation.

```bash
icarion --dump-hdf5-schema
```

Shows the structure of output HDF5 files including:
- `/metadata/` - Config, reproducibility, system, species, reactions, completion
- `/trajectory/` - Ion trajectory data
- `/ions/` - Ion initial conditions
- `/domains/` - Per-domain geometry, environment, and fields

### `--dump-config-schema`
Show available JSON configuration schemas.

```bash
icarion --dump-config-schema
```

Lists all schema files in `./schema/` and validation commands.

### `--list-collision-models`
List available collision models with descriptions.

```bash
icarion --list-collision-models
```

**Available models:**
- **NoCollisions** - No collision modeling (free flight)
- **HardSphere** - Hard sphere collisions (deterministic; same as HSD)
- **Langevin** - Langevin collisions (polarization)
- **Friction** - Frictional drag
- **HSS** - Hard sphere stochastic
- **EHSS** - Enhanced hard sphere statistical

### `--list-integrators`
List available numerical integrators with characteristics.

```bash
icarion --list-integrators
```

**Available integrators:**
- **RK4** - 4th-order Runge-Kutta (fixed timestep, fast, deterministic)
- **RK45** - Adaptive Runge-Kutta (variable timestep, error control)
- **Boris** - Boris pusher (specialized for E×B fields)

### `--validate-schema`
Validate a configuration file against the JSON schema.

```bash
icarion --validate-schema config.json
```

Runs `schema/validator.py` to check:
- JSON syntax correctness
- Required fields presence
- Type consistency
- Value constraints

### `--check-deps`
Verify all dependencies and display their versions.

```bash
icarion --check-deps
```

**Output includes:**
- Core libraries (HDF5, JsonCpp, cxxopts)
- Optional features (CUDA, OpenMP)
- Build configuration details

---

## Examples

### Basic Simulation

```bash
# Run with default settings
icarion examples/ims/ims_basic.json
```

### Validation Only

```bash
# Check configuration without running
icarion --dry-run examples/tof/tof_basic.json

# Detailed validation with schema check
icarion --validate-schema examples/orbitrap/orbitrap_basic.json
```

### Debugging

```bash
# Verbose logging to file
icarion --verbose --log-file debug.log config.json

# Validate and show detailed info
icarion --validate-config --log-level DEBUG config.json
```

### Parameter Studies

```bash
# Run multiple simulations with different seeds
for seed in {1..10}; do
    icarion --seed $seed \
            --output results_seed${seed}.h5 \
            config.json
done

# CPU thread scaling benchmark
for threads in 1 2 4 8 16; do
    icarion --threads $threads \
            --benchmark \
            --output results_threads${threads}.h5 \
            config.json
done

# Test different collision models
for model in NoCollisions HSD Langevin EHSS; do
    icarion --set physics.collision_model=$model \
            --output results_${model}.h5 \
            config.json
done

# Timestep convergence study
for dt in 1e-9 1e-10 1e-11; do
    icarion --set simulation.dt_s=$dt \
            --output results_dt${dt}.h5 \
            config.json
done
```

### Custom Output

```bash
# Custom output location
icarion --output-dir ./results/experiment1 \
        --output trajectory_run1.h5 \
        config.json

# Disable reactions
icarion --no-reactions \
        --set output.trajectory_file=no_reactions.h5 \
        config.json
```

### Combining Options

```bash
# Full featured run with performance optimization
icarion --threads 8 \
        --benchmark \
        --verbose \
        --log-file simulation.log \
        --seed 42 \
        --set simulation.dt_s=5e-11 \
        --set physics.collision_model=Langevin \
        --output-dir ./results \
        --output high_res_trajectory.h5 \
        config.json
```

---

## Exit Codes

ICARION uses standard exit codes:

| Code | Meaning |
|------|---------|
| `0` | Success - simulation completed normally |
| `1` | Error - invalid arguments, missing config, validation failure, or runtime error |

**Checking exit codes in scripts:**

```bash
#!/bin/bash
icarion config.json
if [ $? -eq 0 ]; then
    echo "Simulation successful"
    python analyze_results.py
else
    echo "Simulation failed"
    exit 1
fi
```

---

## Tips & Best Practices

### 1. **Always validate first**
```bash
icarion --validate-config config.json && icarion config.json
```

### 2. **Use meaningful output names**
```bash
icarion --output "run_$(date +%Y%m%d_%H%M%S).h5" config.json
```

### 3. **Keep logs for long simulations**
```bash
icarion --log-file simulation_$(date +%Y%m%d).log config.json
```

### 4. **Document parameter overrides**
```bash
# Save command to reproduce results
echo "icarion --seed 42 --set simulation.dt_s=1e-10 config.json" > run_info.txt
icarion --seed 42 --set simulation.dt_s=1e-10 config.json
```

### 5. **Check build info before reporting issues**
```bash
icarion --dump-build-info > build_info.txt
```

---

## See Also

- [Configuration Guide](CONFIG_GUIDE.md)
- [HDF5 Output Structure](HDF5_OUTPUT_STRUCTURE.md)
- [Developer's Guide](DEVELOPERS_GUIDE.md)
- [Architecture Overview](ARCHITECTURE.md)

---

## Getting Help

- **Command-line help:** `icarion --help`
- **Build information:** `icarion --dump-build-info`
- **Schema documentation:** `icarion --dump-config-schema`
- **GitHub Issues:** https://github.com/ICARION-Project/ICARION/issues
- **Documentation:** https://icarion.readthedocs.io (coming soon)

---

**Last Modified:** December 2025  
**ICARION Version:** 1.0.0

# ICARION Command-Line Interface (CLI) Guide

**Version:** 1.0  
**Last Updated:** November 21, 2025

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
# ICARION v1.0.0 (commit: 4a2cbee)
# Build: Release | GPU: Disabled | OpenMP: Disabled
```

### `--config <FILE>`
Specify the configuration file (alternative to positional argument).

```bash
icarion --config examples/ims_basic.json
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
including warnings and suggestions.

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

- `text` - Human-readable colored console output (default)
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
- `simulation.enable_gpu` - Enable GPU acceleration (true/false)
- `simulation.enable_openmp` - Enable OpenMP parallelization (true/false)

**Physics parameters:**
- `physics.collision_model` - Collision model (NoCollisions, HardSphere, Langevin, Friction, HSS, EHSS)
- `physics.enable_reactions` - Enable chemical reactions (true/false)
- `physics.enable_space_charge` - Enable space charge effects (true/false)
- `physics.enable_ou_thermalization` - Enable Ornstein-Uhlenbeck thermalization (true/false)

**Output parameters:**
- `output.folder` - Output directory path
- `output.trajectory_file` - Trajectory filename
- `output.print_progress` - Print progress updates (true/false)

**Database paths:**
- `species_database_path` - Path to species database
- `reaction_database_path` - Path to reaction database

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

### Future Options (Planned for v1.1)

#### `--benchmark` *(Coming Soon)*
Print detailed timing statistics for each simulation phase.
```bash
icarion --benchmark config.json
```

#### `--profile` *(Coming Soon)*
Enable profiling instrumentation (requires profiler build).
```bash
icarion --profile config.json
```

#### `--check-nan` *(Coming Soon)*
Enable NaN/Inf checks in the integrator for debugging numerical issues.
```bash
icarion --check-nan config.json
```

---

## Information Flags

These flags provide information about ICARION and exit immediately  
without requiring a configuration file.

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
- `/metadata/` - Simulation metadata
- `/trajectory/` - Ion trajectory data
- `/species/` - Species properties
- `/config/` - Complete configuration snapshot

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
- **HardSphere** - Hard sphere collisions (deterministic)
- **Langevin** - Langevin collisions (polarization)
- **Friction** - Frictional drag
- **HSS** - Hard sphere stochastical
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
icarion examples/ims_basic.json
```

### Validation Only

```bash
# Check configuration without running
icarion --dry-run examples/tof_basic.json

# Detailed validation with schema check
icarion --validate-schema examples/orbitrap_basic.json
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

# Test different collision models
for model in NoCollisions HardSphere Langevin EHSS; do
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
# Full featured run
icarion --verbose \
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
| `1` | Error - invalid arguments, missing config, validation failure |
| `2` | Runtime error - simulation crashed, numerical issues |

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
- [HDF5 Output Schema](OUTPUT_SCHEMA.md)
- [Developer's Guide](DEVELOPERS_GUIDE.md)
- [Architecture Overview](ARCHITECTURE.md)

---

## Getting Help

- **Command-line help:** `icarion --help`
- **Build information:** `icarion --dump-build-info`
- **Schema documentation:** `icarion --dump-config-schema`
- **GitHub Issues:** https://github.com/icarion-project/icarion/issues
- **Documentation:** https://icarion.readthedocs.io (coming soon)

---

**Last Modified:** November 21, 2025  
**ICARION Version:** 1.0.0


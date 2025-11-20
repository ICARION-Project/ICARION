# ICARION Configuration Guide

This guide explains how to create and structure JSON configuration files for ICARION simulations.

## Table of Contents

- [Quick Start](#quick-start)
- [Configuration Structure](#configuration-structure)
- [Required Fields](#required-fields)
- [Schema Validation](#schema-validation)
- [Helper Tools](#helper-tools)
- [Examples](#examples)

---

## Quick Start

### Minimal Configuration

```json
{
  "simulation": {
    "total_time_s": 1e-3,
    "dt_s": 1e-9,
    "integrator": "RK4"
  },
  "physics": {
    "collision_model": "NoCollisions"
  },
  "output": {
    "folder": "./results",
    "trajectory_file": "trajectories.h5"
  },
  "domains": [
    {
      "name": "my_domain",
      "instrument": "IMS",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.0],
        "length_m": 0.05,
        "radius_m": 0.01
      },
      "env": {
        "pressure_Pa": 101325.0,
        "temperature_K": 300.0,
        "gas_species": "He"
      }
    }
  ]
}
```

### Creating a Config with the Helper Tool

```bash
# Interactive mode
python3 scripts/create_config.py

# Command-line mode
python3 scripts/create_config.py \
  --template ims \
  --output my_config.json \
  --time 1e-3 \
  --timestep 1e-9
```

---

## Configuration Structure

ICARION uses a hierarchical JSON configuration with four main sections:

```json
{
  "simulation": { ... },    // Time parameters, integrator, execution mode
  "physics": { ... },       // Collision model, reactions, space charge
  "output": { ... },        // Output paths and options
  "domains": [ ... ]        // Array of instrument domains
}
```

### Optional Top-Level Fields

```json
{
  "title": "My Simulation",                          // Human-readable description
  "species_database_path": "data/species.json",      // Path to species database
  "reaction_database_path": "data/reactions.json",   // Path to reaction database
  "ion_cloud_path": "data/initial_cloud.h5"          // Path to initial ion distribution
}
```

---

## Required Fields

### Simulation Section

```json
"simulation": {
  "total_time_s": 1e-3,      // REQUIRED: Total simulation time [seconds]
  "dt_s": 1e-9,              // REQUIRED: Time step [seconds]
  "integrator": "RK4",       // Default integrator (can be overridden per-domain)
  "write_interval": 100,     // Steps between trajectory snapshots
  "enable_gpu": false,       // Enable GPU acceleration
  "enable_openmp": false,    // Enable OpenMP threading
  "rng_seed": 42             // Random number generator seed
}
```

**Available Integrators:**
- `"RK4"` - 4th order Runge-Kutta (default, stable)
- `"RK45"` - Runge-Kutta-Fehlberg (adaptive step size)
- `"Boris"` - Boris integrator (for strong magnetic fields)

### Physics Section

```json
"physics": {
  "collision_model": "HSMC",           // Collision model
  "enable_ou_thermalization": false,   // Enable Ornstein-Uhlenbeck thermalization
  "force_ou_for_stochastic": false     // Force OU for stochastic models
}
```

**Available Collision Models:**

- `"NoCollisions"` - Free flight (no collisions)

- `"HSD"` - Hard Sphere Deterministic (collision frequency does not depend on velocity)
- `"Friction"` - Frictional-drag model (depending on mobility)
- `"Langevin"` - Langevin dynamics (long-range and velocity-dependent collision frequency)

- `"HSS"` - Hard Sphere Stochastic  
- `"EHSS"` - Exact Hard Sphere Scattering (recommended for low pressures)


### Output Section

```json
"output": {
  "folder": "./results",              // Output directory (relative or absolute)
  "trajectory_file": "trajectories.h5", // HDF5 trajectory file name
  "print_progress": true              // Print progress to console
}
```

### Domains Section

Each domain represents a physical region with specific instrument geometry and conditions:

```json
"domains": [
  {
    "name": "drift_region",           // REQUIRED: Domain identifier
    "instrument": "IMS",              // REQUIRED: Instrument type
    "integrator": "RK4",              // OPTIONAL: Override global integrator
    
    "geometry": {                     // REQUIRED
      "origin_m": [0.0, 0.0, 0.0],    // Starting position [x, y, z] in meters
      "length_m": 0.05,               // Domain length along z-axis
      "radius_m": 0.01                // Cylindrical radius (for IMS, TOF, etc.)
    },
    
    "env": {                          // REQUIRED
      "pressure_Pa": 101325.0,        // Gas pressure [Pascal]
      "temperature_K": 300.0,         // Gas temperature [Kelvin]
      "gas_species": "He",            // Gas species (He, N2, Ar, etc.)
      "gas_velocity_m_s": [0, 0, 0]   // Gas flow velocity [m/s]
    },
    
    "fields": {                       // OPTIONAL
      "DC": {
        "axial_V": 250.0,             // Axial DC voltage [V]
        "EN_Td": 10.0                 // Reduced field strength [Townsend]
      },
      "RF": {
        "voltage_V": 0.0,             // RF amplitude [V]
        "frequency_Hz": 1e6,          // RF frequency [Hz]
        "phase_rad": 0.0              // RF phase [radians]
      }
    }
  }
]
```

**Available Instrument Types:**

- `"IMS"` - Ion Mobility Spectrometry
- `"TOF"` - Time-of-Flight
- `"LQIT"` - Linear Quadrupole Ion Trap
- `"Orbitrap"` - Orbitrap mass analyzer
- `"Quadrupole"` - Quadrupole mass filter
- `"FTICR"` - Fourier Transform Ion Cyclotron Resonance
- `"NoFixedInstrument"` - Generic region

---

## Schema Validation

### JSON Schema Files

All configuration schemas are located in `src/core/config/schema/v1.0/`:

- **`icarion-config.schema.json`** - Master schema (top-level)
- **`simulation.schema.json`** - Simulation parameters
- **`physics.schema.json`** - Physics options
- **`output.schema.json`** - Output configuration
- **`domain.schema.json`** - Domain definition
- **`geometry.schema.json`** - Geometry specification
- **`environment.schema.json`** - Gas environment
- **`fields.schema.json`** - Electric/magnetic fields
- **`common-types.schema.json`** - Reusable type definitions

### Validating Your Config

#### Using Python Validator

```bash
# Validate a config file
python3 src/core/config/schema/validator.py my_config.json

# Validate multiple files
python3 src/core/config/schema/validator.py examples/*.json
```

#### Using ICARION CLI

```bash
# Validate config without running simulation
./build/src/icarion_main --validate-config my_config.json

# Dry-run (load config, initialize, but don't simulate)
./build/src/icarion_main --dry-run my_config.json
```

---

## Helper Tools

### Config Creation Script

Located at `scripts/create_config.py`:

```bash
# Interactive mode (guided prompts)
python3 scripts/create_config.py

# Template-based creation
python3 scripts/create_config.py --template ims --output ims_config.json

# Customize template
python3 scripts/create_config.py \
  --template tof \
  --time 5e-4 \
  --timestep 1e-9 \
  --pressure 1e-6 \
  --temperature 300 \
  --output my_tof.json

# List available templates
python3 scripts/create_config.py --list-templates
```

**Available Templates:**

- `ims` - Ion Mobility Spectrometry (drift tube)
- `tof` - Time-of-Flight mass spectrometer
- `lqit` - Linear Quadrupole Ion Trap
- `orbitrap` - Orbitrap mass analyzer
- `minimal` - Minimal valid configuration

---

## Examples

### Ion Mobility Spectrometry (IMS)

```json
{
  "simulation": {
    "total_time_s": 5e-4,
    "dt_s": 1e-8,
    "integrator": "RK4",
    "enable_gpu": true
  },
  "physics": {
    "collision_model": "HSS"
  },
  "output": {
    "folder": "./results/ims",
    "trajectory_file": "ims_trajectories.h5"
  },
  "domains": [
    {
      "name": "drift_region",
      "instrument": "IMS",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.0],
        "length_m": 0.05,
        "radius_m": 0.015
      },
      "env": {
        "pressure_Pa": 200.0,
        "temperature_K": 300.0,
        "gas_species": "He"
      },
      "fields": {
        "DC": {
          "EN_Td": 10.0
        }
      }
    }
  ]
}
```

### Time-of-Flight (TOF)

```json
{
  "simulation": {
    "total_time_s": 1e-4,
    "dt_s": 1e-10,
    "integrator": "RK4"
  },
  "physics": {
    "collision_model": "NoCollisions"
  },
  "output": {
    "folder": "./results/tof",
    "trajectory_file": "tof_trajectories.h5"
  },
  "domains": [
    {
      "name": "acceleration_region",
      "instrument": "TOF",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.0],
        "length_m": 0.02,
        "radius_m": 0.01
      },
      "env": {
        "pressure_Pa": 1e-6,
        "temperature_K": 300.0,
        "gas_species": "He"
      },
      "fields": {
        "DC": {
          "axial_V": 5000.0
        }
      }
    },
    {
      "name": "drift_region",
      "instrument": "TOF",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.02],
        "length_m": 0.5,
        "radius_m": 0.01
      },
      "env": {
        "pressure_Pa": 1e-6,
        "temperature_K": 300.0,
        "gas_species": "He"
      }
    }
  ]
}
```

### Multi-Domain Setup

```json
{
  "simulation": {
    "total_time_s": 1e-3,
    "dt_s": 1e-9,
    "integrator": "RK4"
  },
  "physics": {
    "collision_model": "HSS"
  },
  "output": {
    "folder": "./results/hybrid",
    "trajectory_file": "trajectories.h5"
  },
  "domains": [
    {
      "name": "ims_drift",
      "instrument": "IMS",
      "integrator": "RK4",
      "geometry": { "origin_m": [0, 0, 0], "length_m": 0.05, "radius_m": 0.015 },
      "env": { "pressure_Pa": 200.0, "temperature_K": 300.0, "gas_species": "He" },
      "fields": { "DC": { "EN_Td": 10.0 } }
    },
    {
      "name": "tof_region",
      "instrument": "TOF",
      "integrator": "RK45",
      "geometry": { "origin_m": [0, 0, 0.05], "length_m": 0.5, "radius_m": 0.02 },
      "env": { "pressure_Pa": 1e-6, "temperature_K": 300.0, "gas_species": "He" },
      "fields": { "DC": { "axial_V": 5000.0 } }
    }
  ]
}
```

---

## Common Patterns

### Per-Domain Integrator Override

Use the global `simulation.integrator` as default, but override for specific domains:

```json
{
  "simulation": {
    "integrator": "RK4"  // Global default
  },
  "domains": [
    {
      "name": "low_pressure_region",
      "integrator": "RK45",  // Override: use adaptive stepping
      "env": { "pressure_Pa": 1e-6, ... }
    },
    {
      "name": "high_pressure_region",
      // No integrator specified → uses "RK4" from simulation
      "env": { "pressure_Pa": 101325, ... }
    }
  ]
}
```

### Collision Model Selection

```json
// No collisions (vacuum, short flight times)
"physics": { "collision_model": "NoCollisions" }

// High pressure (IMS, drift tubes)
"physics": { "collision_model": "Friction" }

// High precision (research simulations)
"physics": { "collision_model": "EHSS" }
```

---

## Troubleshooting

### Common Errors

**Error: Missing required field 'dt_s' in simulation config**

- Solution: Add `"dt_s": 1e-9` to simulation section

**Error: Missing required field 'total_time_s' in simulation config**

- Solution: Add `"total_time_s": 1e-3` to simulation section

**Error: Domain missing required 'name' field**

- Solution: Add `"name": "my_domain"` to each domain

**Error: Configuration validation failed**

- Solution: Run `--validate-config` to see detailed error messages

### Best Practices

1. **Always validate configs** before running expensive simulations:

   ```bash
   ./build/src/icarion_main --validate-config my_config.json
   ```

2. **Start with a template** from `examples/` directory

3. **Use appropriate time steps**:

   - IMS (collisional): `dt_s = 1e-10` to `1e-9`
   - TOF (vacuum): `dt_s = 1e-9` to `1e-8`

4. **Check units**: All values in SI units (meters, seconds, Pascal, Kelvin)

5. **Enable GPU** for large simulations (`enable_gpu: true`)

---

## See Also

- **Schema Documentation**: `src/core/config/schema/README.md`
- **API Documentation**: `docs/PUBLIC_CPP_API_v1.0.md` #todo
- **CLI Flags**: `docs/CLI_INTERFACE_v1.0.md`#todo
- **Example Configs**: `examples/*.json`
- **Python Helper**: `scripts/create_config.py --help`

---

## Version

This guide corresponds to ICARION v1.0 configuration schema.

For schema version history, see `src/core/config/schema/CHANGELOG.md`.

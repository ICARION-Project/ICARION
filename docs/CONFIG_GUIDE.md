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
  "species_database": "data/species_database_v1.json",  // Species properties database (optional)
  "reaction_database": "data/reactions_database_v1.json", // Reaction rates database (optional)
  "ion_cloud": "data/initial_cloud.h5"               // Initial ion distribution
}
```

**Note on Database Paths:**

- Both `species_database` and `reaction_database` are **optional**
- If not specified, ICARION will automatically search for global fallback databases:
  - `data/species_database_v1.json` (global species database)
  - `data/reactions_database_v1.json` (global reactions database)
- The fallback search starts from the config file's directory and searches up to 5 levels
- A log message `ℹ No [species/reaction] database specified, using global fallback: ...` indicates fallback usage
- If no database is found (neither specified nor fallback), databases remain empty (which is valid for simulations without reactions)

---

## Species and Reaction Databases (v1.0)

ICARION supports external databases for species properties and reaction rates.

### Species Database

Species databases define physical properties of ions and neutrals in user-friendly units:

```json
{
  "species": {
    "H3O+": {
      "name": "Hydronium ion",
      "mass_amu": 19.02,
      "charge": 1,
      "mobility_cm2Vs": 2.8,
      "CCS_A2": 11.0,
      "reference_temperature_K": 300.0,
      "reference_pressure_Pa": 101325.0,
      "ccs_method": "trajectory"
    },
    "H2O": {
      "name": "Water",
      "mass_amu": 18.015,
      "charge": 0,
      "polarizability_A3": 1.47
    }
  }
}
```

**Required fields:**

- `mass_amu`: Molecular mass in atomic mass units
- `charge`: Charge state (integer: -2, -1, 0, +1, +2, ...)

**Optional fields (ions):**

- `mobility_cm2Vs`: Reduced ion mobility at STP [cm²/(V·s)]
- `CCS_A2`: Collision cross-section [Ų]

**Optional fields (neutrals):**

- `polarizability_A3`: Electric polarizability [ų]

**Optional metadata:**

- `name`: Human-readable name
- `geometry_file`: Path to molecular geometry (e.g., .xyz, .pdb, required for EHSS)
- `reference_temperature_K`: Temperature for mobility/CCS measurements (just for reproducibility)
- `reference_pressure_Pa`: Pressure for mobility measurements (just for reproducibility)
- `ccs_method`: Method used to determine CCS (`"trajectory"`, `"projection"`, etc., for reproducibility)

**Unit conversions:** ICARION automatically converts to SI internally:

- `mass_amu` → kg (× 1.66054e-27)
- `charge` → C (× 1.60218e-19)
- `mobility_cm2Vs` → m²/(V·s) (× 1e-4)
- `CCS_A2` → m² (× 1e-20)
- `polarizability_A3` → m³ (× 1e-30)

### Reaction Database

Reaction databases define ion-molecule reactions with concentration-dependent rates:

```json
{
  "reactions": [
    {
      "id": "rxn_001_h3o_to_h5o2",
      "reactant": "H3O+",
      "product": "H5O2+",
      "rate_constant_m3s": 3.5e-9,
      "order": [
        {
          "species": "H2O",
          "exponent": 1,
          "concentration_m3": 2.5e25
        }
      ],
      "description": "Proton transfer from H3O+ to water cluster",
      "reference": "Smith et al., J. Chem. Phys. (2020)",
      "temperature_K": 300.0
    }
  ]
}
```

**Required fields:**

- `id`: Unique reaction identifier
- `reactant`: Species ID (must exist in species database)
- `product`: Species ID (must exist in species database)
- `rate_constant_m3s`: Base rate constant [m³/s]

**Optional fields:**

- `order`: Array of concentration-dependent terms
  - `species`: Species ID for concentration dependence
  - `exponent`: Concentration exponent (allowed: 0, 1, or 2)
  - `concentration_m3`: Fixed concentration [m⁻³] for pseudo-first-order
- `description`: Human-readable description
- `reference`: Literature reference
- `temperature_K`: Temperature at which rate was measured

**Effective rate calculation:**

The effective rate includes concentration dependencies:
```
k_eff = k * [S₁]^n₁ * [S₂]^n₂ * ...
```

### Using Databases in Configs

**Enable species database:**

```json
{
  "species_database": "data/species_database_v1.json",
  "domains": [ ... ]
}
```

**Enable reactions:**

```json
{
  "species_database": "data/species_database_v1.json",
  "reaction_database": "data/reactions_database_v1.json",
  "physics": {
    "enable_reactions": true
  }
}
```

**Path resolution:** Paths can be:

- Absolute: `/absolute/path/to/database.json`
- Relative to config file: `../data/species.json`
- Relative to working directory: `data/species.json`

**Schema validation:** Use provided JSON schemas to validate databases:

```bash
python3 schema/validate_schema.py schema/species.schema.json data/species_database_v1.json
python3 schema/validate_schema.py schema/reactions.schema.json data/reactions_database_v1.json
```

**Example configs:**

- `examples/ims_with_species_db.json` - IMS with species database
- `examples/reaction_demo.json` - Ion-molecule reactions

---

## Ion Initialization (v1.0)

ICARION provides flexible ion initialization with per-species position boundaries and velocity distributions.

### Configuration Syntax

Ion configuration can be specified in two ways:

#### Option 1: Generate from Species List

Each species can have its own spatial region (boundaries) and velocity distribution:

```json
{
  "ions": {
    "species": [
      {
        "id": "H3O+",
        "count": 1000,
        "position": {
          "type": "gaussian",
          "center": [0.0, 0.0, 0.0],
          "std": [0.001, 0.001, 0.001]
        },
        "velocity": {
          "type": "thermal",
          "temperature_K": 300
        }
      }
    ]
  }
}
```

#### Option 2: Load from Ion Cloud File

```json
{
  "ions": {
    "from_file": "ion_clouds/my_cloud.json"
  }
}
```

### Position Distributions

Each species can be created in its own spatial region (boundaries):

#### Point (Single Location)

All ions at same position:

```json
"position": {
  "type": "point",
  "center": [0.0, 0.0, 0.0]  // [m]
}
```

#### Gaussian (Normal Distribution)

3D Gaussian distribution with independent standard deviations per axis:

```json
"position": {
  "type": "gaussian",
  "center": [0.0, 0.0, 0.0],         // Mean position [m]
  "std": [0.001, 0.002, 0.003]       // Std dev per axis [m]
}
```

#### Uniform Sphere

Uniform distribution within a sphere:

```json
"position": {
  "type": "uniform_sphere",
  "center": [0.0, 0.0, 0.0],   // Sphere center [m]
  "radius": 0.005              // Sphere radius [m]
}
```

#### Uniform Cylinder

Uniform distribution within a cylinder (aligned with z-axis):

```json
"position": {
  "type": "uniform_cylinder",
  "center": [0.0, 0.0, 0.0],   // Cylinder center [m]
  "radius": 0.002,             // Radial extent [m]
  "length": 0.01               // Axial extent (z-direction) [m]
}
```

#### Uniform Box

Uniform distribution within a rectangular box:

```json
"position": {
  "type": "uniform_box",
  "min": [-0.001, -0.001, 0.0],    // Min corner [m]
  "max": [0.001, 0.001, 0.01]      // Max corner [m]
}
```

### Velocity Distributions

#### Fixed (All Ions Same Velocity)

```json
"velocity": {
  "type": "fixed",
  "value": [0.0, 0.0, 100.0]  // Velocity [m/s]
}
```

#### Thermal (Maxwell-Boltzmann)

Random directions sampled from Maxwell-Boltzmann distribution. **No directed drift** (always random):

```json
"velocity": {
  "type": "thermal",
  "temperature_K": 300  // Temperature [K]
}
```

**Note:** For thermal distribution, velocity directions are **always random**. Use `kinetic` for directed velocities.

#### Kinetic (Directed with Energy)

Fixed kinetic energy in a specified direction, with optional angular spread:

```json
"velocity": {
  "type": "kinetic",
  "energy_eV": 0.1,              // Kinetic energy [eV]
  "direction": [0, 0, 1],        // Direction vector (normalized automatically)
  "spread_angle_deg": 10         // Optional: cone angle spread [degrees]
}
```

#### Gaussian

Gaussian distribution for each velocity component:

```json
"velocity": {
  "type": "gaussian",
  "mean": [0.0, 0.0, 100.0],     // Mean velocity [m/s]
  "std": [10.0, 10.0, 20.0]      // Std dev per axis [m/s]
}
```

### Multi-Species Example

Different species in different spatial regions with different velocities:

```json
{
  "species_database": "data/species_database_v1.json",
  "ions": {
    "species": [
      {
        "id": "H3O+",
        "count": 500,
        "position": {
          "type": "gaussian",
          "center": [0.0, 0.0, 0.0],
          "std": [0.001, 0.001, 0.001]
        },
        "velocity": {
          "type": "thermal",
          "temperature_K": 300
        }
      },
      {
        "id": "H5O2+",
        "count": 500,
        "position": {
          "type": "uniform_sphere",
          "center": [0.005, 0.0, 0.0],
          "radius": 0.002
        },
        "velocity": {
          "type": "thermal",
          "temperature_K": 350
        }
      },
      {
        "id": "O2-",
        "count": 100,
        "position": {
          "type": "uniform_cylinder",
          "center": [0.0, 0.0, 0.01],
          "radius": 0.001,
          "length": 0.005
        },
        "velocity": {
          "type": "kinetic",
          "energy_eV": 0.1,
          "direction": [0, 0, 1],
          "spread_angle_deg": 10
        }
      }
    ]
  }
}
```

### Ion Cloud File Format

For loading ions from file:

```json
{
  "ions": [
    {
      "species": "H3O+",
      "pos": [0.0, 0.0, 0.0],      // Position [m]
      "vel": [0.0, 0.0, 100.0],    // Velocity [m/s]
      "birth_time": 0.0            // Birth time [s]
    },
    {
      "species": "H5O2+",
      "pos": [0.001, 0.0, 0.0],
      "vel": [10.0, 0.0, 50.0],
      "birth_time": 0.0
    }
  ]
}
```

**Note:** Species must exist in species database. Properties (mass, charge, mobility, CCS) are looked up automatically.

### Validation

Ion configuration is validated during loading:

- Species IDs must exist in species database
- Position/velocity parameters must be physically reasonable
- Warnings for missing species or invalid configurations

Validation errors are collected in `ValidationResult` (not thrown as exceptions).

**Example validation output:**

```
✓ Generated 1100 ions from configuration
  - 500 H3O+ ions (gaussian position, thermal velocity)
  - 500 H5O2+ ions (sphere position, thermal velocity)
  - 100 O2- ions (cylinder position, kinetic velocity)
```

**Example configs:**

- `examples/ion_clouds/multi_species_example.json` - Multi-species with different boundaries

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

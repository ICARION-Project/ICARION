ICARION v1.0 — Ion Collision And Reaction IntegratiON

A modular C++/CUDA framework for multi-domain ion dynamics simulation.

---

<p align="center"> <img src="https://img.shields.io/badge/C++17-HPC-blue?style=flat-square"> <img src="https://img.shields.io/badge/CUDA-optional-green?style=flat-square"> <img src="https://img.shields.io/badge/HDF5-output-important?style=flat-square"> <img src="https://img.shields.io/badge/License-Apache--2.0-orange?style=flat-square"> </p>

---

# Overview

**ICARION** is a high-performance ion trajectory simulation engine for:

- drift tube IMS
- RF quadrupoles
- Orbitraps
- TOF analyzers
- linear quadrupole ion traps (LQIT)
- multi-domain instruments (IMS → Quadrupole → TOF, …)

It supports:

- Stochastic & deterministic collision models (EHSS, HSS, Friction, Langevin, HardSphere), see docs/COLLISION_MODELS.md
- Ion–neutral reactions with Arrhenius rates
- Modular integrators (fixed-step RK4, adaptive RK45, Boris pusher)
- GPU acceleration via CUDA
- HDF5 output with metadata for full reproducibility
- Strict input validation for species, reactions, domains, and global parameters

The architecture is designed for reproducible science, publication-grade accuracy, and future extensibility (FieldSolver, Optimizer, multi-physics).

# Key Features

## Physics

- Deterministic & stochastic collisions (see docs/COLLISION_MODELS.md)
- Ion–neutral reactions with temperature-dependent rates
- Mason–Schamp CCS/mobility conversion
- Domain-dependent gas environments
- DC, RF, AC, magnetic fields (all electric fields either analytical or initialized from precomputed field arrays)

## Instruments

- IMS
- Quadrupole (RF)
- Orbitrap (hyperlogarithmic geometry)
- TOF
- LQIT with AC excitation
- FT-ICR

## Numerics

- fixed-step RK4, adaptive RK45, Boris pusher
- Deterministic + stochastic collision loop
- OpenMP support (multi-core CPU)
- GPU support (enable_gpu = true)

## Input System

4-file SSOT architecture:

```pgsql
species.json
reactions.json
config.json
molecular_geometry.json
```

## Output

- Open-source HDF5 format with hierarchical structure
- Chunked, extendable datasets
- Species metadata only once (reduces file size)
- Full config and system metadata snapshot stored in output file for reproducibility

---

# Building ICARION

## Requirements

- C++17 compiler (GCC, Clang, MSVC)
- CMake ≥ 3.16
- HDF5 ≥ 1.10
- (Optional) CUDA ≥ 11

**Build**

```bash
git clone https://github.com/ICARION-Project/ICARION.git
cd ICARION
mkdir build && cd build
cmake ..
make -j8
```

**Run**

```bash
./icarion input/config.json
```

---

# Input Files

ICARION uses a 4-file Single-Source-of-Truth (SSOT) input architecture:

1. species.json (default: `data/species_database_v1.json`)  
   Species database with ion properties (mass, charge, mobility, CCS, etc.)

Species database loaded via `SpeciesConfig`
(real implementation: `SpeciesLoader.cpp`)

Example:

```json
{
  "species": {
    "H3O+": {
      "mass_amu": 19.02,
      "charge": 1,
      "mobility_cm2Vs": 24.1,
      "CCS_A2": 24.9,
      "name": "Hydronium ion",
      "reference_temperature_K": 273.15,
      "reference_pressure_Pa": 101325,
      "ccs_method": "trajectory (MobCal-MPI 2.0)",
      "geometry_file": "molecules/H3O+.json"
    }
  }
}
```

2. reactions.json (default: `data/reactions_database_v1.json`)  
   Reaction database with ion–neutral reaction definitions and Arrhenius parameters.

Reaction database loaded via `ReactionConfig`
(real implementation: `ReactionLoader.cpp`)

Example:

```json
{
  "reactions": [
    {
      "id": "rxn_001_h3o_to_pentanalH+",
      "reactant": "H3O+",
      "product": "PentanalH+",
      "rate_constant": 3e-9,
      "order": [
        {
          "species": "Pentanal",
          "exponent": 1,
          "concentration_m3": 2.5e19
        }
      ]
    }
  ]
}
```

3. config.json  
   Main simulation configuration file defining instrument domains, fields, integrator settings, etc.

Minimal working IMS example (from `examples/ims/ims_basic.json`):

```json
{
  "simulation": {
    "total_time_s": 5e-4,
    "dt_s": 1e-9,
    "write_interval": 1000,
    "integrator": "RK4",
    "enable_gpu": false,
    "enable_openmp": true,
    "rng_seed": 123
  },
  "physics": {
    "collision_model": "HSS",
    "enable_space_charge": true
  },
  "output": {
    "folder": "results/ims",
    "trajectory_file": "ims_trajectories.h5",
    "print_progress": true
  },
  "ions": {
    "species": [
      {
        "id": "H3O+",
        "count": 100,
        "position": {
          "type": "gaussian",
          "center": [0.0, 0.0, 0.001],
          "std": [0.001, 0.001, 0.0005]
        },
        "velocity": {
          "type": "thermal",
          "temperature_K": 300.0
        }
      }
    ]
  },
  "domains": [
    {
      "name": "drift_region",
      "instrument": "IMS",
      "geometry": {
        "origin_m": [0.0, 0.0, -0.01],
        "length_m": 0.06,
        "radius_m": 0.5
      },
      "env": {
        "temperature_K": 300.0,
        "pressure_Pa": 200.0,
        "gas_species": "He",
        "gas_velocity_m_s": [0.0, 0.0, 0.0]
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

For more examples, see the `examples/` folder.

4. molecular_geometry.json (optional)  
   Molecular geometry definitions for ions (used in EHSS collision model).

Defined in each species database entry via "geometry_file" field.

Example (`data/molecules/H3O+.json`):

```json
{
  "molecule": {
    "name": "H3O+",
    "diameter_m": 3.3e-10,
    "CCS_m2": 24.9e-19,
    "atoms": [
      {
        "element": "O",
        "pos": [-0.00850062513799088, 0.013502526415391874, -0.03864593471061347],
        "mass_u": 15.999,
        "partial_charge_e": 0.00396,
        "LJ_sigma_angstrom": 2.4344,
        "LJ_epsilon_eV": 0.1034324933362249
      },
      {
        "element": "H",
        "pos": [0.9602283748620088, -0.08413547358460813, -0.0014669347106134722],
        "mass_u": 1.008,
        "partial_charge_e": 0.331991,
        "LJ_sigma_angstrom": 2.261,
        "LJ_epsilon_eV": 0.059579817756640736
      },
      {
        "element": "H",
        "pos": [-0.34485562513799106, 0.7037425264153919, 0.5612030652893866],
        "mass_u": 1.008,
        "partial_charge_e": 0.332032,
        "LJ_sigma_angstrom": 2.261,
        "LJ_epsilon_eV": 0.059579817756640736
      },
      {
        "element": "H",
        "pos": [-0.4804506251379914,-0.833919473584608, 0.05365306528938653],
        "mass_u": 1.008,
        "partial_charge_e": 0.332018,
        "LJ_sigma_angstrom": 2.261,
        "LJ_epsilon_eV": 0.059579817756640736
      }
    ]
  }
}
```

For more examples, see the `data/molecules/` folder.

# Output Format (HDF5)

Generated by `hdf5Writer.cpp`.

**Top-level structure**:

```css
/trajectory
  /time                 [T]
  /positions            [T × N × 3]
  /velocities           [T × N × 3]
  /species_ids          [T × N]
  /domain_indices       [T × N]
/ions (contains initial data)
/domains
  /domain_0
    /metadata
    /geometry
    /environment
    /fields
  /domain_1
    ...
/metadata
  /config (simulation parameters)
  /reproducibility (seed, version, build information)
  /system (basic system info)
  /species (metadata from species.json, only used species are stored)
  /reactions (metadata from reactions.json, only reactions involving used species are stored)
  /completion (performance metrics)

```

## Console Output

Clean, HPC-style runtime messages at simulation start:

```
╔════════════════════════════════════════════════════════════════════════════╗
║                                                                            ║
║             ██╗ ██████╗ █████╗ ██████╗ ██╗ ██████╗ ███╗   ██╗              ║
║             ██║██╔════╝██╔══██╗██╔══██╗██║██╔═══██╗████╗  ██║              ║
║             ██║██║     ███████║██████╔╝██║██║   ██║██╔██╗ ██║              ║
║             ██║██║     ██╔══██║██╔══██╗██║██║   ██║██║╚██╗██║              ║
║             ██║╚██████╗██║  ██║██║  ██║██║╚██████╔╝██║ ╚████║              ║
║             ╚═╝ ╚═════╝╚═╝  ╚═╝╚═╝  ╚═╝╚═╝ ╚═════╝ ╚═╝  ╚═══╝              ║
║                                                                            ║
║                  Ion Collision And Reaction IntegratiON                    ║
║                   High-Performance Trajectory Simulator                    ║
║                                                                            ║
╚════════════════════════════════════════════════════════════════════════════╝

   Version:      1.0.0
   Git Commit:   c5258a6 (Dec  1 2025 06:53:52)
   Build Type:   Release with OpenMP
   Compiler:     GCC 13.3.0

   License:      MIT
   Support:      https://github.com/ICARION-Project/ICARION/issues

────────────────────────────────────────────────────────────────────────────

 Configuration

   File:         examples/lqit/lqit_basic.json
   Log Level:    INFO

────────────────────────────────────────────────────────────────────────────

 System Information

   Hostname:     Workstation
   OS:           Linux 6.6.87.2-microsoft-standard-WSL2 (Ubuntu 24.04.2 LTS)
   CPU:          AMD Ryzen 9 7950X 16-Core Processor (32 cores)
   Memory:       30 GB RAM

   Threads:      16 OpenMP threads

════════════════════════════════════════════════════════════════════════════

 Starting simulation at 2025-12-01 09:02:42 UTC

════════════════════════════════════════════════════════════════════════════
```

**During execution:**
- Progress updates with ETA and performance metrics
- Final summary with walltime, ions/sec throughput
- Output file location and size

---

# Validation & Physics Benchmarks

ICARION v1.0 provides built-in validation configurations:

- IMS mobility vs Mason–Schamp
- Quadrupole a–q stability points
- Orbitrap axial frequency (±1%)
- TOF energy-time scaling
- LQIT secular frequency checks

---

# Directory Structure

```
ICARION/
├── build/                      # CMake build artifacts (generated)
├── cmake/                      # CMake modules and configuration
│   ├── CompilerSettings.cmake  # Compiler flags and optimizations
│   ├── CUDAConfig.cmake        # CUDA detection and configuration
│   ├── Dependencies.cmake      # External dependencies (HDF5, JsonCpp, OpenMP)
│   └── Targets.cmake           # Build targets and linking
├── data/                       # Reference databases
│   ├── molecules/              # Molecular geometry files (JSON)
│   ├── reactions_database_v1.json
│   └── species_database_v1.json
├── docs/                       # User documentation
│   ├── ARCHITECTURE.md         # System design and code structure
│   ├── CLI_USAGE.md            # Command-line interface reference
│   ├── COLLISION_MODELS.md     # Physics of collision models
│   ├── CONFIG_GUIDE.md         # Configuration file format
│   ├── DEVELOPERS_GUIDE.md     # Contributing and development setup
│   ├── HDF5_OUTPUT_STRUCTURE.md # Output file specification
│   ├── JSON_LOGGING.md         # Logging system documentation
│   └── TROUBLESHOOTING.md      # Common issues and solutions
├── examples/                   # Ready-to-run example configurations
│   ├── lqit/                   # Linear quadrupole ion trap examples
│   ├── orbitrap/               # Orbitrap mass analyzer examples
│   ├── quadrupole/             # Quadrupole mass filter examples
│   ├── reactions/              # Reaction kinetics examples
│   ├── ims/                    # Ion mobility spectrometry examples
│   ├── field_arrays/           # Pre-computed field data
│   ├── ion_clouds/             # Initial ion distributions
│   ├── waveforms/              # Custom RF/DC waveforms
│   └── README.md
├── results/                    # Simulation output (HDF5, JSON logs)
├── schema/                     # JSON schema definitions
│   ├── icarion-config.schema.json  # Main config schema
│   ├── species.schema.json
│   ├── reactions.schema.json
│   ├── validator.py            # Schema validation tool
│   └── README.md
├── scripts/                    # Utility scripts
│   ├── compute_ccs_maps.py     # CCS precomputation
│   ├── create_config.py        # Config file generator
│   └── README.md
├── src/                        # C++ source code
│   ├── core/                   # Core simulation engine
│   │   ├── collision/          # Collision model implementations
│   │   ├── domain/             # Multi-domain support
│   │   ├── geometry/           # Geometric boundaries
│   │   ├── integrator/         # ODE integrators (RK4, RK45, Boris)
│   │   ├── ion/                # Ion state management
│   │   └── simulation/         # Main simulation loop
│   ├── fieldsolver/            # Field computation module (future)
│   ├── main/                   # Entry point and CLI
│   ├── optimizer/              # Optimization module (future)
│   └── utils/                  # Utilities (logging, I/O, config)
├── tests/                      # Unit and integration tests
│   ├── unit/                   # Unit tests (Catch2)
│   ├── integration/            # End-to-end tests
│   ├── collision/              # Collision model tests
│   ├── instruments/            # Instrument-specific tests
│   └── README.md
├── tmp/                        # Temporary files and planning docs
│   └── docs/                   # Development planning documents
├── tools/                      # Standalone tools
│   └── ccs_precompute.cpp      # CCS precomputation utility
├── validation/                 # Physics validation suite
│   ├── configs/                # Validation configurations (309 configs)
│   │   ├── instruments/        # Instrument-specific tests
│   │   ├── physics/            # Physics benchmark tests
│   │   └── performance/        # Performance benchmarks
│   ├── scripts/                # Validation automation scripts
│   └── README.md
├── CMakeLists.txt              # Root CMake configuration
├── LICENSE                     # MIT License
└── README.md                   # This file
```

---

# Roadmap (v1.1 → v2.0)

- FieldSolver integration (BEM/FMM)
- Optimizer module (genetic algorithms, gradient descent)
- Advanced collision models (long-range potentials, velocity-dependent CCS)

---

# License

MIT, see LICENSE file.

---

# Citation

A CPC paper will follow.
Until then, please cite the GitHub repository.

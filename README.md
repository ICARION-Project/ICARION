ICARION v1.0 — Ion Collision And Reaction IntegratiON  
Modular C++/CUDA framework for multi-domain ion dynamics simulation.

---

<p align="center">
  <img src="https://img.shields.io/badge/C%2B%2B-17-blue?style=flat-square" />
  <img src="https://img.shields.io/badge/License-MIT-green?style=flat-square" />
</p>

---

# Release & API

- **Versioning:** v1.0.0 (semantic versioning). See `CHANGELOG.md`.
- **Stable surface:** JSON configuration schema is considered stable for v1.x.
- **Internal API:** C++ headers/classes are internal and may evolve between minor releases.
- **License:** MIT (see `LICENSE`); third-party dependencies listed in `CMakeLists.txt` and `cmake/`.
- **Experimental components (off-path for v1.0 results):** GPU EHSS geometry upload, GPU space-charge P³M, and adaptive field interpolation are present but incomplete and should be treated as experimental.

# What & Who

- **What is ICARION?** Modular ion trajectory simulator (C++17) for mass spectrometry, ion mobility devices and ion optics with collision and reaction support.
- **What can it do in v1.0?** IMS, RF quadrupole, Orbitrap, TOF, LQIT, FT-ICR; EHSS/HSS stochastic collision models Friction/Langevin/HardSphere deterministic collision models; Arrhenius reactions; RK4/RK45/Boris integrators; HDF5 with reproducibility metadata and config snapshot.
- **What can it not do yet?** No full-field solver, no optimizer loop, limited GPU coverage (see below).
- **Who is it for?** Researchers/engineers needing reproducible ion mobility / MS simulations with configurable physics and domains.
- **Expectation management:** ICARION prioritizes physical correctness and modularity. Performance optimization and GPU offloading are active development areas.
- **Integrator note:** RK45 now keeps per-ion adaptive state and is OpenMP-safe; batch paths (CPU/GPU) require uniform `dt` across active ions.
- **GPU status:** GPU support is provided as an **experimental backend**. While functional, it is not yet optimized for large-scale production simulations.
- **Output memory guard:** Set `output.buffer_byte_cap` (bytes) to cap in-memory trajectory buffering and fail fast before OOM; `0` disables the cap.

# Minimal Example (IMS)

```bash
cmake -B build -S . && cmake --build build
./build/src/icarion_main examples/ims/ims_basic.json
```

Key settings in `examples/ims/ims_basic.json`:
- `simulation.integrator: "RK4"`; `physics.collision_model: "HSS"`
- `simulation.rng_seed: 123` (reproducibility)
- `output.trajectory_file: "ims_trajectories.h5"` (writes HDF5 with version, git hash, build flags, RNG seed, integrator/collision model)
- `simulation.enable_gpu: false` (GPU is optional/experimental)

Run produces `results/ims/ims_trajectories.h5` (plus a config snapshot `ims_trajectories.config.json` alongside the output).

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

SSOT architecture:

- `config.json` (required): simulation, physics, output, domains, ions
- `species_database` (optional, defaults to `data/species_database_v1.json`)
- `reaction_database` (optional, defaults to `data/reactions_database_v1.json`)
- `molecular_geometry.json` (optional, referenced by species for EHSS)

## Output

- Open-source HDF5 format with hierarchical structure
- Chunked, extendable datasets
- Species metadata only once (reduces file size)
- Selected config/system metadata stored for reproducibility (full JSON not embedded in v1.0; full config snapshot now written alongside output as `<basename>.config.json`)

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

**GPU build:** add `-DUSE_GPU_ACCEL=ON` to enable CUDA paths (requires CUDA toolkit).

**Run**

```bash
./build/src/icarion_main input/config.json
```

The startup banner is intentionally minimal to avoid noisy logs in batch/CI runs.

---

# Input Files

ICARION uses an SSOT input architecture:

1. `config.json` — simulation, physics, output, domains, ions (required)
2. `species_database` — optional, defaults to `data/species_database_v1.json` (species properties)
3. `reaction_database` — optional, defaults to `data/reactions_database_v1.json` (reaction kinetics)
4. `molecular_geometry.json` — optional, referenced from species for EHSS (geometries)

Species database loaded via `SpeciesConfig` (`SpeciesLoader.cpp`).

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
      "environment": {
        "temperature_K": 300.0,
        "pressure_Pa": 200.0,
        "gas_species": "He",
        "gas_velocity_m_s": [0.0, 0.0, 0.0]
      },
      "fields": {
        "dc": {
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
  /config (selected simulation parameters)
  /reproducibility (seed, version, build information)
  /system (basic system info)
  /species (metadata from species.json, only used species are stored)
  /reactions (metadata from reactions.json, only reactions involving used species are stored)
  /completion (performance metrics)

```

## Console Output

Startup logging shows a short header with version/build info followed by progress
updates (ETA) and a final summary (walltime, ions/sec, output path). No ASCII banner
is printed by default.

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
│   │   ├── config/             # FullConfig types, loader, validation
│   │   ├── gpu/                # CUDA helpers and kernels
│   │   ├── integrator/         # SimulationEngine, strategies, domains
│   │   ├── io/                 # HDF5 writer/reader
│   │   ├── log/                # Logging wrappers
│   │   ├── param/              # Legacy/bridge parameters
│   │   ├── physics/            # Forces, collisions, reactions
│   │   ├── types/              # IonState, IonEnsemble, Vec3, enums
│   │   └── utils/              # Math, safety, profiling helpers
│   ├── fieldsolver/            # Field computation module (future)
│   ├── main/                   # Entry point and setup wiring
│   ├── optimizer/              # Optimization module (future)
│   └── utils/                  # Shared utilities (logging, math, CLI helpers)
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

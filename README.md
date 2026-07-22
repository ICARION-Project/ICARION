ICARION v1.1.1 — Ion Collision And Reaction IntegratiON
Modular C++/CUDA framework for multi-domain ion dynamics simulation.

---

<p align="center">
  <img src="https://img.shields.io/badge/C%2B%2B-17-blue?style=flat-square" />
  <img src="https://img.shields.io/badge/License-GPL--3.0--only-blue?style=flat-square" />
  <a href="https://icarion.readthedocs.io/en/latest/?badge=latest">
    <img src="https://readthedocs.org/projects/icarion/badge/?version=latest" alt="Documentation Status" />
  </a>
  <a href="https://github.com/ICARION-Project/ICARION/releases/latest">
    <img src="https://img.shields.io/github/v/release/ICARION-Project/ICARION?style=flat-square" alt="Latest Release" />
  </a>
</p>

---

# Download

Latest release packages are available here:

[Download latest release](https://github.com/ICARION-Project/ICARION/releases/latest)

Start with `ICARION-Launcher-Guide.md` if you use the Windows or Linux launcher.

# Release & API

- **Versioning:** v1.1.1 (semantic versioning). See `CHANGELOG.md`.
- **Stable surface:** JSON configuration schema is considered stable for v1.1.x.
- **Internal API:** C++ headers/classes are internal and may evolve between minor releases.
- **License:** GPL-3.0-only (see `LICENSE`); third-party dependencies listed in `CMakeLists.txt` and `cmake/`.
- **Experimental components (off-path for v1.1.0 results):** GPU EHSS geometry upload, GPU space-charge P³M, and adaptive field interpolation are present but incomplete; the primary runtime GPU path in `SimulationEngine` is disabled for v1.1.0 (helpers remain buildable for dev/testing).

# Documentation

The user documentation is available at:

<https://icarion.readthedocs.io/en/latest/>

It includes installation instructions, first simulations, configuration files,
species and reaction databases, collision models, HDF5 output, analysis
workflows, and validation.

ICARION builds on established trajectory method and ion trajectory simulation traditions; see [Related and Complementary Software](docs/RELATED_SOFTWARE.md). Its [self-describing IPM tables](https://icarion.readthedocs.io/en/latest/ipm-precomputation/) retain offline collision provenance for runtime use.

Additional repository documentation:

- [docs/CONFIG_GUIDE.md](docs/CONFIG_GUIDE.md) — configuration fields, schema, and validation.
- [docs/CLI_USAGE.md](docs/CLI_USAGE.md) — command-line flags and batch usage.
- [docs/LAUNCHER_GUIDE.md](docs/LAUNCHER_GUIDE.md) — minimal Windows/Linux launcher usage and basic analysis.
- [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) — common build/run issues (incl. WSL).
- [docs/HDF5_OUTPUT_STRUCTURE.md](docs/HDF5_OUTPUT_STRUCTURE.md) — output file layout.
- [docs/COLLISION_MODELS.md](docs/COLLISION_MODELS.md) — collision physics background and model guidance.
- [docs/RELATED_SOFTWARE.md](docs/RELATED_SOFTWARE.md) — complementary CCS, mobility, and instrument-simulation software.
- [docs/ARCHITECTURE.md](docs/ARCHITECTURE.md) — high-level design and module overview.
- [docs/GPU_ARCHITECTURE.md](docs/GPU_ARCHITECTURE.md) — GPU code structure and data flow.
- [docs/JSON_LOGGING.md](docs/JSON_LOGGING.md) — structured logging format and examples.
- [validation/README.md](validation/README.md) — validation suite and ready-to-run configs.
- [validation/VALIDATION_REPORT_v1.0.0.md](validation/VALIDATION_REPORT_v1.0.0.md) — detailed physics validation baseline.

# What & Who

- **What is ICARION?** Modular ion trajectory simulator (C++17) for mass spectrometry, ion mobility devices and ion optics with collision and reaction support.
- **What can it do in v1.1.0?** IMS, TIMS, RF quadrupole, Orbitrap, TOF, LQIT, FT-ICR; EHSS/HSS/IPM stochastic collision models; Friction/Langevin/HardSphere deterministic collision models; Arrhenius and equilibrium-linked reactions; RK4/RK45/Boris integrators; HDF5 with reproducibility metadata, config snapshot, compact output, and collision diagnostics.
- **What can it not do yet?** No full-field solver, no optimizer loop, limited GPU coverage (see below), magnetic field map providers not wired (analytical/uniform B only).
- **Who is it for?** Researchers/engineers needing reproducible ion mobility / MS simulations with configurable physics and domains or scientists researching on ion transport regimes/phenomena.
- **Expectation management:** ICARION prioritizes physical correctness and modularity. Performance optimization and GPU offloading are active development areas.
- **Integrator note:** RK45 keeps per-ion adaptive state; OpenMP determinism is covered by tests. Batch paths (CPU/GPU) require uniform `dt` across active ions.
- **GPU status:** GPU codepaths are compiled but the primary runtime GPU path is disabled for v1.1.0; helpers remain experimental for developers.
- **Output memory guard:** Set `output.buffer_byte_cap` (bytes) to cap in-memory trajectory buffering and fail fast before OOM; `0` disables the cap.

# Keywords & Acronyms

- IMS: Ion Mobility Spectrometry
- LQIT: Linear Quadrupole Ion Trap
- TOF: Time-of-Flight analyzer
- FT-ICR: Fourier Transform Ion Cyclotron Resonance
- CCS: Collision Cross Section
- EHSS / HSS: (Exact) Hard Sphere Scattering collision models
- P³M: Particle–Particle Particle–Mesh space-charge solver
- EN_Td / E/N: Reduced electric field (Townsend; 1 Td = 1e-21 V·m²)
- RK4 / RK45: Fixed-step and adaptive Runge–Kutta integrators (fourth order or fourth order with fifth-order error control)
- SoA / AoS: Structure of Arrays / Array of Structures data layout

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
- Domain-dependent gas environments
- DC, RF, AC, magnetic fields (magnetic maps are not wired; analytical/uniform B only; electric fields either analytical or initialized from precomputed field arrays)

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
- GPU support (compiled; primary runtime GPU path disabled in v1.1.0; helpers remain experimental)

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
- Config, species/reaction DBs, and field arrays are embedded in HDF5 (larger files but self-contained) plus a `<basename>.config.json` snapshot; structured metadata groups remain for quick inspection.

---

# Building ICARION

## Requirements (CPU build)

- C++17 compiler (GCC/Clang/MSVC)
- CMake ≥ 3.16
- HDF5 ≥ 1.10, Eigen3 ≥ 3.4, jsoncpp ≥ 1.9, nlohmann_json ≥ 3.10, OpenSSL, BLAS, spdlog, cxxopts
- Catch2 v3 (for tests/ctest)
- OpenMP recommended (auto-detected if available)
- (Optional) CUDA ≥ 11 **only** when `-DUSE_GPU_ACCEL=ON` (default is OFF; CPU works without CUDA).

### Quick start (Ubuntu/WSL 22.04 / 24.04, CPU)

```bash
sudo apt update && sudo apt install -y \
  build-essential cmake git pkg-config \
  libeigen3-dev libjsoncpp-dev nlohmann-json3-dev libhdf5-dev \
  libssl-dev libopenblas-dev libspdlog-dev libcxxopts-dev

# Optional (tests/ctest): Catch2 v3. If your distro only ships Catch2 v2,
# install Catch2 v3 from source.

cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
./build/src/icarion_main examples/ims/ims_basic.json
```

### Build an installable package

```bash
cmake -B build-package -S . \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX=/usr \
  -DUSE_GPU_ACCEL=OFF \
  -DICARION_BUILD_TESTS=OFF
cmake --build build-package -j"$(nproc)"
cmake --build build-package --target package
```

This creates a `.deb` package on Linux plus a `.tar.gz` install archive. See
[docs/PACKAGING.md](docs/PACKAGING.md) for install paths and package usage.
See [docs/LAUNCHER_GUIDE.md](docs/LAUNCHER_GUIDE.md) for the minimal
Windows/Linux launchers and basic analysis buttons.
Pushing a tag like `v1.1.0` also triggers the GitHub Actions release
packaging workflow, including a Windows `.zip` artifact.

### GPU build (optional, dev/experimental)

```bash
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release -DUSE_GPU_ACCEL=ON
cmake --build build -j"$(nproc)"
```

CUDA toolkit and drivers must be installed; runtime falls back to CPU for v1.1.0.

### Tested on (current developer setup)

- Ubuntu 24.04 LTS on WSL2 — CPU build; CUDA 12.0 toolkit build tested with `-DUSE_GPU_ACCEL=ON` (runtime still CPU in v1.1.0)
- Notes: broader platform coverage to follow as we add validation runs.

**From a fresh clone:**

```bash
git clone https://github.com/ICARION-Project/ICARION.git
cd ICARION
cmake -B build -S . -DCMAKE_BUILD_TYPE=Release
cmake --build build -j"$(nproc)"
```

**Troubleshooting:** see [docs/TROUBLESHOOTING.md](docs/TROUBLESHOOTING.md) for missing packages (HDF5, BLAS/OpenMP, CUDA detection on WSL) and cache cleanup tips.

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
      "rate_constant": 3e-15,
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
    /name
    /instrument
    /solver
    /geometry
    /environment
    /fields
      /dc
      /rf
      /ac
      /waveforms (optional)
  /domain_1
    ...
/metadata
  /config (selected simulation parameters)
  /physics (collision/reaction handler metadata)
  /reproducibility (seed, version, build information)
  /system (basic system info)
  /species (metadata from species.json, only used species are stored)
  /reactions (metadata from reactions.json, only reactions involving used species are stored)
  /completion (run status + final time/active ions)

```

## Console Output

Startup logging shows a short header with version/build info followed by progress
updates (ETA) and a final summary (walltime, ions/sec, output path). No ASCII banner
is printed by default.

---

# Validation & Physics Benchmarks

ICARION v1.1.0 provides built-in validation configurations:

- IMS mobility vs Mason–Schamp
- Quadrupole a–q stability points
- Orbitrap axial frequency (±1%)
- TOF energy-time scaling
- LQIT secular frequency checks
- dynamic equilibrium reactions
- IPM offline sample generation and runtime consumption
- TIMS elution ordering
- multi-event collision sampling

Full physics validation baseline results are documented in `validation/VALIDATION_REPORT_v1.0.0.md` (see `validation/README.md` for how to run the suite).

Fast regression tests (CTests) live under `tests/` and are recommended after a build.
See `tests/README.md` for the complete list and tags. Run from `build/`:

```bash
ctest --output-on-failure
```

---

# Directory Structure

```
ICARION/
├── CHANGELOG.md                # Release history
├── CMakeLists.txt              # Root CMake configuration
├── LICENSE                     # GNU GPL v3.0
├── README.md                   # This file
├── build/                      # CMake build artifacts (generated)
├── analysis/                   # Analysis scripts/outputs
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
│   ├── releases/               # Release notes for x.y.0 releases
│   ├── HDF5_OUTPUT_STRUCTURE.md # Output file specification
│   ├── JSON_LOGGING.md         # Logging system documentation
│   └── TROUBLESHOOTING.md      # Common issues and solutions
├── examples/                   # Ready-to-run example configurations
│   ├── fticr/                  # FT-ICR examples
│   ├── lqit/                   # Linear quadrupole ion trap examples
│   ├── orbitrap/               # Orbitrap mass analyzer examples
│   ├── quadrupole/             # Quadrupole mass filter examples
│   ├── reactions/              # Reaction kinetics examples
│   ├── ims/                    # Ion mobility spectrometry examples
│   ├── field_arrays/           # Pre-computed field data
│   ├── ion_clouds/             # Initial ion distributions
│   ├── tof/                    # Time-of-flight examples
│   ├── waveforms/              # Custom RF/DC waveforms
│   └── README.md
├── phase0_analysis/            # Legacy analysis runs
├── schema/                     # JSON schema definitions
│   ├── icarion-config.schema.json  # Main config schema
│   ├── species.schema.json
│   ├── reactions.schema.json
│   ├── validator.py            # Schema validation tool
│   └── README.md
├── scripts/                    # Utility scripts
│   ├── create_config.py        # Config file generator
│   ├── generate_ims_en_sweep.py # IMS E/N sweep generator
│   ├── generate_lqit_ac_sweep.py # LQIT AC sweep generator
│   ├── run_ims_en_sweep.sh      # IMS sweep runner
│   ├── run_lqit_ac_sweep.sh     # LQIT sweep runner
│   └── README.md
├── src/                        # C++ source code
│   ├── core/                   # Core simulation engine
│   │   ├── config/             # FullConfig types, loader, validation
│   │   ├── gpu/                # CUDA helpers and kernels
│   │   ├── integrator/         # SimulationEngine, strategies, domains
│   │   ├── io/                 # HDF5 writer/reader
│   │   ├── log/                # Logging wrappers
│   │   ├── physics/            # Forces, collisions, reactions
│   │   ├── types/              # IonState, IonEnsemble, Vec3, enums
│   │   └── utils/              # Math, safety, profiling helpers
│   ├── fieldsolver/            # Field computation module (future)
│   ├── main/                   # Entry point and setup wiring
│   ├── optimizer/              # Optimization module (future)
│   └── utils/                  # Shared utilities (logging, math, CLI helpers)
├── tests/                      # Unit and integration tests
│   ├── config/                 # Config loader tests
│   ├── gpu/                    # GPU kernel tests
│   ├── helpers/                # Test helpers
│   ├── instruments/            # Instrument-specific tests
│   ├── integration/            # End-to-end tests
│   ├── integrator/             # Integrator/engine tests
│   ├── io/                     # HDF5 writer tests
│   ├── physics/                # Forces/collisions/reactions tests
│   ├── unit/                   # Unit tests (Catch2)
│   ├── utils/                  # Hashing and utility tests
│   └── README.md
├── tools/                      # Standalone tools
│   ├── ccs_precompute.cpp
│   ├── ehss_samples_precompute.cpp
│   ├── ehss_offline_precompute.cpp
│   └── interaction_potential_precompute.cpp
├── validation/                 # Physics validation suite
│   ├── configs/                # Validation configurations (309 configs)
│   │   ├── instruments/        # Instrument-specific tests
│   │   ├── physics/            # Physics benchmark tests
│   │   └── performance/        # Performance benchmarks
│   ├── scripts/                # Validation automation scripts
│   ├── figures/                # Validation plots
│   ├── logs/                   # Validation logs
│   ├── results/                # Validation outputs
│   ├── README.md
│   └── VALIDATION_REPORT_v1.0.0.md
```

---

# Roadmap (Future Releases)

- FieldSolver integration (BEM/FMM)
- Optimizer module (genetic algorithms, gradient descent)
- Advanced collision or reaction models (e.g., long-range potentials)

---

# License

GNU GPL v3.0, see LICENSE file.

---

# Citation

A peer-reviewed journal article describing ICARION is planned. Until the article
is available, please cite the corresponding ICARION
[software release](https://doi.org/10.5281/zenodo.20599037) and
[GitHub repository](https://github.com/ICARION-Project/ICARION).

Once the journal article has been published, please cite the article in addition
to the software release when using ICARION in scientific work.

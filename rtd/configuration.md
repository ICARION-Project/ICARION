# Configuration files

ICARION simulations are controlled by JSON configuration files. A good configuration should answer six questions:

1. **How long should the simulation run?** (`simulation`)
2. **Which physics should be active?** (`physics`)
3. **Where should output be written?** (`output`)
4. **Which instrument regions exist?** (`domains`)
5. **Which ions are initialized?** (`ions`)
6. **Which species and reactions should be known?** (`species_database`, `reaction_database`)

The main repository also contains the detailed reference file `docs/CONFIG_GUIDE.md`; this RTD page is intended as a user-facing starting point.

---

## Minimal mental model

A configuration is a single JSON object with the following common sections:

```json
{
  "title": "short human-readable simulation name",
  "species_database": "data/species_database_v1.json",
  "reaction_database": "data/reactions_database_v1.json",
  "simulation": {},
  "physics": {},
  "output": {},
  "domains": [],
  "ions": {}
}
```

The schema requires `simulation`, `physics`, `output`, and at least one entry in `domains`. A normal trajectory simulation also needs `ions`, even though specialized workflows can load or generate ion ensembles differently. For example, `reaction_database` is only needed if reactions are enabled. In practice, new users should start from an example file and modify one parameter at a time.

---

## Required core sections

### `simulation`

This section defines the numerical time axis and integrator.

| Parameter | Typical value | Unit | Meaning |
|---|---:|---:|---|
| `total_time_s` | `1e-4` to `1` | s | Physical simulation time. Must be long enough for ions to cross the domain or reach the relevant trapping time. |
| `dt_s` | `1e-10` to `1e-7` | s | Base integration step. Choose smaller values for fast RF fields, high collision rates, and high kinetic energy. |
| `integrator` | `"RK4"`, `"RK45"`, `"Boris"` | - | Numerical integrator. Use `Boris` mainly for magnetic-field problems; `RK4`/`RK45` for most electrostatic/RF examples. |
| `write_interval` | `10` to `1000` | steps | Store every n-th integration step. Increase this to reduce HDF5 size. |
| `rng_seed` | integer | - | Global random seed for reproducible stochastic collisions and reactions. |
| `enable_openmp` | `true` | - | Enable OpenMP parallel execution. |
| `enable_gpu` | `false` | - | Request GPU execution where supported. CPU is the validated default for v1.1. |

Example:

```json
{
  "simulation": {
    "total_time_s": 0.001,
    "dt_s": 1e-9,
    "integrator": "RK4",
    "write_interval": 100,
    "rng_seed": 42,
    "enable_openmp":true,
    "enable_gpu": false
  }
}
```

!!! tip
    Start with a small ion count and a short `total_time_s`. Once the configuration loads and produces reasonable output, increase the ensemble size and simulation time.

---

### `physics`

This section selects the active collision, reaction, and space charge models. For more information on the used collision models, see [Collision models](collision-models.md).

!!! note
    `physics.collision_model` is a global setting for the whole run in v1.1. Different domains can use different geometry, fields, gas pressure, temperature, and gas species, but they do not select separate collision models.

| Parameter | Example | Meaning |
|---|---|---|
| `collision_model` | `"NoCollisions"`, `"Friction"`, `"HSS"`, `"EHSS"`, `"InteractionPotentialModel"` | Selects ion-neutral collision handling. |
| `collision_subcycles_per_step` | `1` | Splits each stochastic collision update into smaller substeps; useful when event probabilities become large, but still an approximation. For the most event-resolved calculations, choose `dt_s` so that at most one collision per ion and base time step is expected. |
| `collision_multi_event_mode` | `false` | Enables a practical multi-event collision approximation by enforcing enough subcycles to allow multiple stochastic collision opportunities within one macro step. |
| `collision_max_events_per_step` | `1` | Minimum subcycle count used by `collision_multi_event_mode`. |
| `enable_reactions` | `true` / `false` | Enables stochastic species conversion from the reaction database based on configured reaction kinetics. |
| `enable_space_charge` | `true` / `false` | Enables ion-ion space charge effects if configured. |
| `space_charge_model` | `"auto"`, `"direct"`, `"grid"`, `"gpu"` | Selects the space charge backend. `auto` chooses an appropriate CPU model unless GPU space charge is explicitly requested and available. |
| `ipm_orientation_mode` | `"random"` or `"fixed"` | Orientation sampling for `InteractionPotentialModel`; omitted configs default to `"random"`. |
| `ipm_fixed_orientation_index` | `0` | Orientation index used when `ipm_orientation_mode` is `"fixed"`. |
| `ipm_vrel_log_prefix` | `"debug/ipm_vrel"` | Optional relative-velocity diagnostic CSV prefix for `InteractionPotentialModel`. |
| `ipm_momentum_log_prefix` | `"debug/ipm_dp"` | Optional momentum-kick diagnostic CSV prefix for `InteractionPotentialModel`. |

Example without reactions:

```json
{
  "physics": {
    "collision_model": "HSS",
    "enable_reactions": false,
    "enable_space_charge": false
  }
}
```

Example with reactions:

```json
{
  "species_database": "data/species_database_v1.json",
  "reaction_database": "data/reactions_database_v1.json",
  "physics": {
    "collision_model": "HSS",
    "enable_reactions": true
  }
}
```

---

### `output`

This section controls where ICARION writes HDF5 trajectories and snapshots.

| Parameter | Example | Meaning |
|---|---|---|
| `folder` | `"results/ims"` | Output directory. |
| `trajectory_file` | `"ims_trajectories.h5"` | Main HDF5 output file name. |
| `trajectory_mode` | `"full"` or `"minimal"` | Full trajectories or compact final-state only output. |

Example:

```json
{
  "output": {
    "folder": "results/ims",
    "trajectory_file": "ims_trajectories.h5",
    "trajectory_mode": "full"
  }
}
```

See [Output files](output-files.md) and [Analysis workflow](analysis.md) for what to do with the generated HDF5 file.

---

## Domains

`domains` is an array of instrument regions. A simple IMS simulation may contain one domain; a multi-dommanin (e.g., IMS-quadrupole) simulation contains at least two. See [Multi-domain simulations](multi-domain.md) for more details.

Each domain should define:

| Group | Typical fields | Meaning |
|---|---|---|
| identity | `name`, `instrument` | Human-readable name and instrument/domain type. |
| geometry | `origin_m`, `length_m`, `radius_m` | Spatial extent and geometry boundaries. When ions reach these boundaries, the boundary actions specified by `boundary` will be enforced. |
| environment | `pressure_Pa`, `temperature_K`, `gas_species` | Neutral gas conditions used by collisions, reactions, and transport. |
| fields | instrument-dependent | DC, RF, AC, and/or magnetic fields for analytic fields; or imported field arrays. |
| boundary | `type` | Specifies what happens when ions reach a domain wall. These including absorption or different reflection descriptions, see [Boundary conditions](boundary-conditions.md). |

Available `instrument` types include: 
* ion mobility spectrometers: `IMS`
* trapped ion mobility spectrometers: `TIMS`
* linear (quadrupole) ion traps: `LQIT`
* Orbitraps: `ORBITRAP` 
* RF-Quadrupole (either as ion guide or mass filter, depending on specified voltages): `QUADRUPOLE`
* time-of-flight regions: `TOF`
* FT-ICR: `FT-ICR`, or `FTICR`

Minimal IMS domain skeleton:

```json
{
  "domains": [
    {
      "name": "drift_region",
      "instrument": "IMS",
      "geometry": {
        "origin_m": [0.0, 0.0, 0.0],
        "length_m": 0.1,
        "radius_m": 0.01
      },
      "environment": {
        "pressure_Pa": 2000.0,
        "temperature_K": 300.0,
        "gas_species": "He"
      }
    }
  ]
}
```

### Choosing domain parameters

For a first IMS run:

- Use one domain with `instrument = "IMS"`.
- Use He if you want to compare to the provided species database and validation examples.
- Start at moderate pressure, e.g. `200 Pa`, for stochastic HSS/EHSS runs.
- Use a with appropriate length that the ion cloud leaves the starting region during `total_time_s`.

For RF devices:

- Make sure `dt_s` resolves the RF period. As a rule of thumb, use many integration steps (at least more than 10, ideally more) per RF cycle.
- Keep the first run short and inspect whether ions remain inside the geometry.

For coupled simulations with multiple instrument domains, see [Multi-domain simulations](multi-domain.md).

Field keys, waveform references, field array implementation, and magnetic field options are
covered in [Fields and waveforms](fields-waveforms.md), with all required voltage parameters for each instrument summarized. 
Boundary behavior is covered in [Boundary conditions](boundary-conditions.md).

---

## Ion initialization

The `ions` section defines which ion species are simulated and how their initial positions and velocities are sampled.

The common pattern is:

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
          "temperature_K": 300.0
        }
      }
    ]
  }
}
```

Available position distributions include:

| Type | Required parameters | Meaning |
|---|---|---|
| `point` | `center` | All ions start at one coordinate. |
| `gaussian` | `center`, `std` | Normal distribution with independent x/y/z widths. |
| `uniform_sphere` | `center`, `radius` | Uniform sphere. |
| `uniform_cylinder` | `center`, `radius`, `length` | Uniform cylinder aligned with z. |
| `uniform_box` | `min`, `max` | Uniform rectangular box. |

Available velocity distributions include:

| Type | Required parameters | Meaning |
|---|---|---|
| `thermal` | `temperature_K` | Maxwellian thermal velocities with uniform random orientation sampling. |
| `fixed` | `value` | Same velocity vector for all ions. |
| `kinetic` | `energy_eV`, `direction` | Directed kinetic energy, optionally with `spread_angle_deg`. |
| `gaussian` | `mean`, `std` | Gaussian velocity components. |

!!! warning
    The `id` in `ions.species[].id` must match a species ID in the species database. If the species database is missing or the ID differs by spelling/case, the run may fail validation.

---

## Species and reaction database paths

Use explicit database paths when possible:

```json
{
  "species_database": "data/species_database_v1.json",
  "reaction_database": "data/reactions_database_v1.json"
}
```

If a path is not specified, ICARION can search for fallback database files under `data/`, but explicit paths make published simulations easier to reproduce.

See [Species database](species-database.md) and [Reactions](reactions.md) for examples.

---

## How to validate a config before a long run

Use the CLI dry-run mode first:

```bash
icarion --dry-run config.json
```

or, if the executable is still in the build tree:

```bash
./build/src/icarion_main --dry-run config.json
```

Then run with verbose logging:

```bash
icarion --verbose config.json
```

Useful checks before running a large ensemble:

- Does every ion species exist in the species database?
- Are all units SI in the config file unless explicitly documented otherwise?
- Is `dt_s` small enough for RF fields/collisions?
- Is `write_interval` large enough to avoid huge HDF5 files but small enough to resolve expected behavior?
- Is the output folder unique for this run?
- Are reactions disabled unless a reaction database is actually intended?

After preparing a configuration, run a [first simulation](first-simulation.md).

---

## A practical first workflow

1. Copy an existing example JSON.
2. Change only `output.folder` and `output.trajectory_file`.
3. Run `--dry-run`.
4. Run with 100 to 1000 ions.
5. Inspect the HDF5 file and logs.
6. Only then increase ion count, simulation time, and output density.

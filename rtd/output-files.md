# Output files

ICARION writes simulation results as HDF5 files. The output is inspectable without rerunning the simulation: it stores trajectory arrays, initial ion metadata, domain settings, selected configuration values, reproducibility metadata, and optional compact analysis output.

For the first example run:

```bash
./build/src/icarion_main examples/ims/ims_basic.json
```

the default files are:

```text
results/ims/
├── ims_trajectories.h5
└── ims_trajectories.config.json
```

The `.config.json` file is the resolved configuration snapshot written next to the HDF5 file. The same snapshot is also embedded in the HDF5 metadata when available.

---

## Configure output

The output path is controlled by the `output` section:

```json
{
  "output": {
    "folder": "results/ims",
    "trajectory_file": "ims_trajectories.h5",
    "trajectory_mode": "full",
    "print_progress": true
  }
}
```

| Field | Meaning |
|---|---|
| `folder` | Directory where ICARION writes the output file and config snapshot. |
| `trajectory_file` | HDF5 file name inside `folder`. |
| `trajectory_mode` | `"full"` for time-series trajectories, `"minimal"` for compact final-state output. |
| `print_progress` | Print runtime progress to the console. |
| `buffer_byte_cap` | Optional trajectory buffer memory cap in bytes; `0` means unlimited. |

You can override output paths from the CLI (see also the [CLI Reference](cli-reference.md)):

```bash
./build/src/icarion_main \
  --output-dir results/test_run \
  --output test_run.h5 \
  examples/ims/ims_basic.json
```

---

## High-level HDF5 layout

A normal full-output file contains:

```text
ims_trajectories.h5
├── metadata/
├── trajectory/
├── ions/
├── domains/
└── analysis/        # optional; created by selected output modes or diagnostics
```

Not every subgroup is present in every run. For example, `/metadata/reactions/` is only written when reactions are enabled and matching reactions are loaded.

---

## `/metadata/`

The metadata group stores information needed to understand and reproduce the run.

| Group | Purpose |
|---|---|
| `/metadata/config/` | Key simulation/config values and resolved config JSON. |
| `/metadata/config/integrator_params/` | Integrator-related runtime metadata. |
| `/metadata/config/derived_summary/` | Counts such as number of domains, species, and reactions. |
| `/metadata/physics/` | Selected collision and reaction handler metadata. |
| `/metadata/reproducibility/` | RNG seed, git hash, build information, input hashes, embedded inputs. |
| `/metadata/system/` | Host, OS, CPU, memory, timestamp, and GPU info if available. |
| `/metadata/species/` | Species table for species used by the ion ensemble. |
| `/metadata/reactions/` | Reaction table for reactions matching used species. |
| `/metadata/completion/` | Completion status, final time, active ion count, and timestamp. |

Important datasets:

| Dataset | Meaning |
|---|---|
| `/metadata/config/config_json` | Resolved config snapshot as JSON text. |
| `/metadata/config/format_version` | HDF5 metadata format version. |
| `/metadata/config/dt_s` | Base timestep from the config. |
| `/metadata/config/total_time_s` | Configured physical simulation time. |
| `/metadata/config/write_interval` | Output interval in integration steps. |
| `/metadata/config/collision_model` | Global collision model used by the run. |
| `/metadata/config/trajectory_mode` | Output mode: `full` or `minimal`. |
| `/metadata/reproducibility/global_seed` | RNG seed used for stochastic processes. |
| `/metadata/reproducibility/git_hash` | Git commit hash captured at build/runtime. |
| `/metadata/reproducibility/git_dirty` | Whether the source tree was dirty at build/runtime. |
| `/metadata/reproducibility/input_hash/config_sha256` | Hash of the input config file if available. |
| `/metadata/reproducibility/input_hash/species_db_sha256` | Hash of the species database if available. |
| `/metadata/reproducibility/input_hash/reaction_db_sha256` | Hash of the reaction database if available. |

For publication or debugging, start with `/metadata/config/config_json`, `/metadata/reproducibility/`, and `/metadata/completion/`.

---

## `/trajectory/`

The trajectory group stores time-dependent ion arrays when `output.trajectory_mode = "full"`.

| Dataset | Shape | Meaning |
|---|---:|---|
| `time` | `[T]` | Stored output times. |
| `positions` | `[T, N, 3]` | Ion positions in m. |
| `velocities` | `[T, N, 3]` | Ion velocities in m/s. |
| `domain_indices` | `[T, N]` | Active domain index for each ion at each stored frame. |
| `species_id_indices` | `[T, N]` | Species-pool index for each ion at each stored frame. |
| `species_ids` | `[T, N]` | Compatibility string labels when written. |
| `time_per_ion` | `[T, N]` | Per-ion time when adaptive/per-ion stepping is active. |

`T` is the number of stored output frames, not the number of integration steps. It depends on:

```text
T ≈ total_steps / write_interval + final frame
```

Large files usually come from high ion counts, long simulation times, or small `write_interval`.

!!! tip
    If the HDF5 file is too large, increase `simulation.write_interval`, reduce ion count, shorten `simulation.total_time_s`, or use `trajectory_mode = "minimal"` when only final-state data are needed.

---

## `/ions/`

The ions group stores per-ion metadata captured at initialization and updated at finalization.

| Dataset | Shape | Meaning |
|---|---:|---|
| `initial_species_id` | `[N]` | Initial species label. |
| `initial_pos_x`, `initial_pos_y`, `initial_pos_z` | `[N]` | Initial position components in m. |
| `initial_vel_x`, `initial_vel_y`, `initial_vel_z` | `[N]` | Initial velocity components in m/s. |
| `birth_time_s` | `[N]` | Ion birth time. |
| `death_time_s` | `[N]` | Ion loss/death time; `-1` for ions still active at the end. |
| `charge_C` | `[N]` | Ion charge in Coulomb. |

Use this group to check whether the generated ion cloud matches your intended species, positions, velocities, and birth times.

---

## `/domains/`

The domains group stores a copy of the domain settings used by the run.

```text
domains/
├── domain_0/
│   ├── name
│   ├── instrument
│   ├── solver
│   ├── geometry/
│   ├── environment/
│   └── fields/
└── domain_1/
    └── ...
```

Common domain datasets:

| Dataset | Meaning |
|---|---|
| `domain_<i>/name` | Domain name from the config. |
| `domain_<i>/instrument` | Instrument/domain type. |
| `domain_<i>/solver` | Domain solver/integrator selection. |
| `domain_<i>/geometry/origin_m` | Domain origin in global coordinates. |
| `domain_<i>/geometry/length_m` | Domain axial length. |
| `domain_<i>/geometry/radius_m` | Cylindrical radius where applicable. |
| `domain_<i>/environment/pressure_Pa` | Gas pressure. |
| `domain_<i>/environment/temperature_K` | Gas temperature. |
| `domain_<i>/environment/gas_species` | Gas species label. |
| `domain_<i>/environment/gas_velocity_ms` | Bulk gas flow velocity. |
| `domain_<i>/fields/dc/`, `rf/`, `ac/` | Static or `t=0` field settings. |

This group is useful when comparing two HDF5 files and checking whether the intended geometry, gas, and field settings were actually used.

---

## `/metadata/species/`

The species metadata table is filtered to species used by the ion ensemble.

| Dataset | Meaning |
|---|---|
| `names` | Species IDs. |
| `mass_kg` | Species mass in kg. |
| `charge_C` | Species charge in Coulomb. |
| `mobility_m2Vs` | Mobility converted to SI units. |
| `ccs_m2` | Generic collision cross section converted to m². |

For time-dependent species analysis, combine `/metadata/species/names` or `/analysis/minimal_transport/species_pool` with `trajectory/species_id_indices`.

---

## `/metadata/reactions/`

If reactions are enabled and the reaction database matches the used species, the reactions metadata group records the active reaction table.

| Dataset | Meaning |
|---|---|
| `id` | Reaction ID. |
| `reactant_1` | Reactant species ID. |
| `reactant_2` | Reserved/secondary reactant field. |
| `product_1` | Product species ID. |
| `rate_constant_m3s` | Stored rate constant value. |
| `type` | Internal reaction type code. |

`rate_constant_m3s` is the stored dataset name; interpret the physical units
according to the reaction order described in [Reactions](reactions.md).

The actual species state over time is stored through `trajectory/species_id_indices` in full mode or final `species_id_indices` in minimal mode.

---

## `trajectory_mode = "minimal"`

Minimal mode is useful when you do not need full time series trajectories. It keeps the HDF5 metadata but writes compact final-state data under:

```text
/analysis/minimal_transport/
```

Important datasets:

| Dataset | Shape | Meaning |
|---|---:|---|
| `final_pos_x`, `final_pos_y`, `final_pos_z` | `[N]` | Final position components. |
| `final_vel_x`, `final_vel_y`, `final_vel_z` | `[N]` | Final velocity components. |
| `ion_time_s` | `[N]` | Final per-ion time. |
| `birth_time_s` | `[N]` | Ion birth time. |
| `death_time_s` | `[N]` | Ion loss/death time. |
| `domain_index` | `[N]` | Final domain index. |
| `species_id_indices` | `[N]` | Final species index. |
| `species_pool` | `[S]` | Species ID lookup table when available. |
| `active` | `[N]` | Whether the ion is active at the end. |
| `born` | `[N]` | Whether the ion had been born/activated. |

Minimal mode is not a drop-in replacement for arrival time or diffusion analysis because intermediate trajectory frames are omitted.

---

## Quick inspection with Python

Install `h5py` in your analysis environment:

```bash
python3 -m pip install h5py numpy pandas matplotlib
```

List top-level groups:

```python
import h5py

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    print(list(f.keys()))
```

Print all datasets and shapes:

```python
import h5py

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    def show(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(name, obj.shape, obj.dtype)

    f.visititems(show)
```

Load positions and time in full mode:

```python
import h5py
import numpy as np

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    t = f["trajectory/time"][:]
    pos = f["trajectory/positions"][:]   # [T, N, 3]

z = pos[:, :, 2]
print("time frames:", t.shape[0])
print("ions:", z.shape[1])
print("final z range:", np.nanmin(z[-1]), np.nanmax(z[-1]))
```

Load final states in minimal mode:

```python
import h5py

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    minimal = f["analysis/minimal_transport"]
    z_final = minimal["final_pos_z"][:]
    active = minimal["active"][:].astype(bool)

print("ions:", z_final.size)
print("active fraction:", active.mean())
```

See [Analysis workflow](analysis.md) for more complete post-processing examples.

---

## Practical interpretation checklist

Before interpreting physics, check:

- `/metadata/completion/success` indicates a successful run.
- `/metadata/config/config_json` matches the configuration you intended to run.
- `/metadata/reproducibility/global_seed` is recorded.
- `/metadata/reproducibility/git_hash` and `git_dirty` are known.
- `/metadata/reproducibility/input_hash/` contains the expected config/species/reaction hashes.
- `/domains/` contains the intended geometry, gas, and field settings.
- `/ions/initial_species_id` contains the species you expected.
- `death_time_s` is not showing immediate losses from bad initial positions or boundaries.
- `trajectory/time` spacing is fine enough for the analysis you want to perform.

These checks help distinguish a physics/model effect from a configuration, database, or output-sampling issue.

# Analysis workflow

This page shows a practical first analysis workflow for ICARION HDF5 [output](output-files.md). The goal is to help new users answer basic questions:

- Did the run finish properly and write the expected datasets?
- How many ions survived or were transmitted?
- What were the final positions and species?
- How can arrival times, drift velocities, or simple mobility estimates be extracted?

---

## Recommended Python environment

Create a lightweight analysis environment:

```bash
python3 -m venv .venv-analysis
source .venv-analysis/bin/activate
pip install h5py numpy pandas matplotlib
```

For notebooks, also install:

```bash
pip install jupyter
```

---

## Inspect a file

Start by listing groups and datasets:

```python
from pathlib import Path
import h5py

path = Path("results/ims/ims_trajectories.h5")

with h5py.File(path, "r") as f:
    def show(name, obj):
        if isinstance(obj, h5py.Dataset):
            print(name, obj.shape, obj.dtype)
        else:
            print(name + "/")

    f.visititems(show)
```

For large files, inspect only the top-level groups first:

```python
import h5py

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    print(list(f.keys()))
    print(list(f["trajectory"].keys()))
```

---

## Load trajectory arrays

```python
import h5py
import numpy as np

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    time = f["trajectory/time"][:]             # [T]
    pos = f["trajectory/positions"][:]         # [T, N, 3]
    vel = f["trajectory/velocities"][:]        # [T, N, 3]

x = pos[:, :, 0]
y = pos[:, :, 1]
z = pos[:, :, 2]

print("frames:", time.size)
print("ions:", z.shape[1])
print("final z min/max:", np.nanmin(z[-1]), np.nanmax(z[-1]))
```

For very large files, avoid loading everything at once. Slice only what you need:

```python
with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    time = f["trajectory/time"][:]
    z_final = f["trajectory/positions"][-1, :, 2]
```

---

## Final-state table

A convenient first analysis product is one row per ion:

```python
import h5py
import pandas as pd

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    final_pos = f["trajectory/positions"][-1, :, :]
    final_vel = f["trajectory/velocities"][-1, :, :]

    data = {
        "x_m": final_pos[:, 0],
        "y_m": final_pos[:, 1],
        "z_m": final_pos[:, 2],
        "vx_m_s": final_vel[:, 0],
        "vy_m_s": final_vel[:, 1],
        "vz_m_s": final_vel[:, 2]
    }

    if "domain_indices" in f["trajectory"]:
        data["domain_index"] = f["trajectory/domain_indices"][-1, :]

    if "ions/death_time_s" in f:
        data["death_time_s"] = f["ions/death_time_s"][:]

    df = pd.DataFrame(data)

print(df.head())
print(df.describe())
```

Save the compact table:

```python
df.to_csv("results/ims/final_state.csv", index=False)
```

---

## Transmission or survival estimate

A simple survival estimate can be computed from `death_time_s` if available:

```python
import h5py
import numpy as np

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    death = f["ions/death_time_s"][:]

alive = death < 0
print("surviving fraction:", alive.mean())
print("lost ions:", np.count_nonzero(~alive))
```

For an IMS domain, transmission can also be defined geometrically, e.g. ions that reached a final axial coordinate:

```python
import h5py
import numpy as np

L = 0.10  # drift length in m; use your configured domain length

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    z_final = f["trajectory/positions"][-1, :, 2]

transmitted = z_final >= L
print("geometric transmission:", transmitted.mean())
```

Choose the transmission definition that matches the simulated instrument and boundary conditions.

---

## Arrival-time estimate from trajectories

For IMS simulations, an approximate arrival time is the first stored time at which an ion reaches the end of the drift region.

```python
import h5py
import numpy as np

L = 0.10  # drift length [m]

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    time = f["trajectory/time"][:]
    z = f["trajectory/positions"][:, :, 2]

arrived = z >= L
arrival_time = np.full(z.shape[1], np.nan)

for ion_idx in range(z.shape[1]):
    hits = np.flatnonzero(arrived[:, ion_idx])
    if hits.size:
        arrival_time[ion_idx] = time[hits[0]]

valid = np.isfinite(arrival_time)
print("arrived:", valid.sum(), "/", arrival_time.size)
print("mean arrival time [s]:", np.nanmean(arrival_time))
print("std arrival time [s]:", np.nanstd(arrival_time))
```

This estimate is limited by `write_interval`. If accurate arrival times are important, write output often enough or use a dedicated arrival time diagnostic if available.

---

## Simple mobility estimate

For a uniform field drift region, a first-pass mobility estimate is:

```text
K = L / (E * t_d)
```

where:

- `L` is the drift length,
- `E` is the electric field magnitude,
- `t_d` is the drift/arrival time.

Example:

```python
import numpy as np

L = 0.10          # m
E = 500.0         # V/m, set this from your config

t_mean = np.nanmean(arrival_time)
K = L / (E * t_mean)

print("K [m^2/(V s)]:", K)
print("K [cm^2/(V s)]:", K * 1e4)
```

For reduced mobility `K0`, apply the pressure and temperature correction used in your simulation.

---

## Species-resolved analysis

If species indices are written, you can separate ions by final species. The exact mapping from species index to species name is stored in metadata/species or compatibility datasets, depending on the output version. See also [species database](species-database.md).

A robust workflow is:

1. inspect available datasets under `/metadata/species/`,
2. inspect `trajectory/species_id_indices`,
3. create a mapping table once for the file,
4. join the mapping to the final-state table.

Start by listing species metadata:

```python
import h5py

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    if "metadata/species" in f:
        print(list(f["metadata/species"].keys()))
    if "trajectory/species_id_indices" in f:
        print(f["trajectory/species_id_indices"].shape)
```

For reaction simulations, plotting the number of ions per species versus time is often the first sanity check.

---

## Plot an axial trajectory subset

```python
import h5py
import matplotlib.pyplot as plt

with h5py.File("results/ims/ims_trajectories.h5", "r") as f:
    time = f["trajectory/time"][:]
    z = f["trajectory/positions"][:, :20, 2]  # first 20 ions

plt.figure()
plt.plot(time, z)
plt.xlabel("time / s")
plt.ylabel("z / m")
plt.tight_layout()
plt.show()
```

---

## Practical checklist

Before interpreting physics, check:

- The HDF5 file contains the expected number of ions.
- The stored `config_json` matches the config you intended to run.
- The git hash and dirty flag are known.
- Ions did not immediately hit boundaries because of a geometry or initial-position mistake.
- The arrival time estimate is not limited by a `write_interval` that is too large.
- Reactions did not stay disabled accidentally.
- Species IDs match the species database.

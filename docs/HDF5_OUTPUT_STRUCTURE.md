# HDF5 Output Structure

**Version:** 1.0.0  
**Last Updated:** December 2025  
**Status:** Implemented in v1.0.0; waveform library, config/species/reaction DBs, and field arrays are embedded in HDF5 (config snapshot also written alongside output).  
**Implementation:** Writer is SoA-native (IonEnsemble) in v1.0.0.

---

## Overview

ICARION writes simulation results to HDF5 files with a hierarchical structure designed for:
- **Reproducibility**: Complete metadata for scientific publications
- **Analysis**: Pandas/Python-friendly tabular formats
- **Performance**: Efficient chunked storage with compression
- **Documentation**: Self-describing with units and descriptions

---

## File Structure

```
simulation.h5
├── metadata/
│   ├── config/                    # Selected configuration parameters
│   ├── physics/                   # Physics handler metadata
│   ├── reproducibility/           # Git hash, RNG seed, build info
│   ├── system/                    # System information
│   ├── species/                   # Species database (tabular)
│   ├── reactions/                 # Reaction database (tabular)
│   └── completion/                # Simulation completion status
├── trajectory/
│   ├── time                       # Time points [T]
│   ├── positions                  # Ion positions [T × N × 3]
│   ├── velocities                 # Ion velocities [T × N × 3]
│   ├── domain_indices             # Domain index per ion [T × N]
│   ├── species_id_indices         # Species pool indices [T × N]
│   ├── time_per_ion               # Per-ion times [T × N] (if adaptive dt)
│   └── species_ids                # Species per ion [T × N] (optional compatibility)
├── ions/
│   ├── initial_species_id         # Species name [N]
│   ├── initial_pos_x/y/z          # Initial position [N]
│   ├── initial_vel_x/y/z          # Initial velocity [N]
│   ├── birth_time_s               # Birth time [N]
│   ├── death_time_s               # Death time [N] (-1 if still alive)
│   └── charge_C                   # Charge [N]
└── domains/
    ├── domain_0/
    │   ├── name                   # Domain name
    │   ├── instrument             # Instrument type
    │   ├── solver                 # Solver type
    │   ├── geometry/              # Geometric parameters
    │   ├── environment/           # Gas conditions
    │   ├── fields/                # Electric fields
    │   │   ├── waveforms/         # Waveform library (if present)
    │   │   ├── dc/                # DC field values (t=0 if waveform)
    │   │   ├── rf/                # RF field values (t=0 if waveform)
    │   │   └── ac/                # AC field values (t=0 if waveform)
    └── domain_1/
        └── ...
└── analysis/                      # Optional diagnostics groups
    └── minimal_transport/         # Optional compact final-state output (trajectory_mode=minimal)
```

---

## Metadata Group

### `/metadata/config/`

Configuration parameters extracted from `FullConfig`.

**Datasets:**

| Name | Type | Description | Units |
|------|------|-------------|-------|
| `format_version` | string | HDF5 format version | - |
| `dt_s` | double | Integration timestep | s |
| `total_time_s` | double | Total simulation time | s |
| `total_steps` | int | Number of integration steps | - |
| `write_interval` | int | Output write interval | steps |
| `integrator` | string | Integration method (RK4, RK45, Boris) | - |
| `collision_model` | string | Collision model (EHSS, HSS, Langevin) | - |
| `enable_reactions` | bool | Reactions enabled? | - |
| `enable_space_charge` | bool | Space charge enabled? | - |
| `enable_space_charge_gpu` | bool | Space charge GPU requested? | - |
| `enable_gpu` | bool | GPU acceleration enabled? | - |
| `output_file` | string | Output file path | - |
| `trajectory_mode` | string | Trajectory output mode (`full` or `minimal`) | - |
| `config_json` | string | Embedded resolved config snapshot (validated + CLI overrides) | - |
| `integrator_params/name` | string | Integrator name | - |
| `integrator_params/rk45_min_step_s` | double | Absolute min step for RK45 (if set) | s |
| `integrator_params/openmp_enabled` | bool | OpenMP enabled flag | - |
| `integrator_params/rk45_atol` | double | RK45 absolute tolerance (default if RK45) | - |
| `integrator_params/rk45_rtol` | double | RK45 relative tolerance (default if RK45) | - |
| `integrator_params/rk45_safety_factor` | double | RK45 safety factor (default) | - |
| `integrator_params/rk45_min_step_factor` | double | RK45 minimum step scaling factor | - |
| `integrator_params/rk45_max_step_factor` | double | RK45 maximum step scaling factor | - |
| `integrator_params/rk45_max_step_increase` | double | RK45 max growth per step | - |
| `integrator_params/rk45_max_step_decrease` | double | RK45 max shrink per step | - |
| `integrator_params/rk45_absolute_min_step_s` | double | RK45 absolute min step (config/strategy default) | s |
| `integrator_params/gpu_collision_threshold` | int | Minimum ions for GPU collision dispatch (default 5000) | - |
| `integrator_params/gpu_space_charge_threshold` | int | Minimum ions for GPU space-charge dispatch (default 1000, 0 if GPU disabled) | - |
| `derived_summary/num_domains` | int | Number of domains in config | - |
| `derived_summary/num_species_db` | int | Species entries loaded | - |
| `derived_summary/num_reactions_db` | int | Reactions loaded | - |

**Note:** `format_version` is currently `2.0.0` in the v1.0.0 release. The full resolved config snapshot is embedded as `config_json` (also written as a `.config.json` snapshot alongside the HDF5). Selected key fields remain broken out as individual datasets for fast access.

---

### `/metadata/physics/`

Physics handler metadata stored separately from config.

**Datasets:**

| Name | Type | Description |
|------|------|-------------|
| `collision_handler` | string | Collision handler name |
| `reaction_handler` | string | Reaction handler name (or "None") |
| `reaction_gpu_threshold` | int | Minimum ions for GPU reaction dispatch |
| `collision_gpu_threshold` | int | Minimum ions for GPU collision dispatch |
| `collision_mixture_limit` | string | Max mixture components in GPU helper |

---

### `/metadata/reproducibility/`

Information required to reproduce the simulation exactly.

**Datasets:**

| Name | Type | Description |
|------|------|-------------|
| `global_seed` | uint | RNG seed for entire simulation |
| `rng_algorithm` | string | RNG algorithm (std::mt19937_64) |
| `seed_scheme` | string | Per-ion seeding method |
| `per_ion_rng_scope` | string | Subsystems using per-ion RNG (collisions, reactions, stochastic forces) |
| `git_hash` | string | Git commit hash of ICARION |
| `git_dirty` | bool | Uncommitted changes present? |
| `code_version` | string | ICARION version (e.g., "1.0.0") |
| `build_type` | string | Release or Debug |
| `compiler_cxx` | string | C++ compiler version |
| `build_info` | string | Build flags summary (type, OpenMP, CUDA, mode) |
| `cuda_version` | string | CUDA version (if enabled) |
| `openmp_enabled` | bool | OpenMP enabled at build time? |
| `openmp_threads` | int | Number of OpenMP threads (1 if disabled) |

**Subgroup:** `/metadata/reproducibility/input_hash/`

| Name | Type | Description |
|------|------|-------------|
| `config_sha256` | string | SHA256 hash of config file |
| `species_db_sha256` | string | SHA256 hash of species database (N/A if not provided) |
| `reaction_db_sha256` | string | SHA256 hash of reaction database (N/A if not provided) |
| `field_arrays/files` | string[] | Field array file paths |
| `field_arrays/sha256` | string[] | SHA256 hashes for field array files |

---

**Subgroup:** `/metadata/reproducibility/input_blobs/`

Embedded copies of external inputs (for reruns without external files).

| Name | Type | Description |
|------|------|-------------|
| `config_json` | string | Full config JSON (or `{}` if unavailable) |
| `species_db_json` | string | Embedded species database (or `{}`) |
| `reaction_db_json` | string | Embedded reaction database (or `{}`) |
| `field_arrays/blob_<i>` | uint8 array | Raw bytes of field array file *i* (if embedded) |
| `field_arrays/blob_<i>_filename` | string | Original path of field array file *i* |

When source paths are available, inputs are embedded for reproducibility. Missing inputs are stored as `{}` or empty blobs.

---

### `/metadata/system/`

System information where simulation was executed.

**Datasets:**

| Name | Type | Description |
|------|------|-------------|
| `hostname` | string | Machine hostname |
| `username` | string | User who ran simulation |
| `os` | string | Operating system and version |
| `kernel` | string | Kernel version |
| `cpu_model` | string | CPU model name |
| `cpu_cores` | int | Number of CPU cores |
| `memory_gb` | double | Total RAM in GB |
| `gpu_model` | string | GPU model (if CUDA enabled) |
| `gpu_memory_gb` | double | GPU memory in GB |
| `driver_version` | string | CUDA driver version |
| `timestamp` | string | Execution timestamp (ISO 8601) |

---

### `/metadata/species/`

Species database in tabular format (pandas-compatible).

**Datasets:**

| Name | Type | Shape | Description | Units |
|------|------|-------|-------------|-------|
| `names` | string | [S] | Species names (e.g., "H3O+") | - |
| `mass_kg` | double | [S] | Particle mass | kg |
| `charge_C` | double | [S] | Particle charge | C |
| `mobility_m2Vs` | double | [S] | Reduced mobility | m²/(V·s) |
| `ccs_m2` | double | [S] | Collision cross section | m² |

**S** = Number of unique species

**Note:** Only species referenced by ions in the simulation are written (not the entire database). This reduces file size for simulations using large species databases.

**Python Example:**
```python
import h5py
import pandas as pd

with h5py.File('simulation.h5', 'r') as f:
    species_df = pd.DataFrame({
        'name': f['/metadata/species/names'][:].astype(str),
        'mass_kg': f['/metadata/species/mass_kg'][:],
        'charge_C': f['/metadata/species/charge_C'][:],
        'mobility_m2Vs': f['/metadata/species/mobility_m2Vs'][:],
        'ccs_m2': f['/metadata/species/ccs_m2'][:]
    })
```

---

### `/metadata/reactions/`

Reaction database in tabular format.

**Datasets:**

| Name | Type | Shape | Description | Units |
|------|------|-------|-------------|-------|
| `id` | string | [R] | Reaction identifier | - |
| `reactant_1` | string | [R] | First reactant species | - |
| `reactant_2` | string | [R] | Second reactant species | - |
| `product_1` | string | [R] | First product species | - |
| `rate_constant_m3s` | double | [R] | Reaction rate constant | m³/s |
| `type` | int | [R] | Reaction type enum | - |

**R** = Number of reactions

**Note:** Only reactions involving species present in the simulation are written (reactants must be in initial ions). This reduces file size for large reaction networks (e.g., 1000+ reactions).

**Reaction Types:**
- `0` = ChargeTransfer
- `1` = ProtonTransfer
- `2` = Association (two-body)
- `3` = Dissociation
- `4` = Switching
- `5` = Unknown

**Note:** Current implementation writes type `2` (Association) for all reactions; type-specific IDs are not yet persisted.

### `/metadata/completion/`

Simulation completion status and final state.

**Datasets:**

| Name | Type | Shape | Description | Units |
|------|------|-------|-------------|-------|
| `success` | int | [1] | Simulation completed successfully (1=yes, 0=no) | - |
| `final_time_s` | double | [1] | Final simulation time reached | s |
| `active_ions` | int | [1] | Number of ions active at completion | - |
| `completion_timestamp` | string | [1] | ISO 8601 timestamp when completed | - |

**Example:**
```python
with h5py.File('simulation.h5', 'r') as f:
    success = f['/metadata/completion/success'][0]
    final_time = f['/metadata/completion/final_time_s'][0]
    active_ions = f['/metadata/completion/active_ions'][0]
    timestamp = f['/metadata/completion/completion_timestamp'][0].decode('utf-8')
```

---

## Trajectory Group

### `/trajectory/`

Time-series data for all ions at each output timestep.

When `output.trajectory_mode = "minimal"`, this group is created but left empty; compact per-ion final-state data is written to `/analysis/minimal_transport` instead.

**Datasets:**

| Name | Type | Shape | Description | Units |
|------|------|-------|-------------|-------|
| `time` | double | [T] | Time points | s |
| `positions` | double | [T × N × 3] | Ion positions (x, y, z) | m |
| `velocities` | double | [T × N × 3] | Ion velocities (vx, vy, vz) | m/s |
| `domain_indices` | int | [T × N] | Domain index per ion | - |
| `species_id_indices` | uint32 | [T × N] | Species pool index per ion | - |
| `time_per_ion` | double | [T × N] | Per-ion times (adaptive dt snapshots) | s |
| `species_ids` | string | [T × N] | Species identifier per ion (optional compatibility) | - |

**T** = Number of timesteps written  
**N** = Number of ions

**Storage:**
- Chunked storage for efficient partial reads
- GZIP compression level 2 for buffered batch writes (level 6 for single-step append)
- Typical compression ratio: 3-5x

**Python Example:**
```python
import h5py
import numpy as np

with h5py.File('simulation.h5', 'r') as f:
    time = f['/trajectory/time'][:]
    positions = f['/trajectory/positions'][:]  # Shape: [T, N, 3]
    
    # Get ion 5 trajectory
    ion_5_x = positions[:, 5, 0]
    ion_5_y = positions[:, 5, 1]
    ion_5_z = positions[:, 5, 2]
    
    # Get all positions at timestep 100
    snapshot = positions[100, :, :]  # Shape: [N, 3]
```

---

## Ions Group

### `/ions/`

Per-ion metadata (initial conditions and properties).

**Datasets:**

| Name | Type | Shape | Description | Units |
|------|------|-------|-------------|-------|
| `initial_species_id` | string | [N] | Species name | - |
| `initial_pos_x` | double | [N] | Initial x position | m |
| `initial_pos_y` | double | [N] | Initial y position | m |
| `initial_pos_z` | double | [N] | Initial z position | m |
| `initial_vel_x` | double | [N] | Initial x velocity | m/s |
| `initial_vel_y` | double | [N] | Initial y velocity | m/s |
| `initial_vel_z` | double | [N] | Initial z velocity | m/s |
| `birth_time_s` | double | [N] | Birth time (0 for initial ions) | s |
| `death_time_s` | double | [N] | Death time (-1 if ion still alive at end) | s |
| `charge_C` | double | [N] | Particle charge | C |

**N** = Number of ions

**Python Example:**
```python
# Find all H3O+ ions
with h5py.File('simulation.h5', 'r') as f:
    species_ids = f['/ions/initial_species_id'][:].astype(str)
    h3o_indices = np.where(species_ids == 'H3O+')[0]
    
    h3o_initial_pos_x = f['/ions/initial_pos_x'][h3o_indices]
```

---

## Domains Group

### `/domains/domain_<i>/`

Configuration for each instrument domain.

**Attributes:**

| Name | Type | Description |
|------|------|-------------|
| `name` | string | Domain name |
| `instrument` | string | Instrument type (IMS, LQIT, TOF, etc.) |
| `solver` | string | Solver type (RK4, RK45, Boris) |

**Subgroup:** `/domains/domain_<i>/geometry/`

| Name | Type | Description | Units |
|------|------|-------------|-------|
| `length_m` | double | Domain length | m |
| `radius_m` | double | Domain radius | m |
| `radius_in_m` | double | Inner radius (for cylindrical) | m |
| `radius_out_m` | double | Outer radius (for cylindrical) | m |
| `origin_m` | double[3] | Domain origin (x, y, z) | m |

**Subgroup:** `/domains/domain_<i>/environment/`

| Name | Type | Description | Units |
|------|------|-------------|-------|
| `pressure_Pa` | double | Gas pressure | Pa |
| `temperature_K` | double | Gas temperature | K |
| `gas_species` | string | Buffer gas species | - |
| `particle_density_m3` | double | Number density | m⁻³ |
| `mean_thermal_velocity_ms` | double | Mean thermal velocity | m/s |
| `gas_velocity_ms` | double[3] | Gas flow velocity | m/s |

**Subgroup:** `/domains/domain_<i>/fields/`

Contains subgroups for DC, RF, and AC fields with voltage/frequency parameters. RF includes `phase_rad`; AC has no phase field.

**Waveform Support:** If time-varying waveforms are used, a `/fields/waveforms/` subgroup is created with one group per named or inline waveform containing its type and parameters. DC/RF/AC fields store the t=0 evaluation for backward compatibility.

---

## Analysis Group (Optional)

### `/analysis/minimal_transport/`

Compact per-ion final-state output written at finalize when
`output.trajectory_mode = "minimal"`.

**Datasets:**

| Name | Type | Shape | Description | Units |
|------|------|-------|-------------|-------|
| `final_pos_x` | double | [N] | Final x position | m |
| `final_pos_y` | double | [N] | Final y position | m |
| `final_pos_z` | double | [N] | Final z position | m |
| `final_vel_x` | double | [N] | Final x velocity | m/s |
| `final_vel_y` | double | [N] | Final y velocity | m/s |
| `final_vel_z` | double | [N] | Final z velocity | m/s |
| `domain_index` | int32 | [N] | Final domain index per ion | - |
| `species_id_indices` | uint32 | [N] | Final species pool index per ion | - |
| `species_pool` | string | [S] | Species name pool for `species_id_indices` | - |
| `active` | uint8 | [N] | Final active flag (1/0) | - |
| `born` | uint8 | [N] | Born flag (1/0) | - |
| `ion_time_s` | double | [N] | Final per-ion simulation time | s |
| `birth_time_s` | double | [N] | Per-ion birth time | s |
| `death_time_s` | double | [N] | Per-ion death time (`-1` if alive) | s |

---
## Python Analysis Examples

### Load Complete Dataset

```python
import h5py
import pandas as pd
import numpy as np

class IcarionHDF5:
    """Helper class to read ICARION HDF5 outputs"""
    
    def __init__(self, filename):
        self.filename = filename
        self.file = h5py.File(filename, 'r')
    
    def __enter__(self):
        return self
    
    def __exit__(self, *args):
        self.file.close()
    
    @property
    def config(self):
        """Get configuration parameters"""
        cfg = self.file['/metadata/config']
        return {key: cfg[key][()] for key in cfg.keys()}
    
    @property
    def species(self):
        """Get species database as DataFrame"""
        sp = self.file['/metadata/species']
        return pd.DataFrame({
            'name': sp['names'][:].astype(str),
            'mass_kg': sp['mass_kg'][:],
            'charge_C': sp['charge_C'][:],
            'mobility_m2Vs': sp['mobility_m2Vs'][:],
            'ccs_m2': sp['ccs_m2'][:]
        })
    
    def get_trajectory(self, ion_index=None):
        """Get trajectory for specific ion or all ions"""
        traj = self.file['/trajectory']
        time = traj['time'][:]
        pos = traj['positions'][:]
        vel = traj['velocities'][:]
        
        if ion_index is not None:
            pos = pos[:, ion_index, :]
            vel = vel[:, ion_index, :]
        
        return time, pos, vel

# Usage
with IcarionHDF5('simulation.h5') as sim:
    print(f"Timestep: {sim.config['dt_s']} s")
    print(f"Species: {sim.species}")
    
    time, pos, vel = sim.get_trajectory(ion_index=0)
    print(f"Ion 0 final position: {pos[-1]}")
```

### Arrival Time Analysis

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

with h5py.File('simulation.h5', 'r') as f:
    time = f['/trajectory/time'][:]
    positions = f['/trajectory/positions'][:]
    
    # Detect when ions reach z = 0.05 m (detector)
    z_positions = positions[:, :, 2]
    detector_z = 0.05
    
    arrival_times = []
    for ion_idx in range(z_positions.shape[1]):
        z_traj = z_positions[:, ion_idx]
        crossed = np.where(z_traj >= detector_z)[0]
        if len(crossed) > 0:
            arrival_times.append(time[crossed[0]])
    
    # Plot histogram
    plt.hist(arrival_times, bins=50)
    plt.xlabel('Arrival Time (s)')
    plt.ylabel('Count')
    plt.title('Ion Mobility Spectrum')
    plt.show()
```

---

## File Size Estimates

Typical file sizes for various simulations (numeric-only trajectory datasets; excludes optional `species_ids` strings):

| Simulation | Ions | Timesteps | Raw (numeric only) | Compressed (3-5x) |
|------------|------|-----------|--------------------|-------------------|
| IMS small | 100 | 1,000 | ~6.4 MB | ~1.3-2.1 MB |
| IMS medium | 1,000 | 5,000 | ~320 MB | ~64-107 MB |
| LQIT | 10,000 | 10,000 | ~6.4 GB | ~1.3-2.1 GB |

**Storage formula (trajectory data only, numeric datasets):**
```
size_bytes ≈ N_ions * N_timesteps * (6*8 + 8 + 4 + 4) + N_timesteps * 8
```

where:
- 6 values per ion per timestep: pos[3], vel[3] (doubles)
- time_per_ion per ion per timestep (double)
- domain_indices (int32) and species_id_indices (uint32)
- 1 time value per timestep (double)
- Typical compression ratio: 3-5x for numeric datasets
- `species_ids` (varlen strings) can add significant overhead when written
- **Note:** Metadata (config, species, reactions, domains) adds ~1-2 MB overhead

---

## Compatibility

### HDF5 Version
- **Minimum:** HDF5 1.10.0
- **Recommended:** HDF5 1.12.0+
- **Compression:** GZIP level 2 for buffered batch writes, level 6 for single-step append

### Python Libraries
```bash
pip install h5py pandas numpy matplotlib
```

### Reading in Other Languages

**C++:**
```cpp
#include <H5Cpp.h>

H5::H5File file("simulation.h5", H5F_ACC_RDONLY);
H5::DataSet ds_time = file.openDataSet("/trajectory/time");
```

**MATLAB:**
```matlab
time = h5read('simulation.h5', '/trajectory/time');
positions = h5read('simulation.h5', '/trajectory/positions');
```

**Julia:**
```julia
using HDF5

h5open("simulation.h5", "r") do file
    time = read(file, "/trajectory/time")
    positions = read(file, "/trajectory/positions")
end
```

---

## Best Practices

### For Users

1. **Always check `git_dirty` flag** - If true, results may not be reproducible
2. **Store input files with outputs** - Use SHA256 hashes to verify
3. **Document post-processing** - Add custom groups for analysis results
4. **Use compression** - Saves 70-80% disk space with minimal performance cost

### For Developers

1. **Never remove datasets** - Only add new ones (backward compatibility)
2. **Use attributes for metadata** - Keep datasets for array data
3. **Write incrementally** - Use chunked storage for large datasets
4. **Test with h5dump** - Verify structure before committing

---

## References

- HDF5 Documentation: https://portal.hdfgroup.org/
- h5py User Guide: https://docs.h5py.org/
- ICARION Config Schema: [`CONFIG_GUIDE.md`](CONFIG_GUIDE.md)

---

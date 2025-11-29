# HDF5 Output Structure

**Version:** 1.0 (FullConfig-based)  
**Last Updated:** 2025-11-21  
**Status:** Implemented

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
│   ├── config/                    # Configuration parameters
│   ├── reproducibility/           # Git hash, RNG seed, build info
│   ├── system/                    # System information
│   ├── species/                   # Species database (tabular)
│   ├── reactions/                 # Reaction database (tabular)
│   └── completion/                # Simulation completion status
├── trajectory/
│   ├── time                       # Time points [T]
│   ├── positions                  # Ion positions [T × N × 3]
│   ├── velocities                 # Ion velocities [T × N × 3]
│   ├── species_ids                # Species per ion [T × N]
│   └── domain_indices             # Domain index per ion [T × N]
├── ions/
│   ├── initial_species_id         # Species name [N]
│   ├── initial_pos_x/y/z          # Initial position [N]
│   ├── initial_vel_x/y/z          # Initial velocity [N]
│   ├── birth_time_s               # Birth time [N]
│   └── charge_C                   # Charge [N]
└── domains/
    ├── domain_0/
    │   ├── name                   # Domain name
    │   ├── instrument             # Instrument type
    │   ├── solver                 # Solver type
    │   ├── geometry/              # Geometric parameters
    │   ├── environment/           # Gas conditions
    │   ├── fields/                # Electric fields
    │   │   ├── waveforms/         # Waveform library (v1.1, if present)
    │   │   ├── dc/                # DC field values (t=0 if waveform)
    │   │   ├── rf/                # RF field values (t=0 if waveform)
    │   │   └── ac/                # AC field values (t=0 if waveform)
    └── domain_1/
        └── ...
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
| `enable_gpu` | bool | GPU acceleration enabled? | - |
| `output_file` | string | Output file path | - |

**Note:** Complete JSON configuration may be added in future versions.

---

### `/metadata/reproducibility/`

Information required to reproduce the simulation exactly.

**Datasets:**

| Name | Type | Description |
|------|------|-------------|
| `global_seed` | uint | RNG seed for entire simulation |
| `rng_algorithm` | string | RNG algorithm (std::mt19937_64) |
| `seed_scheme` | string | Per-ion seeding method |
| `git_hash` | string | Git commit hash of ICARION |
| `git_dirty` | bool | Uncommitted changes present? |
| `code_version` | string | ICARION version (e.g., "1.0.0") |
| `build_type` | string | Release or Debug |
| `compiler_cxx` | string | C++ compiler version |
| `build_info` | string | Complete build information |
| `cuda_version` | string | CUDA version (if enabled) |
| `gpu_arch` | string | GPU architecture (if enabled) |
| `openmp_threads` | int | Number of OpenMP threads |

**Subgroup:** `/metadata/reproducibility/input_hash/`

| Name | Type | Description |
|------|------|-------------|
| `config_sha256` | string | SHA256 hash of config file |
| `species_db_sha256` | string | SHA256 hash of species database |
| `reaction_db_sha256` | string | SHA256 hash of reaction database |

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
- `0` = Three-body
- `1` = Charge transfer
- `2` = Proton transfer

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

**Datasets:**

| Name | Type | Shape | Description | Units |
|------|------|-------|-------------|-------|
| `time` | double | [T] | Time points | s |
| `positions` | double | [T × N × 3] | Ion positions (x, y, z) | m |
| `velocities` | double | [T × N × 3] | Ion velocities (vx, vy, vz) | m/s |
| `species_ids` | string | [T × N] | Species identifier per ion | - |
| `domain_indices` | int | [T × N] | Domain index per ion | - |

**T** = Number of timesteps written  
**N** = Number of ions

**Storage:**
- Chunked storage for efficient partial reads
- GZIP compression level 6
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
| `birth_time_s` | double | [N] | Birth time | s |
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

Contains subgroups for DC, RF, and AC fields with voltage, frequency, and phase parameters. 

**v1.1 Extension:** If time-varying waveforms are used, a `/fields/waveforms/` subgroup is created with one group per named waveform containing its type and parameters. DC/RF/AC fields store the t=0 evaluation for backward compatibility.

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

Typical file sizes for various simulations:

| Simulation | Ions | Timesteps | File Size | Compressed |
|------------|------|-----------|-----------|------------|
| IMS small | 100 | 1,000 | 2.4 MB | 0.6 MB |
| IMS medium | 1,000 | 5,000 | 120 MB | 30 MB |
| LQIT | 10,000 | 10,000 | 2.4 GB | 600 MB |
| TOF large | 100,000 | 50,000 | 120 GB | 30 GB |

**Storage formula:**
```
size_bytes ≈ (N_ions × N_timesteps × 7 × 8) / compression_ratio
```

where:
- 7 datasets per timestep (time, pos[3], vel[3])
- 8 bytes per double
- Typical compression ratio: 3-5x

---

## Compatibility

### HDF5 Version
- **Minimum:** HDF5 1.10.0
- **Recommended:** HDF5 1.12.0+
- **Compression:** GZIP level 6

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

## Format Versioning

**Current Version:** 1.0

### Version History

| Version | Date | Changes |
|---------|------|---------|
| 1.0 | 2024-10 | Initial format (legacy GlobalParams) |
| 2.0 | 2025-11 | FullConfig-based, added metadata groups |

### Version Detection

```python
with h5py.File('simulation.h5', 'r') as f:
    if '/metadata/config/format_version' in f:
        version = f['/metadata/config/format_version'][()].decode()
        print(f"HDF5 format version: {version}")
    else:
        print("Legacy format (v1.0)")
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
- ICARION Config Schema: [`docs/INPUT_SCHEMA.md`](docs/INPUT_SCHEMA.md )
- Example Scripts: `examples/analyze_trajectory.py`

---

**Questions or Issues?**  
https://github.com/ICARION-Project/ICARION/issues

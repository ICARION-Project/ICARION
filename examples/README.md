# ICARION Example Configurations

This directory contains production-ready example JSON configuration  
files for **ICARION Core v1.0**.

All examples follow the **v1.0 schema** and are fully functional.

## Quick Start

Run any example with:

```bash
cd /home/chsch95/ICARION
./build/icarion_core examples/<example_name>.json
```

Validate a configuration:

```bash
python schema/validate_config.py examples/<example_name>.json
```

## Available Examples

### Basic Instrument Examples

#### `ims_basic.json` – Ion Mobility Spectrometer

- **Purpose**: Basic IMS drift separation
- **Ion Count**: 50,000 ions (IonA+, 150 amu, charge +1)
- **Physics**: EHSS collisions, N₂ buffer gas (200 Pa)
- **Runtime**: ~15 seconds

```bash
./build/icarion_core examples/ims_basic.json
```

#### `lqit_basic.json` – Linear Quadrupole Ion Trap

- **Purpose**: Basic Paul trap with RF confinement
- **Ion Count**: 2,000 ions (Ca+)
- **Physics**: RF quadrupole fields, no collisions
- **Runtime**: ~10 seconds

```bash
./build/icarion_core examples/lqit_basic.json
```

#### `orbitrap_basic.json` – Orbitrap Mass Analyzer

- **Purpose**: Orbitrap ion confinement and axial oscillation
- **Ion Count**: 100 ions (Reserpine+, 609.28 Da)
- **Physics**: Hyperbolic electrodes (r_in=10mm, r_out=20mm, r_char=15mm),  
5kV radial voltage
- **Key Features**:
  - Hyperbolic boundary geometry: z² = 0.5(r² - R²) - R_m² ln(r/R)
  - Characteristic axial oscillation frequency
  - Optional voltage sweep for resonant excitation
- **Runtime**: ~5 seconds
- **Applications**: High-resolution mass spectrometry, axial frequency measurements

```bash
./build/icarion_core examples/orbitrap_basic.json
```

**Configuration Highlights:**

```json
{
  "instrument": {
    "type": "Orbitrap",
    "radius_in_m": 0.010,     // Inner spindle electrode
    "radius_out_m": 0.020,    // Outer barrel electrode  
    "radius_char_m": 0.015,   // Characteristic radius (defines field shape)
    "length_m": 0.100,        // Axial extent
    "dc_radial_V": 5000.0,    // Radial DC voltage
    "voltage_sweep": {        // Optional: ramped excitation
      "enabled": false,
      "start_time_s": 0.0005,
      "rise_time_s": 0.0005,
      "slope_V_s": 1000000.0
    }
  }
}
```

**Physics Notes:**

- The Orbitrap uses hyperboloidal electrodes to create a quadro-logarithmic potential
- Ions oscillate axially with frequency proportional to √(m/z)
- Radial motion provides orbital stability
- Field equation: E = k/2 [z·ẑ + (r/2)(1 - R_m²/r²)·r̂]

#### `tof_basic.json` – Time-of-Flight Mass Spectrometer

- **Purpose**: Basic TOF mass spectrometry
- **Ion Count**: 1,200 ions (multi-species)
- **Physics**: Field-free drift, extraction pulse
- **Runtime**: ~8 seconds

```bash
./build/icarion_core examples/tof_basic.json
```

#### `quad_basic.json` – Quadrupole Mass Filter

- **Purpose**: RF quadrupole mass filtering
- **Ion Count**: 1,800 ions (Ar+, N2+)
- **Physics**: RF quadrupole confinement, mass-selective transmission
- **Runtime**: ~5 seconds

```bash
./build/icarion_core examples/quad_basic.json
```

#### `fticr_basic.json` – FT-ICR Mass Spectrometer

- **Purpose**: Ion cyclotron resonance in strong magnetic field
- **Ion Count**: Variable
- **Physics**: Boris integrator, magnetic confinement
- **Runtime**: ~20 seconds

```bash
./build/icarion_core examples/fticr_basic.json
```

#### `spacecharge_basic.json` – Space Charge Test Configuration

- **Purpose**: Basic space charge dynamics
- **Ion Count**: 5,000 ions
- **Physics**: Coulomb interactions, self-consistent fields
- **Runtime**: ~12 seconds

```bash
./build/icarion_core examples/spacecharge_basic.json
```

### Advanced Physics Examples

#### `collision_models_comparison.json` – Collision Model Comparison

- **Purpose**: Compare EHSS, HSMC, and Langevin collision models
- **Ion Count**: 2,500 ions (H+, H2O+, N2+)
- **Physics**: Multiple collision models, drift field
- **Applications**: Method validation, cross-section measurements
- **Runtime**: ~15 seconds

```bash
./build/icarion_core examples/collision_models_comparison.json
```

#### `chemical_reactions.json` – Ion-Molecule Reactions

- **Purpose**: Chemical reaction network simulation
- **Ion Count**: 1,700 ions (H3O+, NH4+, CO2+)
- **Physics**: Langevin collisions, reaction kinetics
- **Applications**: Atmospheric chemistry, reaction rates
- **Runtime**: ~20 seconds

```bash
./build/icarion_core examples/chemical_reactions.json
```

### Performance Examples

#### `performance_benchmark.json` – Performance Validation

- **Purpose**: Benchmark ion dynamics performance
- **Ion Count**: 10,000+ ions
- **Physics**: IMS with EHSS collisions
- **Features**: Performance metrics, throughput validation
- **Runtime**: ~45 seconds

```bash
./build/icarion_core examples/performance_benchmark.json
```

#### `gpu_massive_ensemble.json` – GPU Acceleration

- **Purpose**: Large-scale GPU-accelerated simulation
- **Ion Count**: 50,000+ ions (Ar+ species)
- **Physics**: GPU-optimized integration, EHSS collisions
- **Requirements**: CUDA-capable GPU (NVIDIA)
- **Runtime**: ~30 seconds (GPU), ~180 seconds (CPU fallback)

```bash
./build/icarion_core --gpu examples/gpu_massive_ensemble.json
```

## Schema Compliance

All examples are **v1.0 schema compliant** and validated:

```bash
# Validate all examples
for f in examples/*.json; do
    python schema/validate_config.py "$f"
done
```

### Schema Features

- **Required fields**: `simulation`, `instrument`, `ions`, `fields`, `collisions`, `output`
- **Instrument types**: `IMS`, `LQIT`, `TOF`, `Quadrupole`, `FT_ICR`, `Orbitrap`, `spacecharge`
- **Collision models**: `EHSS`, `HSMC`, `Langevin`, `Friction`, `NoCollisions`
- **Integrators**: `rk4`, `rk45`, `boris`, `verlet`

The ICARION Core engine accepts case-insensitive values (e.g., `"ims"` works),  
but the schema validator requires canonical capitalization.

### Orbitrap-Specific Parameters

The Orbitrap instrument type requires unique geometry parameters:

```json
{
  "instrument": {
    "type": "Orbitrap",
    "radius_in_m": 0.010,      // Inner electrode radius [m]
    "radius_out_m": 0.020,     // Outer electrode radius [m]
    "radius_char_m": 0.015,    // Characteristic radius [m]
    "length_m": 0.100,         // Axial length [m]
    "dc_radial_V": 5000.0,     // Radial DC voltage [V]
    "voltage_sweep": {          // Optional voltage ramping
      "enabled": false,
      "start_time_s": 0.0,
      "rise_time_s": 0.001,
      "slope_V_s": 1000000.0
    }
  }
}
```

**Key differences from other instruments:**

- Uses hyperbolic boundary surfaces instead of simple cylinders
- Three radii define electrode geometry (inner, outer, characteristic)
- Voltage sweep enables resonant excitation experiments
- Field shape creates mass-dependent axial oscillation frequencies

## Output Files

Each simulation generates:

- **HDF5 trajectory file**: Ion positions, velocities, times
- **Metadata**: Simulation parameters, version info, random seeds
- **Diagnostics** (if enabled): Performance metrics, statistics

### Example Analysis (Python)

```python
import h5py
import numpy as np
import matplotlib.pyplot as plt

# Load trajectory data
with h5py.File('ims_basic.h5', 'r') as f:
    positions = f['trajectories/positions'][:]  # Shape: (steps, ions, 3)
    times = f['trajectories/times'][:]
    
    # Plot mean trajectory
    mean_x = np.mean(positions[:, :, 0], axis=1)
    plt.plot(times * 1e9, mean_x * 1000)  # Convert to ns, mm
    plt.xlabel('Time (ns)')
    plt.ylabel('Mean X Position (mm)')
    plt.title('IMS Drift Dynamics')
    plt.show()
```

See `docs/OUTPUT_SCHEMA.md` for complete HDF5 structure.

## Performance Expectations

Typical single-threaded CPU performance (Intel Xeon, `-O3` optimization):

| Example | Ion Count | Sim Time | Wall Time | ns/ion-step |
| --------- | ----------- | ---------- | ----------- | ------------- |
| `ims_basic.json` | 1,000¹ | 5 μs | ~15 s | ~30 ns |
| `lqit_basic.json` | 2,000 | Variable | ~10 s | ~25 ns |
| `quad_basic.json` | 1,800 | 100 ns | ~5 s | ~15 ns |
| `performance_benchmark.json` | 10,000 | 1 ms | ~45 s | ~45 ns |
| `gpu_massive_ensemble.json` (GPU) | 1,000¹ | 2.5 μs | ~30 s | ~5 ns |

¹ Engine may apply safety limits or normalization to configured ion counts.

**GPU Performance** (RTX 3060+): ~2-5 ns/ion-step for large ensembles.

Run the built-in performance benchmark to test your system:

```bash
./build/icarion_core --performance-test
./build/icarion_core --performance-test --gpu
```

Or use the detailed example configuration:

```bash
./build/icarion_core examples/performance_benchmark.json
```

## System Requirements

### Minimum

- **CPU**: x86-64 with SSE4.2
- **RAM**: 4 GB
- **Storage**: 1 GB free
- **OS**: Linux, macOS, Windows (WSL)

### Recommended

- **CPU**: Multi-core with AVX2 (Intel Xeon, AMD Ryzen/EPYC)
- **RAM**: 16 GB+
- **Storage**: SSD
- **GPU** (optional): NVIDIA with CUDA 11.0+ (4+ GB VRAM)

## Troubleshooting

### Schema Validation Warnings

**Problem**: `python schema/validate_config.py` shows warnings but example runs fine.

**Solution**: The C++ engine is more permissive (case-insensitive). Use:

```bash
# Engine validation (more permissive)
./build/icarion_core --validate examples/example.json

# Schema validation (strict canonical format)
python schema/validate_config.py examples/example.json
```

### Performance Issues

- Verify `-O3` optimization: `./build/icarion_core --dump-build-info | grep flags`
- Check CPU features: `cat /proc/cpuinfo | grep flags | head -1`
- Reduce `output_interval` in config to decrease I/O overhead
- Enable GPU: `--gpu` flag for supported examples

### Memory Issues

- Reduce ion count in configuration
- Increase `output_interval` (fewer snapshots stored)
- Enable compression: `"compression": true` in config

### GPU Examples Fail

- Check CUDA: `nvidia-smi`
- Verify GPU build: `./build/icarion_core --dump-build-info | grep cuda`
- Ensure 4+ GB GPU memory available

## Contributing Examples

To add new examples:

1. Follow v1.0 modular schemas: `schema/species.schema.json`, `schema/reactions.schema.json`, `schema/ions.schema.json`
2. Validate: `python schema/validate_config.py your_example.json`
3. Test execution: `./build/icarion_core your_example.json`
4. Document physics, runtime, applications

## Documentation

- **Schemas**: `schema/species.schema.json`, `schema/reactions.schema.json`, `schema/ions.schema.json`
- **CLI Reference**: `docs/CLI_INTERFACE_v1.0.md`
- **Input Format**: `docs/INPUT_SCHEMA.md`
- **Output Format**: `docs/OUTPUT_SCHEMA.md`
- **C++ API**: `docs/PUBLIC_CPP_API_v1.0.md`

## License

MIT License. Free use for academic and commercial research.

## Log Analysis

ICARION v1.0++ supports JSON-formatted logs for automated analysis.

### Generate JSON Logs

```bash
# Run with JSON logging
./build/icarion_main --log-format json --log-file simulation.log examples/ims_basic.json

# Or redirect console output
./build/icarion_main --log-format json examples/ims_basic.json 2>&1 | grep '^{' > simulation.log
```

### Analyze Logs with Python

**Basic analysis:**

```bash
python examples/analyze_logs.py simulation.log
```

Output:

- Summary statistics (total entries, time span)
- Log level distribution (info, warn, error)
- Category breakdown (main, config, hdf5, perf, etc.)
- Error and warning messages
- Export to CSV for further analysis

**Advanced analysis:**

```bash
python examples/analyze_logs_advanced.py simulation.log
```

Features:

- Filter by category and level
- Time-based analysis (find slowest operations)
- Pattern matching (extract species, domains, particle counts)
- Export to CSV, Excel, JSON

**Requirements:**

```bash
pip install pandas
pip install openpyxl  # Optional: for Excel export
```

### Example: Query Logs with pandas

```python
import pandas as pd
import json

# Load JSON logs
logs = []
with open('simulation.log') as f:
    for line in f:
        if line.strip().startswith('{'):
            logs.append(json.loads(line))

df = pd.DataFrame(logs)
df['time'] = pd.to_datetime(df['time'])

# Filter by category
config_logs = df[df['cat'] == 'config']
print(config_logs[['time', 'msg']])

# Find errors
errors = df[df['level'] == 'error']
print(f"Found {len(errors)} errors")

# Performance analysis
perf_logs = df[df['cat'] == 'perf']
print(f"Performance measurements: {len(perf_logs)}")
```

## Support

- **Documentation**: Main `README.md` and `docs/` directory
- **Performance Testing**: Use `examples/performance_benchmark.json` for system benchmarks
- **Issues**: Repository issue tracker
- **Development**: `docs/DEVELOPERS_GUIDE.md`

---

**ICARION Core v1.0** – High-Performance Ion Collision and Reaction Integration Framework

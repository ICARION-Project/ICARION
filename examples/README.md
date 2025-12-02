# ICARION Examples

Example configurations demonstrating mass spectrometry instruments and simulation capabilities.

## Quick Start

```bash
./build/src/icarion_main examples/<folder>/<config>.json
```

Results written to `results/`. GPU-accelerated when `enable_gpu: true` and N ≥ 5000 ions (Boris threshold ~2500).

## Instrument Examples

### Linear Quadrupole Ion Trap (LQIT) - `lqit/`
- **lqit_basic.json**: RF trapping at q=0.40 (620V @ 1MHz), N₂ buffer gas (3.75 mTorr), ReserpineH+

### Orbitrap - `orbitrap/`
- **orbitrap_basic.json**: Hyperlogarithmic field, f ∝ 1/√m, R > 100k

### Quadrupole Mass Filter - `quadrupole/`
- **quadrupole_basic.json**: RF-only mode (DC=0V), ion transmission without mass filtering

### Time-of-Flight - `tof/`
- **tof_basic.json**: Linear TOF, 1m flight tube, t ∝ √m
- **tof_reflectron.json**: Orthogonal reflectron, energy focusing

### FTICR - `fticr/`
- **fticr_basic.json**: 7T B-field, f_c = qB/(2πm), Penning trap

### Ion Mobility - `ims/`
- **ims_basic.json**: Drift tube, 5000 H₃O⁺ ions, HSS collisions
- **ims_with_field_array.json**: External field (400 V/cm, scaled by DC.axial_V)
- **ims_field_array_time_varying.json**: Superposition of DC (500 V/cm) + RF (±100V @ 1MHz)
- **ims_field_array_multi_domain.json**: Two-stage drift (600 V/cm → 300 V/cm)

### Reactions - `reactions/`
- **reaction_demo.json**: First-order decay kinetics

## Configuration Format

Basic structure (see `schema/` for full specification):

```json
{
  "simulation": {
    "total_time_s": 1e-4,
    "dt_s": 1e-9,
    "enable_gpu": true,
    "integrator": "RK4"
  },
  "physics": {
    "collision_model": "HSS",
    "enable_space_charge": false
  },
  "ions": {
    "species": [
      {
        "id": "H3O+",
        "count": 5000,
        "position": {"type": "gaussian", "center": [0,0,0], "std": [0.001,0.001,0.001]},
        "velocity": {"type": "thermal", "temperature_K": 300}
      }
    ]
  },
  "domains": [
    {
      "name": "drift",
      "instrument": "IMS",
      "geometry": {"origin_m": [0,0,0], "length_m": 0.05, "radius_m": 0.01},
      "environment": {"temperature_K": 300, "pressure_Pa": 200, "gas_species": "He"},
      "fields": {"dc": {"EN_Td": 10.0}}
    }
  ]
}
```

## Output

- **HDF5**: `results/<name>/<name>_trajectories.h5` (positions, velocities, species, timestamps)
- **Log**: `results/<name>/simulation.log` (if enabled via config)

## Documentation

- Configuration: `docs/CONFIG_GUIDE.md`
- Collision models: `docs/COLLISION_MODELS.md`
- HDF5 output: `docs/HDF5_OUTPUT_STRUCTURE.md`
- Schema: `schema/README.md`

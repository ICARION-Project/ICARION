# ICARION Examples

Example configurations demonstrating mass spectrometry instruments and simulation capabilities.

## Quick Start

```bash
./build/src/icarion_main examples/<folder>/<config>.json
```

Results written to `results/`. GPU runtime path is disabled in v1.0.0; `enable_gpu=true` falls back to CPU.

## Instrument Examples

### Linear Quadrupole Ion Trap (LQIT) - `lqit/`
- **lqit_basic.json**: RF trapping at q~=0.40 (620 V @ 1 MHz), N2 buffer gas (0.5 Pa / 3.75 mTorr), ReserpineH+ (includes a small axial DC term).

### Orbitrap - `orbitrap/`
- **orbitrap_basic.json**: Hyperlogarithmic field, collisionless, two species (ReserpineH+, CaffeineH+), radial_V=3500 V.

### Quadrupole Mass Filter - `quadrupole/`
- **quadrupole_basic.json**: RF quadrupole field (DC quad_V=0, axial_V=500 V), CaffeineH+ beam transport.

### Time-of-Flight - `tof/`
- **tof_basic.json**: Linear TOF, 1 m flight tube, 2 kV acceleration, two species.
- **tof_reflectron.json**: Reflectron mirror with radial offset, energy focusing.

### FTICR - `fticr/`
- **fticr_basic.json**: 7 T B-field, quadrupolar trapping, collisionless ions.

### Ion Mobility - `ims/`
- **ims_basic.json**: Drift tube, 500 H3O+ + 500 26DTBPH+ ions, HSS collisions.
- **ims_with_field_array.json**: Field array example (200 V/cm base array scaled by DC.axial_V to 400 V/cm).
- **ims_field_array_time_varying.json**: Superposition of DC (500 V/cm) + RF (+/-100 V @ 1 MHz)
- **ims_field_array_multi_domain.json**: Two-stage drift (600 V/cm to 300 V/cm)
- **ims_ehss_offline_basic.json**: EHSS run using a precomputed single gas offline sample file.
- **ims_ipm_basic.json**: InteractionPotentialModel setup using bundled example `ipm_offline_samples` tables in `data/molecules/precomputed_ipm/`.

### Reactions - `reactions/`
- **reaction_demo.json**: Ion-molecule reaction demo (H3O+ + Pentanal -> PentanalH+) in a gas mixture.

## Configuration Format

Basic structure (see `schema/` for full specification):

```json
{
  "simulation": {
    "total_time_s": 1e-4,
    "dt_s": 1e-9,
    "enable_gpu": false,
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
      "env": {"temperature_K": 300, "pressure_Pa": 200, "gas_species": "He"},
      "fields": {"DC": {"EN_Td": 10.0}}
    }
  ]
}
```

## Output

- **HDF5**: `output.folder/output.trajectory_file` (positions, velocities, species, timestamps)
- **Log**: `output.folder/simulation.log` (if enabled by the run)

## Documentation

- Configuration: `docs/CONFIG_GUIDE.md`
- Collision models: `docs/COLLISION_MODELS.md`
- HDF5 output: `docs/HDF5_OUTPUT_STRUCTURE.md`
- Schema: `schema/README.md`

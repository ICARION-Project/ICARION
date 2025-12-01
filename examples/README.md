# ICARION Examples

Ready-to-run example configurations demonstrating different instrument types and simulation capabilities.

## Quick Start

Run any example from the ICARION root directory:

```bash
cd /path/to/ICARION
./build/src/icarion_main examples/<instrument>/<example_name>.json
```

**Example:**
```bash
./build/src/icarion_main examples/lqit/lqit_basic.json
```

Results are written to `results/<example_name>/`.

## Examples by Instrument Type

### 📊 Linear Quadrupole Ion Trap (LQIT)

**Directory:** `lqit/`

Demonstrates RF ion trapping in a linear Paul trap.

- **lqit_basic.json**: Stable trapping at q=0.4 with HSS collisions

**Key Physics:**
- Mathieu stability parameters (a, q)
- RF confinement with secular motion
- Collision cooling

[**→ Full documentation**](lqit/README.md)

---

### 🔄 Orbitrap

**Directory:** `orbitrap/`

Demonstrates frequency-based mass measurement in hyperlogarithmic fields.

- **orbitrap_basic.json**: Axial oscillations for multiple species

**Key Physics:**
- Mass-dependent axial frequencies (f ∝ 1/√m)
- Fourier transform mass spectrometry
- High-resolution mass measurement

[**→ Full documentation**](orbitrap/README.md)

---

### 📐 Quadrupole Mass Filter

**Directory:** `quadrupole/`

Demonstrates mass-selective ion transmission.

- **quadrupole_basic.json**: Stable transmission at q=0.85, a=-0.15

**Key Physics:**
- Mathieu stability diagram
- Mass filtering via RF/DC scanning
- Transmission efficiency

[**→ Full documentation**](quadrupole/README.md)

---

### ⚗️ Reactions

**Directory:** `reactions/`

Demonstrates ion-molecule reactions and chemical kinetics.

- **reaction_demo.json**: First-order decay kinetics

**Key Physics:**
- First-order and bimolecular reactions
- Rate constants and collision theory
- Species transformation

[**→ Full documentation**](reactions/README.md)

---

## Additional Examples

### Ion Mobility Spectrometry (IMS)

- **ims_basic.json**: Simple IMS simulation with H₃O⁺ ions
- **ims_with_species_db.json**: IMS using species database
- **ims_with_field_array.json**: IMS with externally defined electric field
- **ims_basic_ion_cloud.json**: IMS with realistic ion cloud initialization

### Collision Models

- **collision_models_comparison.json**: Compare NoCollisions, HSS, EHSS

### Field Arrays

- **field_array_validation.json**: Validate field array implementation
- **field_array_time_varying.json**: Time-varying electric fields
- **field_array_multi_domain.json**: Multiple field domains

### RF and Waveforms

- **test_rf_modulation.json**: RF amplitude modulation
- **test_rf_superposition.json**: Superposition of multiple RF fields

### Multi-Gas Environments

- **multi_gas_air.json**: Simulation with air composition (N₂, O₂, Ar, CO₂)

### Advanced Features

- **test_rk45_adaptive.json**: Adaptive timestep integration (RK45)
- **profile.json**: Configuration for performance profiling

## GPU Examples

If built with CUDA support:

- **test_gpu_boris.json**: Boris integrator on GPU
- **test_gpu_rk45.json**: RK45 integrator on GPU
- **test_gpu_space_charge_force.json**: Space charge calculations on GPU
- **test_gpu_boundary_absorption.json**, **test_gpu_boundary_reflection.json**: GPU boundary handling

## Configuration Structure

All examples follow the JSON schema in `schema/icarion-config.schema.json`:

```json
{
  "simulation": { "total_time_s": 1e-4, "dt_s": 1e-9 },
  "ions": { "species": [...], "count": 100 },
  "environment": { "temperature_K": 300, "gases": [...] },
  "fields": { "electric": [...], "magnetic": [...] },
  "geometry": { "boundaries": [...] },
  "output": { "format": "hdf5", "write_interval": 1000 }
}
```

## Data Resources

- **field_arrays/**: Pre-generated field data (HDF5 format)
- **ion_clouds/**: Pre-defined initial ion distributions
- **waveforms/**: Custom RF/DC waveform definitions

## Output Files

Each simulation produces:

1. **Trajectories**: `results/<name>/<name>_trajectories.h5` (HDF5)
   - Ion positions, velocities, timestamps
   - Species IDs, active flags
   
2. **Metadata**: `results/<name>/<name>_metadata.json`
   - Configuration parameters
   - Performance metrics
   - Simulation statistics

## Analysis Scripts

Python examples for data analysis:

```bash
python examples/analyze_logs.py results/<name>/<name>_metadata.json
python examples/analyze_logs_advanced.py results/<name>/
```

## Documentation

- **Configuration**: `docs/CONFIG_GUIDE.md`
- **CLI Usage**: `docs/CLI_USAGE.md`
- **Collision Models**: `docs/COLLISION_MODELS.md`
- **Field Arrays**: `examples/FIELD_ARRAYS_README.md`
- **JSON Schema**: `schema/README.md`

## Performance Tips

- **Thread control**: Use `--threads N` to set OpenMP threads
- **GPU acceleration**: Build with `-DUSE_CUDA=ON` for GPU support
- **Memory usage**: ~100 MB per 10k ions for 1M timesteps
- **Timestep selection**: See `docs/PERFORMANCE.md` (coming soon)

## Contributing

To add a new example:
1. Create config file in appropriate subdirectory
2. Add description to subdirectory README.md
3. Test with `./build/src/icarion_main examples/<dir>/<config>.json`
4. Commit with descriptive message

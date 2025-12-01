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
  - 100 ReserpineH+ ions in N₂ buffer gas
  - RF: 1 MHz at 500V, produces q=0.4 stability
  - Collision cooling over 100 µs
  - Secular motion frequencies ~50-100 kHz

**Key Physics:**
- Mathieu stability parameters (a, q)
- RF confinement with secular motion
- Hard-sphere scattering (HSS) collision model
- Damping and focusing towards trap center

---

### 🔄 Orbitrap

**Directory:** `orbitrap/`

Demonstrates frequency-based mass measurement in hyperlogarithmic fields.

- **orbitrap_basic.json**: Axial oscillations for multiple species
  - ReserpineH+ (m/z 609) and CaffeineH+ (m/z 195)
  - Hyperlogarithmic field: E_z ∝ ln(r)
  - Axial frequencies: ~180 kHz (Caffeine), ~115 kHz (Reserpine)
  - 100 µs simulation captures multiple oscillation periods

**Key Physics:**
- Mass-dependent axial frequencies: f ∝ 1/√m
- Fourier transform gives mass spectrum
- High-resolution mass measurement (R > 100,000)
- Independent of initial kinetic energy

---

### 📐 Quadrupole Mass Filter

**Directory:** `quadrupole/`

Demonstrates mass-selective ion transmission.

- **quadrupole_basic.json**: Stable transmission at q=0.85, a=-0.15
  - 50 ions each of m/z 195, 609, 1000
  - RF: 800 kHz at 400V, DC: -30V
  - Operating point in stability region (first stability island)
  - Unstable ions ejected, stable ions transmitted

**Key Physics:**
- Mathieu stability diagram (a-q space)
- Mass filtering via RF/DC ratio scanning
- Transmission efficiency depends on (a,q) position
- Resolution vs. transmission trade-off

---

### ⏱️ Time-of-Flight (TOF)

**Directory:** `tof/`

Demonstrates mass separation by flight time.

**tof_basic.json** - Linear TOF with orthogonal acceleration:
- **Acceleration**: 2 kV over 2 cm (100 kV/m field)
- **Flight tube**: 1 m field-free drift region
- **Ions**: Reserpine (m/z 609) and Caffeine (m/z 195)
- **Flight times**: 
  - Caffeine: ~23 µs
  - Reserpine: ~41 µs
  - Ratio: 1.78 ≈ √(609/195) = 1.77 ✓
- **Resolution**: Baseline separation by mass

**tof_reflectron.json** - Orthogonal reflectron for improved resolution:
- **Geometry**: Radially offset reflectron arm (x = 8 cm)
- **Ion trajectory**: Diagonal launch → deflection → reflection → detection
- **Reflectron**: -2.2 kV retarding field in 15 cm mirror
- **Flight times**:
  - Caffeine: ~11.3 µs
  - Reserpine: ~20.0 µs
  - Ratio: 1.77 (perfect mass separation maintained)
- **Total path**: ~1.4 m (includes radial deflection)

**Key Physics:**
- Flight time: t = L × √(m/2qV)
- All ions accelerated to same kinetic energy
- Velocity inversely proportional to √m
- Reflectron corrects initial energy spread (higher energy → penetrates deeper → longer path)
- Orthogonal geometry prevents direct line-of-sight

---

### 🔁 Fourier Transform Ion Cyclotron Resonance (FTICR)

**Directory:** `fticr/`

Demonstrates mass measurement via cyclotron frequency in strong magnetic field.

- **fticr_basic.json**: Penning trap with quadrupolar electric trapping
  - **B-field**: 7 Tesla (axial, defines cyclotron motion)
  - **Trapping**: 20V quadrupolar electric field for axial confinement
  - **Cell**: 5 cm length, 1.5 cm radius cylindrical trap
  - **Ions**: 10 ReserpineH+ (m/z 609) + 10 CaffeineH+ (m/z 195)
  - **Initial conditions**: Point source at center, 1K thermal velocities
  - **Cyclotron frequencies**:
    - Caffeine: f_c = 560 kHz
    - Reserpine: f_c = 179 kHz
  - **Survival times**: Heavy ions 60-80 µs, light ions 35-45 µs
  - **Simulation time**: 500 µs

**Key Physics:**
- Cyclotron frequency: f_c = qB/(2πm) - mass measurement via frequency
- Magnetic field creates circular motion in xy-plane
- Electric quadrupole provides harmonic axial (z) oscillation
- Ultra-high mass resolution: R > 1,000,000 possible
- Ultra-high vacuum required (10⁻⁹ Pa)

**Limitations in this example:**
- Ions escape after 35-80 µs due to thermal motion
- Even at 1K, thermal velocities can exceed shallow trapping potential
- Optimal trapping voltage ~20V (too low = weak confinement, too high = overfocusing)
- Real FTICR uses sophisticated trapping sequences, buffer gas cooling, and optimized cell designs
- This example demonstrates the cyclotron principle but not long-term stable trapping

---

### ⚗️ Reactions

**Directory:** `reactions/`

Demonstrates ion-molecule reactions and chemical kinetics.

- **reaction_demo.json**: First-order decay kinetics

**Key Physics:**
- First-order and bimolecular reactions
- Rate constants and collision theory
- Species transformation

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

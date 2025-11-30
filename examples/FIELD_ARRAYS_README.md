# Field Array Examples

This directory contains examples demonstrating field array loading and usage in ICARION.

## Overview

Field arrays allow using precomputed electric fields from BEM/FEM solvers instead of analytical field models. This enables realistic electrode geometries and complex field distributions.

## Examples

### 1. Basic Field Array Usage
**File:** `field_array_validation.json`

Simple single-domain simulation with constant field array scaling.

```bash
./build/src/icarion_main examples/field_array_validation.json
```

**Features:**
- Single domain (IMS drift tube)
- Uniform 200 V/m field from HDF5 file
- No collisions (ballistic motion)
- Validates interpolation accuracy

**Expected Result:** 
- Ions accelerate ~47 mm in 10 µs
- Accuracy: ~91% vs analytical prediction

---

### 2. Multi-Domain Field Arrays
**File:** `field_array_multi_domain.json`

Demonstrates independent field arrays in multiple domains.

```bash
./build/src/icarion_main examples/field_array_multi_domain.json
```

**Features:**
- Two domains with different field strengths
  - Domain 0: 200 V/m (scale=1.0)
  - Domain 1: 100 V/m (scale=0.5)
- Ions start in different domains
- Tests per-domain force registry

**Expected Result:**
- Domain 0 ions drift ~2× faster than Domain 1 ions
- Ratio: ~1.9:1

---

### 3. Time-Varying Fields (Future)
**File:** `field_array_time_varying.json`

**⚠️ NOT YET IMPLEMENTED** - Framework in place but not active.

This example shows the planned configuration format for:
- RF-modulated field arrays
- DC voltage scaling (axial/quad/radial)
- Field superposition

```json
"field_array_terms": [
  {
    "file": "field_array.h5",
    "scale_type": "DC_Axial",
    "_comment": "Scale by DC.axial_V voltage"
  },
  {
    "file": "rf_field.h5",
    "scale_type": "RF",
    "frequency_Hz": 1e6,
    "phase_rad": 0.0,
    "_comment": "RF modulation at 1 MHz"
  }
]
```

---

## Field Array Format

HDF5 files must contain:

### Required Datasets
- `x`, `y`, `z`: 1D coordinate arrays [m]
- `Ex`, `Ey`, `Ez`: 3D field components [V/m]

### Optional Datasets
- `phi`: Electric potential [V] (3D array)

### Example Structure
```
field_array.h5
├── x: (20,) float64   # X coordinates: [-0.01, ..., 0.01] m
├── y: (20,) float64   # Y coordinates: [-0.01, ..., 0.01] m
├── z: (50,) float64   # Z coordinates: [0.0, ..., 0.05] m
├── Ex: (20,20,50) float64  # X-component of E-field
├── Ey: (20,20,50) float64  # Y-component of E-field
└── Ez: (20,20,50) float64  # Z-component of E-field
```

---

## Creating Field Arrays

### Python Example
```python
import h5py
import numpy as np

# Create grid
x = np.linspace(-0.01, 0.01, 20)  # 10mm radius
y = np.linspace(-0.01, 0.01, 20)
z = np.linspace(0.0, 0.05, 50)     # 50mm length

# Create uniform field (200 V/m in z-direction)
Ex = np.zeros((len(x), len(y), len(z)))
Ey = np.zeros((len(x), len(y), len(z)))
Ez = np.full((len(x), len(y), len(z)), 200.0)

# Save to HDF5
with h5py.File('ims_drift_200Vcm.h5', 'w') as f:
    f.create_dataset('x', data=x)
    f.create_dataset('y', data=y)
    f.create_dataset('z', data=z)
    f.create_dataset('Ex', data=Ex)
    f.create_dataset('Ey', data=Ey)
    f.create_dataset('Ez', data=Ez)
```

---

## Configuration Reference

### Field Array Term Structure
```json
{
  "file": "path/to/field_array.h5",
  "scale_type": "Constant",
  "constant_V": 1.0,
  "frequency_Hz": 0.0,
  "phase_rad": 0.0
}
```

### Scale Types

| Type | Description | Status |
|------|-------------|--------|
| `Constant` | Static scaling by `constant_V` | ✅ Working |
| `DC_Axial` | Scale by DC.axial_V voltage | ⚠️ Parsed, not applied |
| `DC_Quad` | Scale by DC.quad_V voltage | ⚠️ Parsed, not applied |
| `DC_Radial` | Scale by DC.radial_V voltage | ⚠️ Parsed, not applied |
| `RF` | RF modulation: `V*cos(ωt+φ)` | ⚠️ Parsed, not applied |

---

## Implementation Status

### ✅ Working
- HDF5 field array loading
- Trilinear interpolation
- Multi-domain support
- GridFieldProvider integration
- ElectricFieldForce field provider mode

### ⚠️ Planned (Config Parsed, Not Active)
- Field superposition (multiple arrays per domain)
- Time-varying scaling (RF, DC voltage modulation)
- Waveform-driven field arrays

### 📋 Future
- Magnetic field arrays
- Field array caching/preloading
- Adaptive field refinement
- GPU texture upload for field arrays (infrastructure exists)

---

## Validation

### Accuracy Test
```bash
# Run validation
./build/src/icarion_main examples/field_array_validation.json

# Analyze results
python3 << EOF
import h5py, numpy as np
with h5py.File('results/field_array_validation/trajectories.h5') as f:
    z = f['trajectory/positions'][:,:,2]
    print(f"Mean displacement: {(z[-1]-z[0]).mean()*1000:.2f} mm")
EOF
```

**Expected:** ~47 mm (ballistic acceleration in 200 V/m field)

### Multi-Domain Test
```bash
./build/src/icarion_main examples/field_array_multi_domain.json
# Should see 2:1 drift ratio between domains
```

---

## Troubleshooting

### Issue: "Unable to open file"
**Cause:** Output directory doesn't exist  
**Solution:**
```bash
mkdir -p results/field_array_validation
```

### Issue: Ions not accelerating
**Possible causes:**
1. Field array not loaded (check log for "Loading field arrays")
2. Ions outside field grid bounds
3. Collision model too aggressive (try `"collision_model": "None"`)

### Issue: Field interpolation returning zero
**Debug:**
```bash
# Check field array contents
python3 -c "import h5py; f=h5py.File('field.h5'); print(f['Ez'][:].mean())"

# Check ion positions vs field bounds
# Ions must be inside [x.min, x.max] × [y.min, y.max] × [z.min, z.max]
```

---

## See Also
- `docs/HDF5_OUTPUT_STRUCTURE.md`: HDF5 file formats
- `src/core/io/fieldArrayLoader.h`: Field loading implementation
- `src/fieldsolver/utils/GridFieldProvider.h`: Interpolation logic

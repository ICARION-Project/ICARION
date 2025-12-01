# ICARION v1.0 Known Limitations

**Version:** 1.0.0  
**Last Updated:** December 1, 2025

This document outlines the known limitations, constraints, and planned improvements for ICARION v1.0.

---

## Table of Contents

- [Physics Limitations](#physics-limitations)
- [Numerical Limitations](#numerical-limitations)
- [Performance Limitations](#performance-limitations)
- [Feature Limitations](#feature-limitations)
- [v1.1 Roadmap](#v11-roadmap)

---

## Physics Limitations

### 1. Space Charge Solver

**Limitation:**
- Direct N-body solver is O(N²), practical only for N < 10,000 ions
- Grid-based solver available but with fixed grid resolution
- No adaptive mesh refinement (AMR)

**Workaround:**
```json
{
  "physics": {
    "space_charge": {
      "method": "Grid",  // Use for N > 10,000
      "grid_points": [64, 64, 64]
    }
  }
}
```

**Planned (v1.1):** Tree-based fast multipole method (FMM) for O(N log N) scaling

---

### 2. Collision Models

**Limitation:**
- EHSS (Enhanced Hard-Sphere Statistical) requires detailed molecular geometry
- Not all species have validated CCS data in database
- Gas mixtures supported but collision cross-section averaging may be approximate

**Workaround:**
- Use HSS (Hard-Sphere Statistical) for species without detailed geometry
- Add custom species to `data/species_database_v1.json`

**Example:**
```json
{
  "species_id": "MyCustomIon",
  "mass_amu": 100.0,
  "charge": 1,
  "CCS_m2": 150e-20  // Measured or estimated
}
```

---

### 3. Reactions

**Supported:**
- ✅ First-order decay (A → B)
- ✅ Bimolecular reactions (A + B → C + D)

**Not Supported (v1.0):**
- ❌ Three-body reactions
- ❌ Dissociative electron attachment (DEA)
- ❌ Photoionization
- ❌ Charge transfer with energy-dependent cross-sections

**Planned (v1.1):** Extended reaction types with energy-dependent rate constants

---

### 4. Electric Fields

**Supported:**
- ✅ Uniform DC fields
- ✅ Cylindrical gradient fields
- ✅ Quadrupole (Paul trap) RF fields
- ✅ Time-varying waveforms
- ✅ Field arrays (imported from external solvers)

**Limitation:**
- No automatic field solver included (use external FEM/FDM tools)
- Field arrays use trilinear interpolation only (no higher-order)

**Workaround:** Use COMSOL, SIMION, or Bempp to generate field maps, then import:
```json
{
  "fields": {
    "type": "FieldArray",
    "source": "path/to/field_map.h5"
  }
}
```

---

## Numerical Limitations

### 1. Time Integration

**Limitation:**
- RK4 uses fixed timestep (no adaptivity)
- RK45 labeled "adaptive" but currently fixed-step in v1.0
- No implicit solvers for stiff systems

**Impact:**
- User must manually choose appropriate timestep
- Very stiff systems may require impractically small timesteps

**Rule of Thumb:**
```
dt < 0.1 × min(T_collision, T_oscillation, T_cyclotron)
```

**Planned (v1.1):** True adaptive RK45 with error control

---

### 2. Boundary Conditions

**Supported:**
- ✅ Dirichlet (absorbing walls)
- ✅ Cylindrical/spherical boundaries with accurate distance calculation

**Not Supported:**
- ❌ Neumann (reflecting walls)
- ❌ Robin (partially absorbing)
- ❌ Periodic boundary conditions

**Impact:** Ions are absorbed upon hitting boundaries (cannot simulate reflecting electrodes)

**Planned (v1.1):** Specular and diffuse reflection boundary conditions

---

### 3. Field Interpolation

**Limitation:**
- Only trilinear interpolation for field arrays
- No cubic or spline interpolation

**Impact:**
- Field gradients may be inaccurate near steep features
- Higher-order methods needed for precise trajectory near electrodes

**Workaround:** Use finer field array resolution

**Planned (v1.2):** Cubic spline interpolation option

---

## Performance Limitations

### 1. CPU Parallelization

**Observed Scaling:**
| Threads | Speedup | Efficiency |
|---------|---------|------------|
| 1       | 1.0×    | 100%       |
| 2       | 1.89×   | 94%        |
| 4       | 3.54×   | 88%        |
| 8       | 4.74×   | 59%        |
| 16      | 5.19×   | 32%        |

**Limitation:**
- Memory bandwidth saturates at 4-8 cores (typical DDR4 system)
- Beyond 8 cores: diminishing returns

**Recommendation:** Use 4-8 threads for best efficiency

---

### 2. GPU Acceleration

**Status (v1.0):**
- ✅ GPU integrators (RK4, Boris) functional
- ✅ GPU collision handlers (HSS, Langevin) functional
- ⚠️ GPU space charge experimental (P³M solver)
- ❌ No multi-GPU support (single GPU only)

**Limitation:**
- GPU memory limits max ions (typically ~1M ions on 8GB GPU)
- Data transfer overhead for small simulations (N < 10,000)

**Recommendation:**
- Use GPU for N > 10,000 ions
- Use CPU for smaller systems or when space charge is critical

**Planned (v1.1):**
- Multi-GPU support with domain decomposition
- Unified memory for larger simulations

---

### 3. Memory Usage

**Limitation:**
- All trajectory data stored in RAM (no streaming to disk)
- Memory scales linearly with (N_ions × N_steps × output_interval)

**Estimate:**
```
Memory (GB) ≈ N_ions × N_steps / write_interval × 100 bytes
```

**Example:**
- 100,000 ions, 1M steps, write_interval=1000
- Memory ≈ 100,000 × 1,000 × 100 bytes ≈ 10 GB

**Workaround:** Increase `write_interval` to reduce output frequency

**Planned (v1.1):** Streaming output mode to reduce memory footprint

---

## Feature Limitations

### 1. Output Format

**Supported:**
- ✅ HDF5 trajectory files

**Not Supported:**
- ❌ VTK format
- ❌ CSV export
- ❌ Real-time visualization

**Workaround:** Use Python/MATLAB to convert HDF5 to desired format:
```python
import h5py
import pandas as pd

with h5py.File('trajectories.h5', 'r') as f:
    df = pd.DataFrame({
        'x': f['trajectory/x'][:].flatten(),
        'y': f['trajectory/y'][:].flatten(),
        'z': f['trajectory/z'][:].flatten()
    })
    df.to_csv('trajectories.csv')
```

---

### 2. Input Flexibility

**Limitation:**
- Configuration must be JSON (no YAML, TOML, or XML)
- No GUI or interactive mode

**Planned (v1.2):** Python API for programmatic configuration

---

### 3. Logging and Diagnostics

**Limitation:**
- No built-in energy/momentum conservation checks
- No automatic timestep stability warnings

**Workaround:** Use validation scripts to post-process results:
```bash
python validation/scripts/check_energy_conservation.py results/trajectories.h5
```

**Planned (v1.1):** Optional runtime diagnostics for conservation laws

---

## v1.1 Roadmap

**Planned Improvements:**

### High Priority
- ✅ True adaptive RK45 with error control
- ✅ Multi-GPU support (domain decomposition)
- ✅ Streaming output mode (reduce memory)
- ✅ Neumann/Robin boundary conditions
- ✅ FMM for O(N log N) space charge

### Medium Priority
- 🔄 Extended reaction types (three-body, DEA)
- 🔄 Cubic spline field interpolation
- 🔄 Energy/momentum conservation checks
- 🔄 Python API (programmatic interface)
- 🔄 YAML/TOML configuration support

### Low Priority
- 🔜 VTK output format
- 🔜 Real-time visualization (experimental)
- 🔜 GUI for configuration (community-driven)
- 🔜 MPI for multi-node clusters

---

## Reporting Issues

If you encounter a limitation not documented here:

1. **Check documentation:** [`docs/`](../docs/)
2. **Search issues:** https://github.com/YOUR_ORG/ICARION/issues
3. **Report new limitation:** Use issue template "Limitation/Feature Request"

---

## Contributing

Contributions to address these limitations are welcome!

See [`CONTRIBUTING.md`](../CONTRIBUTING.md) for guidelines.

---

**Last Updated:** December 1, 2025  
**ICARION Version:** 1.0.0

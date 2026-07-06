# Physics Validation Tests

This directory contains **high-fidelity scientific validation scripts** for physics modules.

These are **NOT CTests** - they are long-running (up to 30 minutes), high-accuracy validations for scientific publications.

## Current Validations

### Gas Flow Transport (`validate_gas_flow_transport.py`)
- **Runtime**: ~5-10 minutes
- **Physics**: Ion transport by gas flow without electric field (SIFT-MS)
- **Tests**: 
  - Pressure dependence (100 Pa, 1000 Pa, 5000 Pa)
  - Terminal velocity: v_terminal = v_gas
  - Thermalization time: τ ∝ 1/P
- **Ensemble**: 1000 ions per condition
- **Output**: 
  - `validation/figures/physics/gas_flow_transport_validation.png`
  - `validation/logs/GAS_FLOW_TRANSPORT_VALIDATION.txt`

### Combined Drift (`validate_combined_drift.py`)
- **Runtime**: ~5 minutes
- **Physics**: Superposition of electric mobility and gas flow (IMS baseline)
- **Tests** (all capped at ≤10 Td):
  - `E=0`, `P=100 Pa` (gas-flow baseline)
  - `E=0`, `P=1000 Pa` (gas-flow baseline)
  - `E=1000 V/m`, `P=1000 Pa` (μ·E + gas flow)
- **Ensemble**: 1000 ions per condition, 1 µs runtime
- **Output**:
  - `validation/figures/physics/combined_drift_validation.png`
  - `validation/logs/COMBINED_DRIFT_VALIDATION.txt`

### TIMS Elution (`validate_tims_elution.py`)
- **Runtime**: ~15-20 seconds
- **Physics**: Mobility-sorted TIMS elution under axial gas flow and a ramped axial field
- **Tests**:
  - Three-species elution order by reduced mobility
  - Peak separation in elution time
  - Approximate agreement with the analytical ramp-fraction estimate
- **Ensemble**: 60 ions per species, 3 species
- **Output**:
  - `validation/figures/physics/tims_elution_validation.png`
  - `validation/logs/TIMS_ELUTION_VALIDATION.txt`

### Adaptive SC Parity (`validate_space_charge_adaptive_parity.py`)
- **Runtime**: ~1 minute
- **Physics**: Stage-synchronous adaptive RK45 space charge vs. fixed-step RK4 (Direct SC)
- **Tests**:
  - Small cloud (200 ions), no external fields
  - Reference RK4 (dt=5e-11 s, t=5e-8 s) vs. adaptive RK45 (dt=5e-9 s)
- **Output**:
  - `validation/logs/ADAPTIVE_SC_PARITY.txt`

### Adaptive SC Microbench (`bench_space_charge_adaptive.py`)
- **Runtime**: ~2–3 minutes
- **Physics**: SC rebuild/rejection overhead for adaptive RK45 (CPU Direct/Grid)
- **Tests**:
  - N = 1k, 5k, 10k ions; no external fields
  - dt = 5e-9 s, t = 5e-7 s
- **Output**:
  - `validation/logs/ADAPTIVE_SC_BENCH.txt`

## Running Validations

From repository root:

```bash
# Single validation
python validation/scripts/physics/validate_gas_flow_transport.py
python validation/scripts/physics/validate_tims_elution.py

# All physics validations (when more are added)
# python validation/scripts/physics/run_all_physics_validations.py
```

## Planned Validations

- [ ] Collision cross-section validation (HSS vs EHSS)
- [ ] Reaction rate validation (H3O+ + VOC → products)
- [ ] Space charge effects on transport
- [ ] Multi-gas mixture thermalization

## Output Structure

```
validation/
├── figures/
│   └── physics/
│       ├── gas_flow_transport_validation.png
│       └── combined_drift_validation.png
├── logs/                       # Detailed validation logs
│   ├── GAS_FLOW_TRANSPORT_VALIDATION.txt
│   └── COMBINED_DRIFT_VALIDATION.txt
└── results/                    # HDF5 trajectories
    ├── gas_flow_transport/
    │   ├── config_100Pa.json
    │   ├── gas_flow_100Pa.h5
    │   ├── config_1000Pa.json
    │   ├── gas_flow_1000Pa.h5
    │   └── ...
    └── combined_drift/
        ├── config_E0_100Pa.json
        ├── combined_E0_100Pa.h5
        ├── config_E5000_1000Pa.json
        ├── combined_E5000_1000Pa.h5
        └── ...
```

## Adding New Validations

1. Create script: `validate_<physics_module>.py`
2. Follow template structure:
   - Large ensemble (≥1000 ions)
   - Long simulation time (sufficient for convergence)
   - Multiple parameter sweeps
   - Quantitative comparison to theory
   - Publication-quality plots
3. Document in this README
4. Add to `run_all_physics_validations.py` (when created)

## Quality Standards

All physics validations must:
- ✅ Run for sufficient time to reach statistical convergence
- ✅ Use large ensembles (≥1000 ions for statistics)
- ✅ Compare quantitatively to theoretical predictions
- ✅ Generate publication-quality figures
- ✅ Document all assumptions and limitations
- ✅ Include error analysis and tolerance checks

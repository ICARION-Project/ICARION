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

## Running Validations

From repository root:

```bash
# Single validation
python validation/scripts/physics/validate_gas_flow_transport.py

# All physics validations (when more are added)
# python validation/scripts/physics/run_all_physics_validations.py
```

## Planned Validations

- [ ] Ion transport with E-field + gas flow (combined drift)
- [ ] Collision cross-section validation (HSS vs EHSS)
- [ ] Reaction rate validation (H3O+ + VOC → products)
- [ ] Space charge effects on transport
- [ ] Multi-gas mixture thermalization

## Output Structure

```
validation/
├── figures/
│   └── physics/
│       └── gas_flow_transport_validation.png
├── logs/                       # Detailed validation logs
│   └── GAS_FLOW_TRANSPORT_VALIDATION.txt
└── results/                    # HDF5 trajectories
    └── gas_flow_transport/
        ├── config_100Pa.json
        ├── gas_flow_100Pa.h5
        ├── config_1000Pa.json
        ├── gas_flow_1000Pa.h5
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

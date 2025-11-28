# IMS Drift Velocity Validation

## Test Matrix

**Objective:** Validate Mason-Schamp mobility equation across 3 collision models

```
v_drift = K₀ × E × (N₀/N)
```

### Parameters:
- **Ion:** H3O+ (K₀ = 3.2 cm²/(V·s) in N2)
- **E-fields:** 1000, 5000, 10000 V/m
- **Pressures:** 100, 1000, 10000 Pa
- **Collision Models:**
  - HSS (Hard-Sphere Scattering) - Stochastic
  - EHSS (Elastic HSS) - Stochastic with molecular geometry
  - Friction - Deterministic mobility-based

### Total: 27 configurations

## Expected Results

Lower pressure → Higher drift velocity (fewer collisions)
Higher E-field → Higher drift velocity (stronger acceleration)

### Model-Specific Tolerances:
- **HSS**: ±25% (spherical approximation)
- **EHSS**: ±50% (molecular geometry effects)
- **Friction**: ±10% (mobility-based, should be most accurate)

## Usage

```bash
# Run single test
cd /home/chsch95/ICARION/validation
../build/src/icarion_main configs/instruments/ims/ims_friction_5000Vm_1000Pa.json

# Run all IMS tests
./scripts/run_instrument_tests.sh ims
```

## Analysis

Extract drift velocities and compare to Mason-Schamp predictions:

```bash
python3 scripts/analyze_ims.py results/instruments/ims/
```

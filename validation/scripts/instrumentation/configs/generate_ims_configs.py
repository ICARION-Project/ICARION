#!/usr/bin/env python3
"""
Generate IMS drift velocity validation configurations

Tests 5 collision models against Mason-Schamp mobility equation:
  v_drift = K₀ × E × (N₀/N)

Collision Models:
  1. HSS (Hard-Sphere Scattering) - Stochastic, full randomization
  2. EHSS (Elastic Hard-Sphere Scattering) - Stochastic, molecular geometry
  3. Friction - Deterministic, mobility-based damping

Test Matrix:
  - 3 E-field strengths (1000, 5000, 10000 V/m)
  - 3 pressures (100, 1000, 10000 Pa)
  - 1 ion species (H3O+)
  - 3 collision models
  Total: 3×3×3 = 27 configurations
"""

import json
import os
from pathlib import Path


def find_validation_dir() -> Path:
    """Locate the validation directory regardless of invocation cwd."""
    for parent in Path(__file__).resolve().parents:
        if parent.name == "validation":
            return parent
    raise RuntimeError("Could not locate 'validation' directory relative to script.")


VALIDATION_DIR = find_validation_dir()
CONFIG_ROOT = VALIDATION_DIR / "configs"

# IMS parameters
DRIFT_LENGTH_M = 0.06  # 5 cm drift tube
DRIFT_TUBE_RADIUS_M = 0.05  # 5 cm radius (large to minimize wall losses)

# Parameter sweeps using EN_Td (reduced field) - detailed E/N mapping
EN_TD_VALUES = [1.0, 2.0, 3.0, 5.0, 7.0, 10.0]  # Td - detailed E/N mapping for validation
PRESSURES_PA = [100.0, 200.0, 500.0, 1000.0, 2000.0, 5000.0]  # Pa - including 20 Pa for better resolution
ION_SPECIES = "H3O+"

# Collision models depend on pressure regime:
# - Very low pressure (10 Pa = 0.1 mbar): Only stochastic models - individual collisions
# - Low/medium pressure (100-1000 Pa): All models (stochastic + Friction for comparison)
# - High pressure (10000 Pa = 100 mbar): Only Friction model is physically meaningful
#   (stochastic models become computationally intractable and less accurate)
def get_models_for_pressure(pressure_Pa):
    """Return appropriate collision models for given pressure"""
    if pressure_Pa >= 10000:  # High pressure regime (>50 mbar)
        return ["Friction"]  # Only deterministic damping
    elif pressure_Pa <= 20:  # Very low pressure (molecular flow regime)
        return ["HSS", "EHSS"]  # Only stochastic - Friction not valid
    else:  # Low/medium pressure (transition regime)
        return ["HSS", "EHSS", "Friction"]  # All models for comparison

# Physical constants
T_K = 300.0  # Temperature
K0_H3Op_SI = 24.1e-4  # m²/(V·s) - H3O+ in He reduced mobility (much higher than N2)
P0 = 101325.0  # Pa
T0 = 273.15  # K

# Simulation parameters
N_IONS = 1000  # Many ions for good statistics and validation
RNG_SEED = 42

def calc_EN_to_E(EN_Td, pressure_Pa):
    """Convert reduced field E/N [Td] to E-field [V/m]"""
    k_B = 1.380649e-23  # J/K
    N = pressure_Pa / (k_B * T_K)  # Number density [m^-3]
    E_Vm = EN_Td * 1e-21 * N  # 1 Td = 1e-21 V·m²
    return E_Vm

def calc_drift_time(EN_Td, pressure_Pa):
    """Calculate expected drift time from Mason-Schamp equation"""
    E_Vm = calc_EN_to_E(EN_Td, pressure_Pa)
    N_ratio = (P0 / pressure_Pa) * (T_K / T0)
    v_drift = K0_H3Op_SI * E_Vm * N_ratio
    t_drift = DRIFT_LENGTH_M / v_drift
    return t_drift, v_drift, E_Vm

def calc_dt(collision_freq_Hz):
    """Calculate timestep from collision frequency
    
    Target: P_collision = ν * dt < 0.1 (10% collision probability per step)
    Therefore: dt < 0.1 / ν
    
    Use 5% collision probability for safety: dt = 0.05 / ν
    """
    target_prob = 0.05  # 5% collision probability per timestep
    dt = target_prob / collision_freq_Hz
    return dt

def calc_collision_freq(pressure_Pa):
    """Estimate collision frequency from kinetic theory"""
    # ν = n * σ * v_rel
    # Simplified: ν ≈ P / (k_B * T) * σ * v_thermal
    k_B = 1.380649e-23  # J/K
    m_gas = 4.0 * 1.66054e-27  # He mass in kg (much lighter than N2)
    m_ion = 19.0 * 1.66054e-27  # H3O+ mass
    mu = (m_ion * m_gas) / (m_ion + m_gas)  # Reduced mass
    
    v_rel = ((8 * k_B * T_K) / (3.14159 * mu))**0.5  # Mean relative velocity
    n = pressure_Pa / (k_B * T_K)  # Number density
    sigma = 24.9e-20  # CCS in m² for H3O+ in He (much smaller than N2)
    
    nu = n * sigma * v_rel
    return nu

def generate_config(EN_Td, pressure_Pa, collision_model):
    """Generate a single IMS validation config"""
    
    t_drift, v_drift, E_Vm = calc_drift_time(EN_Td, pressure_Pa)
    collision_freq = calc_collision_freq(pressure_Pa)
    
    # Simulation time: 2x drift time is enough for drift velocity measurement
    # At high E/N, ions cross quickly - we don't need 3-5x
    total_time_s = 2.0 * t_drift
    
    # Timestep: 5% collision probability per step (P = ν * dt < 0.05)
    dt_s = calc_dt(collision_freq)
    
    # Write interval: Target ~200 snapshots - enough for analysis, reasonable file size
    total_steps = int(total_time_s / dt_s)
    write_interval = max(1, total_steps // 200)
    
    config = {
        "simulation": {
            "total_time_s": total_time_s,
            "dt_s": dt_s,
            "write_interval": write_interval,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "num_threads": 4,
            "rng_seed": RNG_SEED
        },
        "physics": {
            "collision_model": collision_model,
            "enable_space_charge": False,
            "enable_reactions": False
        },
        "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
        "output": {
            "folder": "results/v1.0_test/instruments/ims",
            "trajectory_file": f"ims_{collision_model.lower()}_{E_Vm:.0f}Vm_{pressure_Pa:.0f}Pa.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": ION_SPECIES,
                    "count": N_IONS,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.0],
                        "std": [0.0005, 0.0005, 0.0001]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": T_K
                    }
                }
            ]
        },
        "domains": [
            {
                "domain_index": 0,
                "name": "IMS Drift Tube",
                "instrument": "IMS",
                "geometry": {
                    "origin_m": [0.0, 0.0, -0.01],
                    "length_m": DRIFT_LENGTH_M,
                    "radius_m": DRIFT_TUBE_RADIUS_M
                },
                "env": {
                    "gas_species": "He",
                    "temperature_K": T_K,
                    "pressure_Pa": pressure_Pa,
                    "gas_velocity_m_s": [0.0, 0.0, 0.0]
                },
                "fields": {
                    "DC": {
                        "EN_Td": EN_Td
                    }
                }
            }
        ]
    }
    
    return config

def main():
    """Generate all IMS validation configs"""
    
    output_dir = CONFIG_ROOT / "instruments" / "ims"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Generating IMS validation configs in {output_dir}")
    print(f"Test matrix: Pressure-dependent model selection (stochastic models only at low/medium pressure)")
    
    count = 0
    for EN_Td in EN_TD_VALUES:
        for pressure_Pa in PRESSURES_PA:
            models = get_models_for_pressure(pressure_Pa)
            for collision_model in models:
                
                config = generate_config(EN_Td, pressure_Pa, collision_model)
                
                # Filename
                filename = f"ims_{collision_model.lower()}_{EN_Td:.0f}Td_{pressure_Pa:.0f}Pa.json"
                filepath = output_dir / filename
                
                # Write config
                with open(filepath, 'w') as f:
                    json.dump(config, f, indent=2)
                
                count += 1
                
                # Print summary
                t_drift, v_drift, E_Vm = calc_drift_time(EN_Td, pressure_Pa)
                collision_freq = calc_collision_freq(pressure_Pa)
                print(f"  {count:2d}. {collision_model:9s} | E/N={EN_Td:4.0f} Td ({E_Vm:5.0f} V/m) | P={pressure_Pa:5.0f} Pa | "
                      f"v_drift={v_drift:6.1f} m/s | t_drift={t_drift*1e6:5.1f} µs | ν_coll={collision_freq/1e6:.2f} MHz")
    
    print(f"\n✅ Generated {count} IMS validation configurations")
    
    # Create README
    readme_path = output_dir / "README.md"
    with open(readme_path, 'w') as f:
        f.write("""# IMS Drift Velocity Validation

## Test Matrix

**Objective:** Validate Mason-Schamp mobility equation across 5 collision models

```
v_drift = K₀ × E × (N₀/N)
```

### Parameters:
- **Ion:** H3O+ (K₀ = 10.5 cm²/(V·s) in He)
- **Reduced fields (E/N):** 10, 40, 100 Td (realistic IMS range)
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
- **Langevin**: ±20% (polarization corrections)
- **HSD**: ±25% (deterministic hard-sphere)

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
""")
    
    print(f"✅ Created README: {readme_path}")

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Generate transport physics validation suite: Drift Velocity

Validates drift velocity vs E/N relationship:
  v_drift = K₀ × E
  
where K₀ is the reduced mobility (literature values available).

Systematic parameter sweep:
- E/N (reduced field): 10-200 Td
- Temperature: 150K, 300K, 500K
- Collision models: HSS, Langevin, Friction
- Test species: H3O+, PentanalH+, ReserpineH+
"""

import json
import numpy as np
from pathlib import Path


def find_validation_dir() -> Path:
    for parent in Path(__file__).resolve().parents:
        if parent.name == "validation":
            return parent
    raise RuntimeError("Unable to locate 'validation' directory relative to script")


VALIDATION_DIR = find_validation_dir()
OUT_DIR = VALIDATION_DIR / "configs" / "physics" / "transport"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Test parameters
TEMPERATURES_K = [150.0, 300.0, 500.0]
EN_VALUES_TD = [10, 30, 50, 100, 150, 200]  # Townsend units (1 Td = 10⁻¹⁷ V·cm²)
COLLISION_MODELS = ["HSS", "Langevin", "Friction"]

# Test species (different masses and CCS)
TEST_SPECIES = [
    {
        "name": "H3O+",
        "mass_amu": 19.0,
        "ccs_m2": 3.2e-19,  # ~32 Ų
        "K0_literature_cm2_Vs": 2.8,  # Literature reduced mobility in He
    },
    {
        "name": "PentanalH+",
        "mass_amu": 87.0,
        "ccs_m2": 6.5e-19,  # ~65 Ų
        "K0_literature_cm2_Vs": 1.9,
    },
    {
        "name": "ReserpineH+",
        "mass_amu": 609.0,
        "ccs_m2": 16.5e-19,  # ~165 Ų
        "K0_literature_cm2_Vs": 0.85,
    },
]

# Buffer gas (Helium)
BUFFER_GAS = "He"
BUFFER_GAS_MASS_AMU = 4.0

# Drift tube geometry
DRIFT_LENGTH_M = 0.050  # 50 mm
DRIFT_RADIUS_M = 0.010  # 10 mm

# Simulation parameters
SIM_TIME_FACTOR = 5.0  # Simulate 5× drift time
N_IONS = 200  # Statistical ensemble
DT_FACTOR = 100  # Timestep: τ_collision / 100
WRITE_INTERVAL = 100

print("\n" + "="*80)
print("Transport Physics: Drift Velocity Validation Suite Generator")
print("="*80)
print(f"Temperatures: {TEMPERATURES_K} K")
print(f"E/N values: {EN_VALUES_TD} Td")
print(f"Collision models: {COLLISION_MODELS}")
print(f"Test species: {[s['name'] for s in TEST_SPECIES]}")
print(f"Buffer gas: {BUFFER_GAS}")
print(f"Drift length: {DRIFT_LENGTH_M*1000:.1f} mm")
print("="*80 + "\n")

configs_generated = 0

for temp_K in TEMPERATURES_K:
    for en_td in EN_VALUES_TD:
        for collision_model in COLLISION_MODELS:
            for species in TEST_SPECIES:
                
                # Calculate pressure from E/N
                # E/N [Td] = E/(n·10⁻¹⁷) where n is number density [m⁻³]
                # Ideal gas: n = P/(k_B·T)
                # Therefore: P = E/(E/N)·k_B·T / 10⁻¹⁷
                
                # Choose electric field (V/m) to give reasonable drift time
                # Target: ~100 µs drift time
                # v_drift ≈ K₀ × E ≈ 1-3 cm²/(V·s) × E
                # For 50 mm drift in 100 µs: v = 500 m/s
                # E ≈ 500 / (K₀ in m²/(V·s)) ≈ 500 / (K₀·10⁻⁴) = 5e6 / K₀
                
                K0_estimate = species["K0_literature_cm2_Vs"]  # cm²/(V·s)
                E_field = 5e6 / K0_estimate  # V/m (target ~100 µs drift)
                
                # Calculate pressure from E/N
                k_B = 1.380649e-23  # J/K
                en_si = en_td * 1e-21  # Convert Td to V·m²
                pressure_Pa = E_field / en_si * k_B * temp_K
                
                # Calculate number density
                n_density = pressure_Pa / (k_B * temp_K)  # m⁻³
                
                # Estimate collision time for timestep
                # τ_coll ≈ 1 / (n·σ·v_rel)
                # v_rel ≈ √(8kT/(πμ)) where μ is reduced mass
                mass_ion_kg = species["mass_amu"] * 1.66054e-27
                mass_gas_kg = BUFFER_GAS_MASS_AMU * 1.66054e-27
                mu = (mass_ion_kg * mass_gas_kg) / (mass_ion_kg + mass_gas_kg)
                v_rel = np.sqrt(8 * k_B * temp_K / (np.pi * mu))
                
                sigma = species["ccs_m2"]
                tau_coll = 1.0 / (n_density * sigma * v_rel)
                
                # Timestep
                dt_s = tau_coll / DT_FACTOR
                
                # Drift time estimate
                v_drift_est = K0_estimate * 1e-4 * E_field  # m/s
                t_drift_est = DRIFT_LENGTH_M / v_drift_est
                
                # Simulation time
                sim_time_s = SIM_TIME_FACTOR * t_drift_est
                
                # Configuration
                config = {
                    "simulation": {
                        "total_time_s": sim_time_s,
                        "dt_s": dt_s,
                        "write_interval": WRITE_INTERVAL,
                        "integrator": "RK4",
                        "enable_gpu": False,
                        "enable_openmp": True,
                        "rng_seed": 42
                    },
                    "physics": {
                        "collision_model": collision_model,
                        "enable_space_charge": False,
                        "enable_reactions": False
                    },
                    "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
                    "output": {
                        "folder": "validation/results/v1.0_test/physics/transport/drift",
                        "trajectory_file": f"drift_{collision_model}_{species['name']}_{temp_K:.0f}K_{en_td}Td.h5",
                        "print_progress": True
                    },
                    "ions": {
                        "species": [{
                            "species_id": species["name"],
                            "count": N_IONS,
                            "position": {
                                "type": "gaussian",
                                "center": [0.0, 0.0, -DRIFT_LENGTH_M/2],  # Start at beginning
                                "std": [0.001, 0.001, 0.001]  # Tight initial cloud
                            },
                            "velocity": {
                                "type": "thermal",
                                "temperature_K": temp_K
                            }
                        }]
                    },
                    "domains": [{
                        "domain_index": 0,
                        "name": f"Drift Tube {collision_model} {species['name']} T={temp_K}K E/N={en_td}Td",
                        "instrument": "IMS",
                        "geometry": {
                            "type": "cylinder",
                            "radius_m": DRIFT_RADIUS_M,
                            "length_m": DRIFT_LENGTH_M,
                            "origin_m": [0.0, 0.0, -DRIFT_LENGTH_M/2]
                        },
                        "env": {
                            "gas_species": BUFFER_GAS,
                            "temperature_K": temp_K,
                            "pressure_Pa": pressure_Pa,
                            "gas_velocity_m_s": [0.0, 0.0, 0.0]
                        },
                        "fields": {
                            "DC": {
                                "uniform_field_V_m": [0.0, 0.0, E_field]  # Axial E-field
                            }
                        }
                    }]
                }
                
                # Write config
                filename = f"drift_{collision_model}_{species['name']}_{temp_K:.0f}K_{en_td}Td.json"
                out_file = OUT_DIR / filename
                with open(out_file, 'w') as f:
                    json.dump(config, f, indent=2)
                
                configs_generated += 1
                
                # Print summary for first few configs
                if configs_generated <= 3:
                    print(f"Generated: {filename}")
                    print(f"  E-field: {E_field:.1f} V/m")
                    print(f"  Pressure: {pressure_Pa:.2e} Pa")
                    print(f"  E/N: {en_td:.1f} Td")
                    print(f"  τ_coll: {tau_coll*1e6:.2f} µs")
                    print(f"  dt: {dt_s*1e9:.2f} ns")
                    print(f"  v_drift (est): {v_drift_est:.1f} m/s")
                    print(f"  t_drift (est): {t_drift_est*1e6:.1f} µs")
                    print(f"  Sim time: {sim_time_s*1e6:.0f} µs")
                    print()

print("="*80)
print(f"Total configurations: {configs_generated}")
print(f"Output: {OUT_DIR}")
print("\nBreakdown:")
print(f"  {len(TEMPERATURES_K)} temperatures × {len(EN_VALUES_TD)} E/N × {len(COLLISION_MODELS)} models × {len(TEST_SPECIES)} species")
print(f"  = {len(TEMPERATURES_K) * len(EN_VALUES_TD) * len(COLLISION_MODELS) * len(TEST_SPECIES)} configs")
print("\nNext steps:")
print("  1. Run: for cfg in configs/physics/transport/drift_*.json; do")
print("             ../build/src/icarion_main $cfg; done")
print("  2. Analyze: python3 scripts/analyze_transport_drift.py")
print("  3. Validate K₀ extraction vs literature values")
print("="*80 + "\n")

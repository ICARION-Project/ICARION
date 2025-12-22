#!/usr/bin/env python3
"""
Generate reaction kinetics validation suite

Validates:
1. First-order decay: A+ → B+ (exponential decay)
2. Bimolecular: A+ + N2 → B+ + N2 (rate equation)

Tests reaction rate constant extraction and species evolution.
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
OUT_DIR = VALIDATION_DIR / "configs" / "physics" / "reactions"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Environment (typical IMS conditions)
TEMPERATURE_K = 300.0
PRESSURE_PA = 1000.0  # 1 kPa

# Reaction rate constants (realistic values)
# First-order: k in s^-1
# Bimolecular: k in cm³/s (will convert to SI)

FIRST_ORDER_RATES_S = [1e3, 5e3, 1e4]  # 1k, 5k, 10k /s (fast reactions)
BIMOLECULAR_RATES_CM3_S = [1e-10, 5e-10, 1e-9]  # Typical ion-molecule rates

# Simulation parameters
N_IONS = 1000  # Moderate ensemble for statistics
RNG_SEED = 42

# Domain (small box, reactions in gas phase)
DOMAIN_SIZE_M = 0.01  # 1 cm cube (no fields, pure kinetics)

print("\n" + "="*80)
print("Reaction Kinetics Validation Suite Generator")
print("="*80)
print(f"Temperature: {TEMPERATURE_K:.1f} K")
print(f"Pressure: {PRESSURE_PA:.1f} Pa")
print(f"Domain: {DOMAIN_SIZE_M*100:.1f} cm cube (field-free)")
print("="*80 + "\n")

configs_generated = 0

# ============================================================================
# FIRST-ORDER DECAY: A+ → B+ (unimolecular, no reagent gas)
# ============================================================================

print("Generating First-Order Decay Configurations")
print("-" * 60)

for k_s in FIRST_ORDER_RATES_S:
    # Time constant: τ = 1/k
    tau_s = 1.0 / k_s
    
    # Simulate for 5 time constants (> 99% conversion)
    total_time_s = 5.0 * tau_s
    
    # Timestep: τ/100 (fine resolution)
    dt_s = tau_s / 100.0
    
    # Write interval: every τ/10 (10 points per time constant)
    write_interval = max(1, int((tau_s / 10.0) / dt_s))
    
    config = {
        "simulation": {
            "total_time_s": total_time_s,
            "dt_s": dt_s,
            "write_interval": write_interval,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": RNG_SEED
        },
        "physics": {
            "collision_model": "None",  # No collisions, pure reactions
            "enable_space_charge": False,
            "enable_reactions": True,
            "reactions": {
                "database_path": "/home/chsch95/ICARION/data/reaction_database.json",
                "reactions": [
                    {
                        "reactant_ion": "ReactantA+",
                        "product_ion": "ProductB+",
                        "rate_constant_cm3_s": k_s * 1e6 / 2.5e19,  # Convert to pseudo-bimolecular
                        "reaction_type": "unimolecular"
                    }
                ]
            }
        },
        "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
        "output": {
            "folder": "validation/results/v1.0_test/physics/reactions",
            "trajectory_file": f"first_order_decay_k{k_s:.0e}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": "ReactantA+",
                    "count": N_IONS,
                    "position": {
                        "type": "uniform",
                        "min": [0.0, 0.0, 0.0],
                        "max": [DOMAIN_SIZE_M, DOMAIN_SIZE_M, DOMAIN_SIZE_M]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": TEMPERATURE_K
                    }
                }
            ]
        },
        "domains": [{
            "domain_index": 0,
            "name": f"First-order decay k={k_s:.0e} /s",
            "instrument": "Generic",
            "geometry": {
                "type": "box",
                "dimensions": [DOMAIN_SIZE_M, DOMAIN_SIZE_M, DOMAIN_SIZE_M],
                "origin_m": [0.0, 0.0, 0.0]
            },
            "env": {
                "gas_species": "He",
                "temperature_K": TEMPERATURE_K,
                "pressure_Pa": PRESSURE_PA,
                "gas_velocity_m_s": [0.0, 0.0, 0.0]
            },
            "fields": {}  # No fields
        }]
    }
    
    # Write config
    out_file = OUT_DIR / f"first_order_decay_k{k_s:.0e}.json"
    with open(out_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    configs_generated += 1
    
    print(f"✅ k = {k_s:.1e} /s")
    print(f"   τ = {tau_s*1e6:.1f} µs")
    print(f"   Simulation: {total_time_s*1e6:.1f} µs (5τ)")
    print(f"   Output: {out_file.name}")

print(f"\nGenerated {len(FIRST_ORDER_RATES_S)} first-order configs")

# ============================================================================
# BIMOLECULAR: A+ + N2 → B+ + N2
# ============================================================================

print("\n" + "="*80)
print("Generating Bimolecular Reaction Configurations")
print("-" * 60)

for k_cm3_s in BIMOLECULAR_RATES_CM3_S:
    # Reaction rate: dn/dt = -k * n * [N2]
    # Number density N2 at pressure P and temperature T
    k_B = 1.380649e-23  # J/K
    n_N2 = PRESSURE_PA / (k_B * TEMPERATURE_K)  # m^-3
    
    # Pseudo-first-order rate constant
    k_eff_s = k_cm3_s * 1e-6 * n_N2  # Convert cm³/s to m³/s, multiply by [N2]
    
    # Time constant
    tau_s = 1.0 / k_eff_s
    
    # Simulate for 5 time constants
    total_time_s = min(5.0 * tau_s, 1e-3)  # Cap at 1 ms
    
    # Timestep: τ/100
    dt_s = tau_s / 100.0
    
    # Write interval
    write_interval = max(1, int((tau_s / 10.0) / dt_s))
    
    config = {
        "simulation": {
            "total_time_s": total_time_s,
            "dt_s": dt_s,
            "write_interval": write_interval,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": RNG_SEED
        },
        "physics": {
            "collision_model": "HSS",  # Enable collisions for realistic scenario
            "enable_space_charge": False,
            "enable_reactions": True,
            "reactions": {
                "database_path": "/home/chsch95/ICARION/data/reaction_database.json",
                "reactions": [
                    {
                        "reactant_ion": "ReactantA+",
                        "reagent_gas": "N2",
                        "product_ion": "ProductB+",
                        "rate_constant_cm3_s": k_cm3_s,
                        "reaction_type": "bimolecular"
                    }
                ]
            }
        },
        "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
        "output": {
            "folder": "validation/results/v1.0_test/physics/reactions",
            "trajectory_file": f"bimolecular_k{k_cm3_s:.0e}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": "ReactantA+",
                    "count": N_IONS,
                    "position": {
                        "type": "uniform",
                        "min": [0.0, 0.0, 0.0],
                        "max": [DOMAIN_SIZE_M, DOMAIN_SIZE_M, DOMAIN_SIZE_M]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": TEMPERATURE_K
                    }
                }
            ]
        },
        "domains": [{
            "domain_index": 0,
            "name": f"Bimolecular A+ + N2 → B+ (k={k_cm3_s:.0e} cm³/s)",
            "instrument": "Generic",
            "geometry": {
                "type": "box",
                "dimensions": [DOMAIN_SIZE_M, DOMAIN_SIZE_M, DOMAIN_SIZE_M],
                "origin_m": [0.0, 0.0, 0.0]
            },
            "env": {
                "gas_species": "N2",
                "temperature_K": TEMPERATURE_K,
                "pressure_Pa": PRESSURE_PA,
                "gas_velocity_m_s": [0.0, 0.0, 0.0]
            },
            "fields": {}
        }]
    }
    
    # Write config
    out_file = OUT_DIR / f"bimolecular_k{k_cm3_s:.0e}.json"
    with open(out_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    configs_generated += 1
    
    print(f"✅ k = {k_cm3_s:.1e} cm³/s")
    print(f"   [N2] = {n_N2:.2e} m⁻³")
    print(f"   k_eff = {k_eff_s:.2e} /s")
    print(f"   τ = {tau_s*1e6:.1f} µs")
    print(f"   Simulation: {total_time_s*1e6:.1f} µs")
    print(f"   Output: {out_file.name}")

print(f"\nGenerated {len(BIMOLECULAR_RATES_CM3_S)} bimolecular configs")

print("\n" + "="*80)
print(f"Total configurations: {configs_generated}")
print(f"Output: {OUT_DIR}")
print("\nExpected behavior:")
print("  First-order:  N_A(t) = N_0 * exp(-k*t)")
print("  Bimolecular:  N_A(t) = N_0 * exp(-k*[N2]*t)")
print("\nNext steps:")
print("  1. Run: for cfg in configs/physics/reactions/*.json; do")
print("             icarion_main $cfg; done")
print("  2. Analyze: python3 analyze_reactions.py")
print("  3. Validate rate constants from exponential fits")
print("="*80 + "\n")

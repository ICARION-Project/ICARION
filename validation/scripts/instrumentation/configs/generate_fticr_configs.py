#!/usr/bin/env python3
"""
Generate FT-ICR validation suite

Validates cyclotron frequency vs m/z relationship:
  f_c = (q*B) / (2π*m)

Tests multiple species to verify exact 1/m scaling.
"""

import json
import numpy as np
from pathlib import Path

# Output directory
OUT_DIR = Path("../configs/instruments/fticr")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# FT-ICR geometry (cylindrical Penning trap)
GEOMETRY = {
    "radius_m": 0.025,      # 25 mm trap radius
    "length_m": 0.100,      # 100 mm trap length
    "origin_m": [0.0, 0.0, -0.050]  # Center trap at z=0
}

# Environment (UHV)
ENVIRONMENT = {
    "gas_species": "He",
    "temperature_K": 300.0,
    "pressure_Pa": 1e-9,  # Ultra-high vacuum
    "gas_velocity_m_s": [0.0, 0.0, 0.0]
}

# Magnetic field (Tesla)
# Typical FT-ICR: 7-15 T (high-field systems)
# We'll use moderate field for validation
B_FIELD_T = 7.0  # 7 Tesla

# Quadrupole trapping voltage (radial confinement)
# V_trap creates radial restoring force
TRAP_VOLTAGE_V = 5.0  # 5V quadrupole trapping

# Test species (wide mass range)
TEST_SPECIES = [
    {"name": "H3O+", "mass_amu": 19.0, "count": 50},
    {"name": "PentanalH+", "mass_amu": 87.0, "count": 50},
    {"name": "CaffeineH+", "mass_amu": 195.0, "count": 50},
    {"name": "ReserpineH+", "mass_amu": 609.0, "count": 50},
]

# Simulation parameters
SIM_TIME_S = 200e-6  # 200 µs (multiple cyclotron orbits)
DT_S = 1e-10  # 100 ps timestep (need fine resolution for fast cyclotron motion)
WRITE_INTERVAL = 1000  # Write every 100 ns

print("\n" + "="*80)
print("FT-ICR Validation Suite Generator")
print("="*80)
print(f"Magnetic field: {B_FIELD_T:.1f} T")
print(f"Trap voltage: {TRAP_VOLTAGE_V:.1f} V")
print(f"Geometry: r={GEOMETRY['radius_m']*1000:.1f} mm, L={GEOMETRY['length_m']*1000:.1f} mm")
print(f"Test species: {len(TEST_SPECIES)}")
print(f"Simulation: {SIM_TIME_S*1e6:.0f} µs")
print("="*80 + "\n")

# Calculate expected cyclotron frequencies
print("Expected cyclotron frequencies:")
print("-" * 60)
q = 1.602176634e-19  # Elementary charge
for species in TEST_SPECIES:
    m_kg = species["mass_amu"] * 1.66054e-27
    f_c_theory = (q * B_FIELD_T) / (2 * np.pi * m_kg)
    print(f"{species['name']:15s} (m={species['mass_amu']:6.1f} u): "
          f"f_c = {f_c_theory/1e6:.3f} MHz")
print("-" * 60 + "\n")

configs_generated = 0

# Generate multi-species configuration
config = {
    "simulation": {
        "total_time_s": SIM_TIME_S,
        "dt_s": DT_S,
        "write_interval": WRITE_INTERVAL,
        "integrator": "RK4",
        "enable_gpu": False,
        "enable_openmp": True,
        "rng_seed": 42
    },
    "physics": {
        "collision_model": "HSS",
        "enable_space_charge": False,
        "enable_reactions": False
    },
    "species_database_path": "/home/chsch95/ICARION/data/species_database_v1.json",
    "output": {
        "folder": "validation/results/v1.0_test/instruments/fticr",
        "trajectory_file": f"fticr_multi_species_B{B_FIELD_T:.1f}T.h5",
        "print_progress": True
    },
    "ions": {
        "species": []
    },
    "domains": [{
        "domain_index": 0,
        "name": f"FT-ICR B={B_FIELD_T:.1f}T V={TRAP_VOLTAGE_V:.1f}V",
        "instrument": "FTICR",
        "geometry": GEOMETRY,
        "env": ENVIRONMENT,
        "fields": {
            "magnetic": {
                "enabled": True,
                "field_strength_T": [0.0, 0.0, B_FIELD_T]
            },
            "DC": {
                "radial_V": TRAP_VOLTAGE_V  # Quadrupole trapping
            }
        }
    }]
}

# Add all test species
for species in TEST_SPECIES:
    config["ions"]["species"].append({
        "species_id": species["name"],
        "count": species["count"],
        "position": {
            "type": "gaussian",
            "center": [0.0, 0.0, 0.0],  # Center of trap
            "std": [0.002, 0.002, 0.005]  # Small spread
        },
        "velocity": {
            "type": "thermal",
            "temperature_K": 300.0
        }
    })

# Write multi-species config
out_file = OUT_DIR / f"fticr_multi_species_B{B_FIELD_T:.1f}T.json"
with open(out_file, 'w') as f:
    json.dump(config, f, indent=2)

configs_generated += 1
print(f"✅ Generated multi-species configuration:")
print(f"   {out_file}")
print(f"   Total ions: {sum(s['count'] for s in TEST_SPECIES)}")
print(f"   Species: {', '.join(s['name'] for s in TEST_SPECIES)}")

# Also generate individual single-species configs
for species in TEST_SPECIES:
    config_single = config.copy()
    config_single["output"]["trajectory_file"] = f"fticr_{species['name']}_B{B_FIELD_T:.1f}T.h5"
    config_single["ions"] = {
        "species": [{
            "species_id": species["name"],
            "count": 100,  # More ions for single-species
            "position": {
                "type": "gaussian",
                "center": [0.0, 0.0, 0.0],
                "std": [0.002, 0.002, 0.005]
            },
            "velocity": {
                "type": "thermal",
                "temperature_K": 300.0
            }
        }]
    }
    
    out_file = OUT_DIR / f"fticr_{species['name']}_B{B_FIELD_T:.1f}T.json"
    with open(out_file, 'w') as f:
        json.dump(config_single, f, indent=2)
    
    configs_generated += 1

print(f"\n✅ Generated {len(TEST_SPECIES)} single-species configurations")
print(f"\n{'='*80}")
print(f"Total configurations: {configs_generated}")
print(f"Output: {OUT_DIR}")
print("\nNext steps:")
print("  1. Run: for cfg in configs/instruments/fticr/*.json; do")
print("             icarion_main $cfg; done")
print("  2. Analyze: python3 analyze_fticr_frequencies.py")
print("  3. Validate f_c vs m/z scaling (exact 1/m relationship)")
print("="*80 + "\n")

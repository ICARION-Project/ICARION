#!/usr/bin/env python3
"""
Generate Orbitrap validation suite

Validates axial oscillation frequency vs m/z relationship:
  f_z = (1/2π) * sqrt(q*k/m)

Tests multiple species to verify frequency scaling with mass.
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
OUT_DIR = VALIDATION_DIR / "configs" / "instruments" / "orbitrap"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Orbitrap geometry (realistic Thermo parameters)
# Based on orbitrap_basic.json example
# Note: r_char > r_out is required for positive k parameter!
GEOMETRY = {
    "radius_in_m": 0.006,      # 6 mm inner electrode
    "radius_out_m": 0.015,     # 15 mm outer electrode
    "radius_char_m": 0.022,    # 22 mm characteristic radius (> r_out!)
    "length_m": 0.040,         # 40 mm trap length (z-direction)
    "origin_m": [0.0, 0.0, -0.020]  # Center trap at z=0
}

# Environment (UHV)
ENVIRONMENT = {
    "gas_species": "He",
    "temperature_K": 300.0,
    "pressure_Pa": 1e-6,  # Ultra-high vacuum
    "gas_velocity_m_s": [0.0, 0.0, 0.0]
}

# Electric field parameters
# Orbitrap uses hyperlogarithmic field: U(r,z) = k/2 · (z² - r²/2 + r_char²·ln(r/r_char))
# The k parameter is calculated from radial voltage V: k = 2V / (r_char² ln(r_out/r_in) - 0.5(r_out² - r_in²))
# Using radial_V = 3500V from orbitrap_basic.json example
RADIAL_VOLTAGE_V = 3500.0  # DC radial voltage (realistic Thermo Orbitrap value)

# Test species (wide mass range)
TEST_SPECIES = [
    {"name": "H3O+", "mass_amu": 19.0, "count": 50},
    {"name": "PentanalH+", "mass_amu": 87.0, "count": 50},
    {"name": "CaffeineH+", "mass_amu": 195.0, "count": 50},
    {"name": "ReserpineH+", "mass_amu": 609.0, "count": 50},
]

# Simulation parameters
SIM_TIME_S = 100e-6  # 100 µs (enough for multiple oscillations)
DT_S = 1e-9  # 1 ns timestep
WRITE_INTERVAL = 100  # Write every 100 ns

print("\n" + "="*80)
print("Orbitrap Validation Suite Generator")
print("="*80)
# Calculate actual k parameter from voltage and geometry
r_in = GEOMETRY['radius_in_m']
r_out = GEOMETRY['radius_out_m']
r_char = GEOMETRY['radius_char_m']
r_char_sq = r_char**2
k_actual = 2.0 * RADIAL_VOLTAGE_V / (r_char_sq * np.log(r_out/r_in) - 0.5 * (r_out**2 - r_in**2))

print(f"Radial voltage: {RADIAL_VOLTAGE_V:.2f} V → k ≈ {k_actual:.0f} V/m²")
print(f"Geometry: r_in={r_in*1000:.1f} mm, r_out={r_out*1000:.1f} mm, "
      f"r_char={r_char*1000:.1f} mm, L={GEOMETRY['length_m']*1000:.1f} mm")
print(f"Test species: {len(TEST_SPECIES)}")
print(f"Simulation: {SIM_TIME_S*1e6:.0f} µs")
print("="*80 + "\n")

configs_generated = 0

# Calculate expected frequencies for each species
print("Expected axial frequencies:")
print("-" * 60)
q = 1.602176634e-19  # Elementary charge
for species in TEST_SPECIES:
    m_kg = species["mass_amu"] * 1.66054e-27
    f_z_theory = (1/(2*np.pi)) * np.sqrt(q * k_actual / m_kg)
    print(f"{species['name']:15s} (m={species['mass_amu']:6.1f} u): "
          f"f_z = {f_z_theory/1000:.2f} kHz")
print("-" * 60 + "\n")

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
        "folder": "validation/results/v1.0.0_test/instruments/orbitrap",
        "trajectory_file": f"orbitrap_multi_species_V{RADIAL_VOLTAGE_V:.2f}.h5",
        "print_progress": True
    },
    "ions": {
        "species": []
    },
    "domains": [{
        "domain_index": 0,
        "name": f"Orbitrap V={RADIAL_VOLTAGE_V:.2f}V (k≈{k_actual:.0f} V/m²)",
        "instrument": "Orbitrap",
        "geometry": GEOMETRY,
        "env": ENVIRONMENT,
        "fields": {
            "DC": {
                "radial_V": RADIAL_VOLTAGE_V
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
            "std": [0.0005, 0.0005, 0.001]  # Small spread
        },
        "velocity": {
            "type": "thermal",
            "temperature_K": 300.0
        }
    })

# Write multi-species config
out_file = OUT_DIR / f"orbitrap_multi_species_V{RADIAL_VOLTAGE_V:.2f}.json"
with open(out_file, 'w') as f:
    json.dump(config, f, indent=2)

configs_generated += 1
print(f"✅ Generated multi-species configuration:")
print(f"   {out_file}")
print(f"   Total ions: {sum(s['count'] for s in TEST_SPECIES)}")
print(f"   Species: {', '.join(s['name'] for s in TEST_SPECIES)}")

# Also generate individual single-species configs for detailed analysis
for species in TEST_SPECIES:
    config_single = config.copy()
    config_single["output"]["trajectory_file"] = f"orbitrap_{species['name']}_V{RADIAL_VOLTAGE_V:.2f}.h5"
    config_single["ions"] = {
        "species": [{
            "species_id": species["name"],
            "count": 100,  # More ions for single-species
            "position": {
                "type": "gaussian",
                "center": [0.0, 0.0, 0.0],
                "std": [0.0005, 0.0005, 0.001]
            },
            "velocity": {
                "type": "thermal",
                "temperature_K": 300.0
            }
        }]
    }
    
    out_file = OUT_DIR / f"orbitrap_{species['name']}_V{RADIAL_VOLTAGE_V:.2f}.json"
    with open(out_file, 'w') as f:
        json.dump(config_single, f, indent=2)
    
    configs_generated += 1

print(f"\n✅ Generated {len(TEST_SPECIES)} single-species configurations")
print(f"\n{'='*80}")
print(f"Total configurations: {configs_generated}")
print(f"Output: {OUT_DIR}")
print("\nNext steps:")
print("  1. Run: for cfg in configs/instruments/orbitrap/*.json; do")
print("             icarion_main $cfg; done")
print("  2. Analyze: python3 analyze_orbitrap_frequencies.py")
print("  3. Validate f_z vs m/z scaling")
print("="*80 + "\n")

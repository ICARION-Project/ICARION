#!/usr/bin/env python3
"""
Generate LQIT mass scan configuration

Simulates a mass-selective scan by linearly sweeping AC excitation frequency.
Tests 3 species with different m/z ratios to validate resonant ejection.

Theory:
- Secular frequency: f_sec ≈ q·Ω / (2√2·2π) where Ω = 2π·f_RF
- For H3O+ (m/z=19), f_sec ≈ 142 kHz
- For larger m/z, f_sec decreases (inversely with √m)
"""

import json
import numpy as np
from pathlib import Path

# Output directory
OUT_DIR = Path("../configs/instruments/lqit")
OUT_DIR.mkdir(parents=True, exist_ok=True)

# LQIT geometry
GEOMETRY = {
    "radius_m": 0.005,  # 5 mm
    "length_m": 0.05,   # 50 mm
    "origin_m": [0.0, 0.0, -0.025]
}

# Environment (UHV)
ENVIRONMENT = {
    "gas_species": "He",
    "temperature_K": 300.0,
    "pressure_Pa": 0.1,
    "gas_velocity_m_s": [0.0, 0.0, 0.0]
}

# RF parameters (FIXED voltage for all species → different q values)
RF_FREQUENCY_HZ = 1.0e6  # 1 MHz
RF_VOLTAGE = 20.0  # Fixed RF voltage (lower q for all species)

# Multi-species test (separate domains, same RF, different masses → different q)
SPECIES = [
    {"id": "H3O+", "m": 19.02, "count": 300},
    {"id": "CaffeineH+", "m": 195.0, "count": 300},
    {"id": "ReserpineH+", "m": 609.0, "count": 300}
]

# AC sweep parameters
AC_VOLTAGE = 0.3  # V (gentle excitation)
F_START_HZ = 10e3   # 10 kHz (start lower)
F_END_HZ = 300e3    # 300 kHz
T_SCAN_S = 500e-6   # 500 µs scan time (slower sweep)

def calc_rf_voltage(q_target, m_amu, radius_m, f_rf):
    """Calculate RF voltage for target q parameter"""
    m_kg = m_amu * 1.66054e-27
    q_C = 1.602176634e-19
    omega = 2 * np.pi * f_rf
    V_rf = (q_target * m_kg * omega**2 * radius_m**2) / (4 * q_C)
    return V_rf

def calc_secular_freq(q_param, f_rf):
    """Calculate theoretical secular frequency"""
    omega_rf = 2 * np.pi * f_rf
    omega_sec = (q_param * omega_rf) / (2 * np.sqrt(2))
    return omega_sec / (2 * np.pi)

print("\n" + "="*70)
print("LQIT Mass Scan Configuration (Fixed RF Voltage)")
print("="*70)
print(f"RF: {RF_FREQUENCY_HZ/1e6:.2f} MHz, {RF_VOLTAGE:.1f} V (fixed)")
print(f"AC Sweep: {F_START_HZ/1e3:.0f} - {F_END_HZ/1e3:.0f} kHz over {T_SCAN_S*1e6:.0f} µs")
print(f"\nExpected parameters (fixed V → q ∝ 1/m → f_sec ∝ 1/√m):")

for species in SPECIES:
    # Calculate q for this mass at fixed RF voltage
    m_kg = species["m"] * 1.66054e-27
    q_C = 1.602176634e-19
    omega = 2 * np.pi * RF_FREQUENCY_HZ
    r0 = GEOMETRY["radius_m"]
    
    q_param = (4 * q_C * RF_VOLTAGE) / (m_kg * omega**2 * r0**2)
    f_sec = calc_secular_freq(q_param, RF_FREQUENCY_HZ)
    
    species["q_param"] = q_param
    species["f_sec_theory"] = f_sec
    
    status = "STABLE" if q_param < 0.908 else "UNSTABLE"
    print(f"  {species['id']:25s} m/z={species['m']:6.1f}, q={q_param:.3f} {status:8s}, f_sec={f_sec/1e3:6.1f} kHz")

print("="*70 + "\n")

# Build configuration
config = {
    "simulation": {
        "total_time_s": T_SCAN_S,
        "dt_s": 1e-9,  # 1 ns timestep (balance speed/accuracy)
        "write_interval": 100,  # Every 100 ns for tracking
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
        "folder": "results/v1.0_test/instruments/lqit",
        "trajectory_file": "lqit_mass_scan.h5",
        "print_progress": True
    },
    "ions": {
        "species": []
    },
    "domains": []
}

# Create separate domain for each species (same RF voltage → mass-dependent q)
for i, species in enumerate(SPECIES):
    # All domains use SAME RF voltage (mass scan condition)
    domain = {
        "domain_index": i,
        "name": f"LQIT {species['id']}",
        "instrument": "LQIT",
        "geometry": GEOMETRY,
        "env": ENVIRONMENT,
        "fields": {
            "RF": {
                "voltage_V": RF_VOLTAGE,
                "frequency_Hz": RF_FREQUENCY_HZ,
                "phase_rad": 0.0
            },
            "DC": {
                "axial_V": 10.0,
                "quad_V": 0.0
            },
            "AC": {
                "voltage_V": AC_VOLTAGE,
                "frequency_Hz": {
                    "type": "linear",
                    "start": F_START_HZ,
                    "end": F_END_HZ,
                    "start_time_s": 0.0,
                    "end_time_s": T_SCAN_S,
                    "clamp": True
                }
            }
        }
    }
    config["domains"].append(domain)
    
    # Add ions for this domain
    config["ions"]["species"].append({
        "species_id": species["id"],
        "count": species["count"],
        "domain_index": i,
        "position": {
            "type": "gaussian",
            "center": [0.0, 0.0, 0.0],  # All at trap center
            "std": [0.0003, 0.0003, 0.001]
        },
        "velocity": {
            "type": "thermal",
            "temperature_K": 300.0
        }
    })

# Write config
out_file = OUT_DIR / "lqit_mass_scan.json"
with open(out_file, 'w') as f:
    json.dump(config, f, indent=2)

print(f"✅ Generated: {out_file.name}")
print(f"   Species: {', '.join(s['id'] for s in SPECIES)}")
print(f"   Total ions: {sum(s['count'] for s in SPECIES)}")
print(f"   Scan rate: {(F_END_HZ - F_START_HZ) / T_SCAN_S / 1e9:.2f} GHz/s\n")

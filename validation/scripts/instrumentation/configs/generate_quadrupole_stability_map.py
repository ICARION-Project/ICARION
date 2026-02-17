#!/usr/bin/env python3
"""
Generate Quadrupole Stability Map validation suite

Creates a systematic sweep over Mathieu (a,q) parameter space to validate
stability boundaries and compare with theory.

Theory (first stability region):
  - 0 < q < 0.908
  - |a| < 0.2 (for q near 0.908)
  - Stability boundary: βx,βy = 0 and βx,βy = 1
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
OUT_DIR = VALIDATION_DIR / "configs" / "instruments" / "quadrupole_stability"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# Quadrupole geometry (LQIT parameters)
GEOMETRY = {
    "radius_m": 0.005,
    "length_m": 0.05,
    "origin_m": [0.0, 0.0, -0.025]
}

# Environment (UHV for clean dynamics)
ENVIRONMENT = {
    "gas_species": "He",
    "temperature_K": 300.0,
    "pressure_Pa": 0.01,  # Very low pressure (collisionless)
    "gas_velocity_m_s": [0.0, 0.0, 0.0]
}

# Ion parameters
ION_SPECIES = "CaffeineH+"  # m/z = 195 (good mid-range mass)
ION_MASS_AMU = 195.0
ION_COUNT = 100  # Smaller for faster sims
RF_FREQUENCY_HZ = 1.0e6  # 1 MHz

# Simulation time (need enough cycles to see instability)
SIM_TIME_S = 50e-6  # 50 µs = 50 RF cycles

# (a,q) parameter grid
# Focus on first stability region boundaries
Q_VALUES = [
    0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.88, 0.90, 0.91, 0.92, 0.95, 1.0
]

A_VALUES = [
    -0.20, -0.15, -0.10, -0.05, 0.0, 0.05, 0.10, 0.15, 0.20
]

def calc_voltages(a, q, m_amu, r0, omega_rf):
    """Calculate RF and DC voltages for target (a,q) parameters
    
    q = 4QV_rf / (m·Ω²·r₀²)
    a = -8QV_dc / (m·Ω²·r₀²)
    """
    m_kg = m_amu * 1.66054e-27
    Q = 1.602176634e-19
    
    V_rf = (q * m_kg * omega_rf**2 * r0**2) / (4 * Q)
    V_dc = -(a * m_kg * omega_rf**2 * r0**2) / (8 * Q)
    
    return V_rf, V_dc

def is_stable_theory(a, q):
    """Theoretical stability prediction (first region)"""
    # First stability region approximation
    if q <= 0:
        return False
    if q > 0.908:
        return False
    
    # Approximate boundaries (simplified)
    # More accurate would require solving Mathieu equations
    if abs(a) > 0.2:
        return False
    
    # Region near q=0.908 requires a≈0
    if q > 0.85 and abs(a) > 0.1:
        return False
    
    return True

print("\n" + "="*80)
print("Quadrupole Stability Map Validation Suite")
print("="*80)
print(f"Ion: {ION_SPECIES} (m/z = {ION_MASS_AMU})")
print(f"RF frequency: {RF_FREQUENCY_HZ/1e6:.2f} MHz")
print(f"Grid: {len(Q_VALUES)} q × {len(A_VALUES)} a = {len(Q_VALUES)*len(A_VALUES)} points")
print(f"Simulation: {SIM_TIME_S*1e6:.0f} µs ({int(SIM_TIME_S*RF_FREQUENCY_HZ)} RF cycles)")
print("="*80 + "\n")

omega_rf = 2 * np.pi * RF_FREQUENCY_HZ
r0 = GEOMETRY["radius_m"]

configs_generated = 0
stable_theory = 0
unstable_theory = 0

for q in Q_VALUES:
    for a in A_VALUES:
        # Calculate required voltages
        V_rf, V_dc = calc_voltages(a, q, ION_MASS_AMU, r0, omega_rf)
        
        # Skip if voltages are unreasonable
        if abs(V_rf) > 2000 or abs(V_dc) > 500:
            continue
        
        # Theoretical prediction
        stable = is_stable_theory(a, q)
        if stable:
            stable_theory += 1
        else:
            unstable_theory += 1
        
        # Build config
        config = {
            "simulation": {
                "total_time_s": SIM_TIME_S,
                "dt_s": 1e-9,  # 1 ns
                "write_interval": 50,  # Every 50 ns
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
                "folder": "validation/results/v1.0.0_test/instruments/quadrupole_stability",
                "trajectory_file": f"quad_q{q:.3f}_a{a:+.3f}.h5",
                "print_progress": False
            },
            "ions": {
                "species": [{
                    "species_id": ION_SPECIES,
                    "count": ION_COUNT,
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
            },
            "domains": [{
                "domain_index": 0,
                "name": f"Quadrupole q={q:.3f} a={a:+.3f}",
                "instrument": "LQIT",
                "geometry": GEOMETRY,
                "env": ENVIRONMENT,
                "fields": {
                    "RF": {
                        "voltage_V": V_rf,
                        "frequency_Hz": RF_FREQUENCY_HZ,
                        "phase_rad": 0.0
                    },
                    "DC": {
                        "quad_V": V_dc,
                        "axial_V": 10.0  # Axial confinement
                    },
                    "AC": {
                        "voltage_V": 0.0,  # No AC excitation
                        "frequency_Hz": 0.0
                    }
                }
            }]
        }
        
        # Write config
        status = "stable" if stable else "unstable"
        out_file = OUT_DIR / f"quad_q{q:.3f}_a{a:+.3f}_{status}.json"
        with open(out_file, 'w') as f:
            json.dump(config, f, indent=2)
        
        configs_generated += 1

print(f"✅ Generated {configs_generated} configurations")
print(f"   Theoretically stable:   {stable_theory}")
print(f"   Theoretically unstable: {unstable_theory}")
print(f"\nOutput: {OUT_DIR}")
print("\nNext steps:")
print("  1. Run batch: for cfg in configs/instruments/quadrupole_stability/*.json; do")
print("                  icarion_main $cfg > /dev/null 2>&1 & done; wait")
print("  2. Analyze: python3 analyze_quadrupole_stability_map.py")
print("  3. Compare with theoretical (a,q) stability diagram")
print()

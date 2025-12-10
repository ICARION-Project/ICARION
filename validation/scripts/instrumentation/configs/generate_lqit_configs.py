#!/usr/bin/env python3
"""
Generate LQIT validation test configurations

Tests Mathieu stability, axial frequency, and ion confinement in a Linear Quadrupole Ion Trap.
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
OUT_DIR = VALIDATION_DIR / "configs" / "instruments" / "lqit"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# LQIT geometry (typical 2D quadrupole trap)
GEOMETRY = {
    "radius_m": 0.005,  # 5 mm inscribed radius (r0)
    "length_m": 0.05,  # 50 mm length
    "origin_m": [0.0, 0.0, -0.025]  # z origin at -L/2 (trap center at z=0)
}

# Environment (UHV conditions)
ENVIRONMENT = {
    "gas_species": "He",
    "temperature_K": 300.0,
    "pressure_Pa": 0.1,  # 10^-3 mbar UHV
    "gas_velocity_m_s": [0.0, 0.0, 0.0]
}

# Ion parameters
ION_SPECIES = "H3O+"
ION_COUNT = 1000
ION_MASS_AMU = 19.02

# RF drive frequency (typical for LQIT)
RF_FREQUENCY_HZ = 1.0e6  # 1 MHz

def calc_rf_voltage_for_q(q_target, m_amu, q_charge, radius_m, omega_rf):
    """Calculate RF voltage needed to achieve target q parameter"""
    m_kg = m_amu * 1.66054e-27
    q_C = q_charge * 1.602176634e-19
    omega = 2 * np.pi * omega_rf
    
    # q = 4·Q·V_rf / (m·Ω²·r₀²)
    # V_rf = q·m·Ω²·r₀² / (4·Q)
    V_rf = (q_target * m_kg * omega**2 * radius_m**2) / (4 * q_C)
    
    return V_rf

def calc_mathieu_params(m_amu, q_charge, radius_m, V_rf, V_dc, omega_rf):
    """
    Calculate Mathieu stability parameters a and q
    
    For LQIT (2D quadrupole):
    q = 4·Q·V_rf / (m·Ω²·r₀²)
    a = -8·Q·V_dc / (m·Ω²·r₀²)
    """
    m_kg = m_amu * 1.66054e-27
    q_C = q_charge * 1.602176634e-19
    omega = 2 * np.pi * omega_rf
    
    q_param = (4 * q_C * V_rf) / (m_kg * omega**2 * radius_m**2)
    a_param = (-8 * q_C * V_dc) / (m_kg * omega**2 * radius_m**2)
    
    return a_param, q_param

def calc_secular_frequency(q_param, omega_rf):
    """
    Calculate secular (pseudopotential) frequency
    
    For small q: ω_sec ≈ q·ω_rf / (2√2)
    """
    return (q_param * omega_rf) / (2 * np.sqrt(2))

def generate_config(rf_voltage_V, dc_voltage_V, collision_model, description):
    """Generate a single LQIT test configuration"""
    
    # Calculate Mathieu parameters
    a, q = calc_mathieu_params(
        m_amu=ION_MASS_AMU,
        q_charge=1,
        radius_m=GEOMETRY["radius_m"],
        V_rf=rf_voltage_V,
        V_dc=dc_voltage_V,
        omega_rf=RF_FREQUENCY_HZ
    )
    
    # Calculate secular frequency
    omega_sec = calc_secular_frequency(q, 2*np.pi*RF_FREQUENCY_HZ)
    f_sec_Hz = omega_sec / (2*np.pi)
    
    # Simulation time: 20 RF periods for stability check
    t_rf_period = 1.0 / RF_FREQUENCY_HZ
    t_total = 20 * t_rf_period
    
    # Timestep: 100 points per RF period
    dt = t_rf_period / 100
    
    # Write interval: every RF period
    write_interval = 100
    
    config = {
        "simulation": {
            "total_time_s": t_total,
            "dt_s": dt,
            "write_interval": write_interval,
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
            "folder": "results/v1.0_test/instruments/lqit",
            "trajectory_file": f"lqit_{collision_model.lower()}_q{q:.3f}_a{a:.3f}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [{
                "species_id": ION_SPECIES,
                "count": ION_COUNT,
                "position": {
                    "type": "gaussian",
                    "center": [0.0, 0.0, 0.0],  # Center at trap origin
                    "std": [0.0005, 0.0005, 0.001]  # 0.5mm radial, 1mm axial
                },
                "velocity": {
                    "type": "thermal",
                    "temperature_K": 300.0
                }
            }]
        },
        "domains": [{
            "domain_index": 0,
            "name": f"LQIT {description}",
            "instrument": "LQIT",
            "geometry": GEOMETRY,
            "env": ENVIRONMENT,
            "fields": {
                "RF": {
                    "frequency_Hz": RF_FREQUENCY_HZ,
                    "voltage_V": rf_voltage_V,
                    "phase_rad": 0.0
                },
                "DC": {
                    "quad_V": dc_voltage_V,
                    "axial_V": 10.0  # Endcap voltage for axial confinement
                }
            }
        }]
    }
    
    return config, a, q, f_sec_Hz

def generate_config_with_ac(rf_voltage_V, dc_voltage_V, ac_voltage_V, ac_freq_Hz, collision_model, description):
    """Generate LQIT config with AC excitation waveform"""
    
    # Get base config
    config, a, q, f_sec_Hz = generate_config(rf_voltage_V, dc_voltage_V, collision_model, description)
    
    # Add AC waveform for resonant excitation
    config["domains"][0]["fields"]["AC"] = {
        "frequency_Hz": ac_freq_Hz,
        "voltage_V": ac_voltage_V,
        "phase_rad": 0.0
    }
    
    # Longer simulation for AC excitation (100 secular periods)
    t_sec_period = 1.0 / f_sec_Hz
    config["simulation"]["total_time_s"] = 100 * t_sec_period
    
    # Update output filename
    config["output"]["trajectory_file"] = f"lqit_{collision_model.lower()}_q{q:.3f}_ac{ac_freq_Hz/1e3:.0f}kHz.h5"
    
    return config, a, q, f_sec_Hz

def main():
    """Generate LQIT validation test matrix"""
    
    print(f"Generating LQIT validation configs in {OUT_DIR}")
    
    # Calculate RF voltages for target q values
    V_rf_q04 = calc_rf_voltage_for_q(0.4, ION_MASS_AMU, 1, GEOMETRY["radius_m"], RF_FREQUENCY_HZ)
    V_rf_q07 = calc_rf_voltage_for_q(0.7, ION_MASS_AMU, 1, GEOMETRY["radius_m"], RF_FREQUENCY_HZ)
    V_rf_q095 = calc_rf_voltage_for_q(0.95, ION_MASS_AMU, 1, GEOMETRY["radius_m"], RF_FREQUENCY_HZ)
    
    # Test matrix: Stability + AC Excitation
    tests = [
        # (RF_V, DC_V, collision_model, description)
        # Stable region (q=0.4, a=0.0)
        (V_rf_q04, 0.0, "HSS", "Stable q=0.4"),
        (V_rf_q04, 0.0, "EHSS", "Stable q=0.4"),
        (V_rf_q04, 0.0, "Friction", "Stable q=0.4"),
        
        # Stable region (q=0.7, a=0.0) - close to boundary
        (V_rf_q07, 0.0, "HSS", "Stable q=0.7"),
        
        # Unstable region (q=0.95, a=0.0) - beyond stability
        (V_rf_q095, 0.0, "HSS", "Unstable q=0.95"),
        
        # Stable with DC (q=0.4, a=0.05) - calculate DC voltage for a=0.05
        (V_rf_q04, 0.24, "HSS", "Stable with DC"),  # V_dc chosen to give a≈0.05
    ]
    
    # AC Excitation tests (resonant excitation at secular frequency)
    ac_tests = [
        # (RF_V, DC_V, AC_V, AC_freq_factor, collision_model, description)
        # On-resonance: AC freq = secular freq
        (V_rf_q04, 0.0, 1.0, 1.0, "HSS", "On-resonance AC"),
        (V_rf_q04, 0.0, 1.0, 1.0, "Friction", "On-resonance AC"),
        
        # Off-resonance: AC freq = 0.5 × secular freq
        (V_rf_q04, 0.0, 1.0, 0.5, "HSS", "Off-resonance AC (0.5x)"),
        
        # Off-resonance: AC freq = 2.0 × secular freq
        (V_rf_q04, 0.0, 1.0, 2.0, "HSS", "Off-resonance AC (2.0x)"),
    ]
    
    configs_generated = []
    
    print(f"\nTest matrix: {len(tests)} stability tests + {len(ac_tests)} AC excitation tests")
    print(f"Ion: {ION_SPECIES} ({ION_MASS_AMU} u)")
    print(f"RF frequency: {RF_FREQUENCY_HZ/1e6:.1f} MHz")
    print(f"Trap radius: {GEOMETRY['radius_m']*1000:.1f} mm")
    print()
    
    print("=== Stability Tests ===")
    for i, (rf_V, dc_V, model, desc) in enumerate(tests, 1):
        config, a, q, f_sec = generate_config(rf_V, dc_V, model, desc)
        
        # Determine stability
        stable = (q < 0.908) and (abs(a) < 0.2)  # Simplified stability boundary
        status = "STABLE" if stable else "UNSTABLE"
        
        filename = OUT_DIR / f"lqit_{model.lower()}_q{q:.3f}_a{abs(a):.3f}.json"
        with open(filename, 'w') as f:
            json.dump(config, f, indent=2)
        
        configs_generated.append(filename)
        
        print(f"{i:2d}. {model:8s} | q={q:.3f}, a={a:+.3f} | f_sec={f_sec/1e3:.1f} kHz | {status:8s} | {desc}")
    
    print(f"\n=== AC Excitation Tests ===")
    for i, (rf_V, dc_V, ac_V, ac_factor, model, desc) in enumerate(ac_tests, 1):
        # First calculate secular frequency
        _, a, q, f_sec = generate_config(rf_V, dc_V, model, "temp")
        ac_freq = ac_factor * f_sec
        
        # Generate config with AC
        config, a, q, f_sec = generate_config_with_ac(rf_V, dc_V, ac_V, ac_freq, model, desc)
        
        filename = OUT_DIR / f"lqit_{model.lower()}_q{q:.3f}_ac{ac_freq/1e3:.0f}kHz.json"
        with open(filename, 'w') as f:
            json.dump(config, f, indent=2)
        
        configs_generated.append(filename)
        
        resonance = "ON-RES" if abs(ac_factor - 1.0) < 0.01 else "OFF-RES"
        print(f"{i:2d}. {model:8s} | AC={ac_freq/1e3:.1f} kHz ({ac_factor:.1f}×f_sec) | V_AC={ac_V:.1f} V | {resonance:7s} | {desc}")
    
    # Generate README
    readme_path = OUT_DIR / "README.md"
    with open(readme_path, 'w') as f:
        f.write(f"""# LQIT Validation Test Suite

## Overview
Tests Linear Quadrupole Ion Trap (LQIT) stability, secular motion, and resonant excitation.

## Test Parameters
- **Ion Species:** {ION_SPECIES} ({ION_MASS_AMU} u, +1 e)
- **Ion Count:** {ION_COUNT}
- **Trap Geometry:** r₀ = {GEOMETRY['radius_m']*1000:.1f} mm, L = {GEOMETRY['length_m']*1000:.1f} mm
- **RF Frequency:** {RF_FREQUENCY_HZ/1e6:.1f} MHz
- **Environment:** {ENVIRONMENT['gas_species']} @ {ENVIRONMENT['pressure_Pa']:.1f} Pa (UHV)

## Test Suite

### Part 1: Mathieu Stability ({len(tests)} tests)
1. **Stable Region (q=0.4):** All collision models
2. **Stability Boundary (q=0.7):** Near-critical stability
3. **Unstable Region (q=0.95):** Expected ion ejection
4. **DC Offset:** Stability with axial confinement

### Part 2: AC Resonant Excitation ({len(ac_tests)} tests)
1. **On-Resonance (f_AC = f_sec):** Strong heating/ejection
2. **Off-Resonance (f_AC = 0.5×f_sec):** Minimal effect
3. **Off-Resonance (f_AC = 2.0×f_sec):** Minimal effect
4. **Collision Model Comparison:** HSS vs Friction damping

## Expected Results

### Stability Tests
- **Stable (q=0.4):** All ions confined, oscillate at f_sec ≈ 141 kHz
- **Stable (q=0.7):** Ions confined, f_sec ≈ 247 kHz  
- **Unstable (q=0.95):** Ions ejected within 5-10 RF periods

### AC Excitation Tests
- **On-Resonance:** Ion amplitude grows exponentially → ejection
- **Off-Resonance:** Ion amplitude stable (no parametric resonance)
- **Collision Damping:** HSS shows stronger damping than Friction in UHV

## Analysis
```bash
# Stability analysis
python3 ../scripts/analyze_lqit_stability.py results/v1.0_test/instruments/lqit/lqit_*_q*.h5

# AC excitation analysis
python3 ../scripts/analyze_lqit_excitation.py results/v1.0_test/instruments/lqit/lqit_*_ac*.h5
```

## Physics Background

### Mathieu Stability
For 2D quadrupole (LQIT):
```
q = 4·Q·V_rf / (m·Ω²·r₀²)
a = -8·Q·V_dc / (m·Ω²·r₀²)
```
Stability: q < 0.908, |a| < 0.2

### Secular Frequency
Pseudopotential approximation (q << 1):
```
ω_sec ≈ q·Ω / (2√2)
```

### Resonant Excitation
AC field at f_sec causes parametric resonance → amplitude growth → ejection

## References
- March & Todd, "Quadrupole Ion Trap Mass Spectrometry" (2005)
- Dehmelt pseudopotential approximation
- Mathieu equation solutions
""")
    
    print(f"\n✅ Generated {len(configs_generated)} LQIT validation configurations")
    print(f"✅ Created README: {readme_path}")

if __name__ == "__main__":
    main()

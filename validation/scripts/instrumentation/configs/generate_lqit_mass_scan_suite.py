#!/usr/bin/env python3
"""
Generate realistic LQIT Mass Scan validation suite

Tests mass-selective ejection at different scan rates matching real QIT instruments.
Validates:
  - Mass resolution vs scan rate
  - Ejection frequency accuracy
  - Peak shapes and widths
  - Scan rate dependence (1-100 kDa/s range)
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
OUT_DIR = VALIDATION_DIR / "configs" / "instruments" / "lqit_scans"
OUT_DIR.mkdir(parents=True, exist_ok=True)

# LQIT geometry
GEOMETRY = {
    "radius_m": 0.005,
    "length_m": 0.05,
    "origin_m": [0.0, 0.0, -0.025]
}

# Environment (UHV)
ENVIRONMENT = {
    "gas_species": "He",
    "temperature_K": 300.0,
    "pressure_Pa": 0.1,
    "gas_velocity_m_s": [0.0, 0.0, 0.0]
}

# RF parameters (optimized for q ~ 0.2-0.5 range for all species)
RF_FREQUENCY_HZ = 1.0e6  # 1 MHz
RF_VOLTAGE = 50.0  # V (keeps H3O+ stable with q~1.0)

# Species for mass range validation (m/z 87-609, factor 7× range)
SPECIES = [
    {"id": "PentanalH+", "m": 87.0},
    {"id": "2,6-DTBPH+", "m": 192.0},  # Good intermediate mass
    {"id": "CaffeineH+", "m": 195.0},
    {"id": "ReserpineH+", "m": 609.0}
]

# Scan rate configurations (matching real QIT performance)
SCAN_CONFIGS = [
    {
        "name": "slow_high_res",
        "description": "Slow scan (high resolution)",
        "scan_time_ms": 500,  # 500 ms
        "scan_rate_kDa_s": 1.2,  # (609-19)/500ms ≈ 1.2 kDa/s
        "ac_voltage": 0.3,  # V (gentle for clean peaks)
        "ion_count": 200
    },
    {
        "name": "normal_balanced",
        "description": "Normal scan (balanced speed/resolution)",
        "scan_time_ms": 100,  # 100 ms
        "scan_rate_kDa_s": 5.9,
        "ac_voltage": 0.5,  # V
        "ion_count": 200
    },
    {
        "name": "fast_low_res",
        "description": "Fast scan (survey mode)",
        "scan_time_ms": 20,  # 20 ms
        "scan_rate_kDa_s": 29.5,
        "ac_voltage": 0.8,  # V (stronger for fast ejection)
        "ion_count": 200
    },
    {
        "name": "ultra_fast",
        "description": "Ultra-fast scan (MS/MS precursor)",
        "scan_time_ms": 5,  # 5 ms
        "scan_rate_kDa_s": 118,
        "ac_voltage": 1.0,  # V (maximum)
        "ion_count": 150
    }
]

def calc_secular_freq(m_amu, rf_voltage, radius_m, f_rf):
    """Calculate theoretical secular frequency"""
    m_kg = m_amu * 1.66054e-27
    q_C = 1.602176634e-19
    omega = 2 * np.pi * f_rf
    q_param = (4 * q_C * rf_voltage) / (m_kg * omega**2 * radius_m**2)
    omega_sec = (q_param * omega) / (2 * np.sqrt(2))
    return omega_sec / (2 * np.pi), q_param

def generate_config(scan_cfg):
    """Generate configuration for one scan rate"""
    
    scan_name = scan_cfg["name"]
    t_scan = scan_cfg["scan_time_ms"] * 1e-3  # Convert to seconds
    ac_voltage = scan_cfg["ac_voltage"]
    ion_count = scan_cfg["ion_count"]
    
    # Calculate frequency range based on species mass range
    # Use safety margin to ensure full coverage
    f_min_theory = calc_secular_freq(SPECIES[-1]["m"], RF_VOLTAGE, GEOMETRY["radius_m"], RF_FREQUENCY_HZ)[0]
    f_max_theory = calc_secular_freq(SPECIES[0]["m"], RF_VOLTAGE, GEOMETRY["radius_m"], RF_FREQUENCY_HZ)[0]
    
    # Add 50% margin
    f_start = max(5e3, f_min_theory * 0.5)  # At least 5 kHz
    f_end = f_max_theory * 1.5
    
    # Calculate timestep (need ~10 points per RF cycle minimum)
    dt = 1.0 / (RF_FREQUENCY_HZ * 10)  # 0.1 µs
    write_interval = max(10, int(t_scan / (dt * 2000)))  # ~2000 output points
    
    print(f"\n{'='*70}")
    print(f"Scan: {scan_cfg['description']}")
    print(f"{'='*70}")
    print(f"Scan time:  {t_scan*1e3:.1f} ms")
    print(f"Scan rate:  {scan_cfg['scan_rate_kDa_s']:.1f} kDa/s")
    print(f"AC voltage: {ac_voltage:.1f} V")
    print(f"Frequency:  {f_start/1e3:.0f} - {f_end/1e3:.0f} kHz")
    print(f"Sweep rate: {(f_end-f_start)/t_scan/1e6:.2f} MHz/s")
    print(f"\nSpecies parameters:")
    
    # Build config
    config = {
        "simulation": {
            "total_time_s": t_scan,
            "dt_s": dt,
            "write_interval": write_interval,
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
            "folder": "results/v1.0_test/instruments/lqit_scans",
            "trajectory_file": f"lqit_scan_{scan_name}.h5",
            "print_progress": True
        },
        "ions": {"species": []},
        "domains": []
    }
    
    # Create domain for each species
    for i, species in enumerate(SPECIES):
        f_sec, q = calc_secular_freq(species["m"], RF_VOLTAGE, GEOMETRY["radius_m"], RF_FREQUENCY_HZ)
        status = "STABLE" if q < 0.908 else "UNSTABLE"
        print(f"  {species['id']:15s} m/z={species['m']:6.1f}, q={q:.3f} {status:8s}, f_sec={f_sec/1e3:6.1f} kHz")
        
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
                    "voltage_V": ac_voltage,
                    "frequency_Hz": {
                        "type": "linear",
                        "start": f_start,
                        "end": f_end,
                        "start_time_s": 0.0,
                        "end_time_s": t_scan,
                        "clamp": True
                    }
                }
            }
        }
        config["domains"].append(domain)
        
        # Add ions for this species
        config["ions"]["species"].append({
            "species_id": species["id"],
            "count": ion_count,
            "domain_index": i,
            "position": {
                "type": "gaussian",
                "center": [0.0, 0.0, 0.0],
                "std": [0.0003, 0.0003, 0.001]
            },
            "velocity": {
                "type": "thermal",
                "temperature_K": 300.0
            }
        })
    
    # Write config
    out_file = OUT_DIR / f"lqit_scan_{scan_name}.json"
    with open(out_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    print(f"\n✅ Generated: {out_file.name}")
    print(f"   Total ions: {len(SPECIES) * ion_count}")
    print(f"   Timestep: {dt*1e9:.1f} ns, Output points: ~{int(t_scan/dt/write_interval)}")

# Generate all configurations
print("\n" + "="*70)
print("LQIT Mass Scan Validation Suite")
print("="*70)
print(f"RF: {RF_FREQUENCY_HZ/1e6:.2f} MHz, {RF_VOLTAGE:.0f} V (fixed)")
print(f"Species: {', '.join(s['id'] for s in SPECIES)}")
print(f"Mass range: {SPECIES[0]['m']:.0f} - {SPECIES[-1]['m']:.0f} u")
print(f"Configurations: {len(SCAN_CONFIGS)}")

for scan_cfg in SCAN_CONFIGS:
    generate_config(scan_cfg)

print("\n" + "="*70)
print(f"✅ Generated {len(SCAN_CONFIGS)} mass scan configurations")
print("="*70)
print("\nNext steps:")
print("  1. Run simulations: for cfg in configs/instruments/lqit_scans/*.json; do icarion_main $cfg; done")
print("  2. Analyze results: python3 analyze_lqit_mass_scan_suite.py")
print("  3. Validate:")
print("     - Mass resolution vs scan rate")
print("     - Ejection frequency accuracy")
print("     - Peak width scaling")
print()

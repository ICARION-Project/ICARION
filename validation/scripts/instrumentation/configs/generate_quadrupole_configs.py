#!/usr/bin/env python3
"""
Generate Quadrupole Stability Map validation configurations

Tests ion transmission through quadrupole mass filter across (a,q) parameter space
to validate against Mathieu equation stability diagram.

Mathieu Stability Parameters:
    a = 8eU / (m r₀² Ω²)
    q = 4eV / (m r₀² Ω²)

where:
    U = DC voltage
    V = RF amplitude (zero-to-peak)
    r₀ = field radius
    Ω = 2πf (angular frequency)
    m = ion mass
    e = elementary charge
"""

import json
import numpy as np
from pathlib import Path

# Quadrupole geometry
R0_M = 0.005  # 5 mm field radius (typical for analytical quad)
LENGTH_M = 0.05  # 50 mm rod length
RADIUS_M = 0.010  # 10 mm chamber radius

# RF parameters
RF_FREQ_HZ = 2.0e6  # 2 MHz (typical)
OMEGA = 2 * np.pi * RF_FREQ_HZ

# Ion parameters
ION_SPECIES = "CaffeineH+"  # m = 195.08 amu
MASS_AMU = 195.08
MASS_KG = MASS_AMU * 1.66053906660e-27

# Physical constants
E_CHARGE = 1.602176634e-19  # C

# Simulation parameters
N_IONS = 100  # ions per test
SIMULATION_TIME_S = 50e-6  # 50 µs (25 RF cycles at 2 MHz)
DT_S = 1e-9  # 1 ns timestep
WRITE_INTERVAL = 500  # every 500 steps

# Output directory
CONFIG_DIR = Path(__file__).parent.parent.parent / "configs" / "instruments" / "quadrupole"

def calc_voltages_from_aq(a, q, m_amu, r0_m, omega):
    """Calculate U (DC) and V (RF) from Mathieu parameters a and q"""
    m_kg = m_amu * 1.66053906660e-27
    
    # From q = 4eV / (m r₀² Ω²)
    V = q * m_kg * r0_m**2 * omega**2 / (4 * E_CHARGE)
    
    # From a = 8eU / (m r₀² Ω²)
    U = a * m_kg * r0_m**2 * omega**2 / (8 * E_CHARGE)
    
    return U, V

def is_stable_mathieu(a, q):
    """Check if (a,q) point is in first stability region (approximate)"""
    # First stability region boundaries (simplified)
    # More accurate would use actual Mathieu stability boundaries
    
    if q < 0 or q > 0.908:
        return False
    
    # Upper boundary: a ≈ 0.237 - 0.26*q + ...(higher order terms)
    a_upper = 0.237 - 0.26 * q
    
    # Lower boundary
    a_lower = -0.1
    
    return a_lower < a < a_upper

def generate_config(a, q, config_num, total_configs):
    """Generate a single quadrupole stability test configuration"""
    
    U, V = calc_voltages_from_aq(a, q, MASS_AMU, R0_M, OMEGA)
    
    # Determine if this should be stable
    stable_theory = is_stable_mathieu(a, q)
    
    config = {
        "simulation": {
            "total_time_s": SIMULATION_TIME_S,
            "dt_s": DT_S,
            "write_interval": WRITE_INTERVAL,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "num_threads": 4,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "NoCollisions",
            "enable_space_charge": False,
            "enable_reactions": False
        },
        "output": {
            "folder": "results/v1.0_test/instruments/quadrupole",
            "trajectory_file": f"quad_a{a:.4f}_q{q:.4f}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "species_id": ION_SPECIES,
                    "count": N_IONS,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, -LENGTH_M/2],
                        "std": [0.0002, 0.0002, 0.0005]  # Small initial spread
                    },
                    "velocity": {
                        "type": "kinetic",
                        "energy_eV": 5.0,  # 5 eV kinetic energy
                        "direction": [0.0, 0.0, 1.0],
                        "spread_angle_deg": 2.0  # Small angular spread
                    }
                }
            ]
        },
        "domains": [
            {
                "domain_index": 0,
                "name": "Quadrupole Mass Filter",
                "instrument": "Quadrupole",
                "geometry": {
                    "origin_m": [0.0, 0.0, -LENGTH_M],
                    "length_m": LENGTH_M,
                    "radius_m": R0_M  # USE FIELD RADIUS, NOT CHAMBER RADIUS!
                },
                "env": {
                    "gas_species": "He",
                    "temperature_K": 300.0,
                    "pressure_Pa": 0.001,  # ~1e-5 mbar (high vacuum)
                    "gas_velocity_m_s": [0.0, 0.0, 0.0]
                },
                "fields": {
                    "RF": {
                        "voltage_V": V,
                        "frequency_Hz": RF_FREQ_HZ,
                        "phase_rad": 0.0
                    },
                    "DC": {
                        "quad_V": U,
                        "axial_V": 50.0  # Small DC offset for ion transport
                    }
                }
            }
        ]
    }
    
    return config, stable_theory

def main():
    """Generate quadrupole stability map configurations"""
    
    CONFIG_DIR.mkdir(parents=True, exist_ok=True)
    
    print("Generating Quadrupole Stability Map validation configs")
    print(f"Output: {CONFIG_DIR}")
    print(f"Ion: {ION_SPECIES} (m = {MASS_AMU} amu)")
    print(f"Quadrupole: r₀ = {R0_M*1000} mm, f = {RF_FREQ_HZ/1e6} MHz")
    print()
    
    # Define (a, q) grid
    # First stability region + instability verification: q ∈ [0, 1.0], a ∈ [-0.1, 0.237]
    q_values = np.linspace(0.05, 1.00, 11)  # 11 points from 0.05 to 1.00 (includes unstable region!)
    a_values = np.linspace(-0.05, 0.23, 8)  # 8 points from -0.05 to 0.23
    
    configs = []
    stable_count = 0
    unstable_count = 0
    
    for q in q_values:
        for a in a_values:
            U, V = calc_voltages_from_aq(a, q, MASS_AMU, R0_M, OMEGA)
            stable_theory = is_stable_mathieu(a, q)
            
            if stable_theory:
                stable_count += 1
            else:
                unstable_count += 1
            
            configs.append((a, q, U, V, stable_theory))
    
    total_configs = len(configs)
    
    print("Stability Map Grid:")
    print(f"  q-values: {len(q_values)} points from {q_values[0]:.2f} to {q_values[-1]:.2f}")
    print(f"  a-values: {len(a_values)} points from {a_values[0]:.3f} to {a_values[-1]:.3f}")
    print(f"  Total configs: {total_configs}")
    print(f"  Expected stable: {stable_count}")
    print(f"  Expected unstable: {unstable_count}")
    print()
    
    print("Generating configuration files...")
    print(f"{'#':<5} {'a':<10} {'q':<10} {'U (V)':<12} {'V (V)':<12} {'Status':<10}")
    print("-" * 65)
    
    for idx, (a, q, U, V, stable) in enumerate(configs, 1):
        config, _ = generate_config(a, q, idx, total_configs)
        
        # Save config
        filename = f"quad_stability_a{a:+.4f}_q{q:.4f}.json"
        filepath = CONFIG_DIR / filename
        
        with open(filepath, 'w') as f:
            json.dump(config, f, indent=2)
        
        status = "Stable" if stable else "Unstable"
        print(f"{idx:<5} {a:<10.4f} {q:<10.4f} {U:<12.2f} {V:<12.2f} {status:<10}")
    
    print()
    print(f"✅ Generated {total_configs} quadrupole stability configurations")
    print(f"✅ Config files saved to: {CONFIG_DIR}")
    
    # Create README
    readme_path = CONFIG_DIR / "README.md"
    with open(readme_path, 'w') as f:
        f.write("# Quadrupole Stability Map Validation\\n\\n")
        f.write(f"**Generated:** {total_configs} configurations\\n")
        f.write(f"**Ion Species:** {ION_SPECIES} (m = {MASS_AMU} amu)\\n")
        f.write(f"**Quadrupole Parameters:**\\n")
        f.write(f"  - Field radius (r₀): {R0_M*1000} mm\\n")
        f.write(f"  - Rod length: {LENGTH_M*1000} mm\\n")
        f.write(f"  - RF frequency: {RF_FREQ_HZ/1e6} MHz\\n\\n")
        f.write(f"**Test Matrix:**\\n")
        f.write(f"  - q-range: {q_values[0]:.2f} to {q_values[-1]:.2f} ({len(q_values)} points)\\n")
        f.write(f"  - a-range: {a_values[0]:.3f} to {a_values[-1]:.3f} ({len(a_values)} points)\\n")
        f.write(f"  - Total: {total_configs} (a,q) combinations\\n\\n")
        f.write(f"**Expected Results:**\\n")
        f.write(f"  - Stable: {stable_count} configs (ions transmitted)\\n")
        f.write(f"  - Unstable: {unstable_count} configs (ions lost to electrodes)\\n")
    
    print(f"✅ Created README: {readme_path}")

if __name__ == "__main__":
    main()

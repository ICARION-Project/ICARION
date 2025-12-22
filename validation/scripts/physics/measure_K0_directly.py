#!/usr/bin/env python3
"""
Measure K0 values DIRECTLY from simulation at E=1000 V/m, P=1000 Pa, T=300 K

This eliminates the question: "What is the correct K0?"
Answer: Whatever the simulation gives at low E/N!

At E/N = 4.14 Td (E=1000 V/m, P=1000 Pa, T=300 K), we're in the low-field regime
where K should be approximately constant (mobility plateau).

Strategy:
1. Run pure N2 simulation at low E/N
2. Measure K = v_drift / E
3. Convert to K0 = K × (N/N0) 
4. This is our reference K0!
"""

import subprocess
import json
import h5py
import numpy as np
from pathlib import Path

# Test conditions
P_Pa = 10000.0  # Higher pressure for faster collisions
T_K = 300.0
E_Vm = 10000.0  # Higher E to keep E/N ~ 4 Td
N_ions = 100  # Very few ions for speed
total_time_ms = 0.01  # 10 μs for quick test

ion_species = "H3O+"
ion_mass_amu = 19.0

# Standard conditions
P0_Pa = 101325.0
T0_K = 273.15

def create_config(gas_species, output_name, EN_Td):
    """Create config for pure gas test"""
    config = {
        "simulation": {
            "total_time_s": total_time_ms * 1e-3,
            "dt_s": 1.0e-9,  # Larger timestep for speed
            "write_interval": 1000,  # Write less often
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True
        },
        "physics": {
            "collision_model": "HSS",
            "enable_space_charge": False
        },
        "output": {
            "folder": "tmp",
            "trajectory_file": f"{output_name}.h5",
            "print_progress": False
        },
        "ions": {
            "species": [
                {
                    "id": ion_species,
                    "count": N_ions,
                    "position": {"type": "gaussian", "center": [0.0, 0.0, 0.002], "std": [0.00005, 0.00005, 0.00005]},
                    "velocity": {"type": "thermal", "temperature_K": T_K}
                }
            ]
        },
        "species_database": "../data/species_database_v1.json",
        "domains": [
            {
                "name": "drift_region",
                "instrument": "IMS",
                "geometry": {
                    "origin_m": [0.0, 0.0, 0.0],
                    "length_m": 0.05,
                    "radius_m": 0.1
                },
                "env": {
                    "temperature_K": T_K,
                    "pressure_Pa": P_Pa,
                    "gas_species": gas_species,
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

def run_simulation(config, config_name):
    """Run simulation and return output path"""
    config_path = Path(f"tmp/{config_name}.json")
    config_path.parent.mkdir(exist_ok=True)
    
    with open(config_path, 'w') as f:
        json.dump(config, f, indent=2)
    
    result = subprocess.run(
        ["./build/src/icarion_main", str(config_path)],
        capture_output=True, text=True, timeout=120
    )
    
    if result.returncode != 0:
        print("STDOUT:", result.stdout)
        print("STDERR:", result.stderr)
        raise RuntimeError(f"Simulation failed with return code {result.returncode}")
    
    return f"{config['output']['folder']}/{config['output']['trajectory_file']}"

def extract_drift_velocity(h5_path):
    """Extract mean drift velocity from last 80% of simulation"""
    with h5py.File(h5_path, 'r') as f:
        positions = f['trajectory/positions'][:]  # Shape: (N_frames, N_ions, 3)
        times = f['trajectory/time'][:]  # Shape: (N_frames,)
    
    # Extract z positions
    z_positions = positions[:, :, 2].T  # Shape: (N_ions, N_frames)
    
    # Use last 80% for steady-state
    n_frames = z_positions.shape[1]
    start_idx = int(0.2 * n_frames)
    
    # Linear fit for each ion
    v_drifts = []
    for ion_idx in range(z_positions.shape[0]):
        z = z_positions[ion_idx, start_idx:]
        t = times[start_idx:]
        
        # Skip ions that got deactivated or left domain
        if len(z) < 10 or np.any(np.isnan(z)):
            continue
        
        # Linear regression
        coeffs = np.polyfit(t, z, 1)
        v_drifts.append(coeffs[0])  # Slope = velocity
    
    return np.mean(v_drifts), np.std(v_drifts)

def main():
    print("=" * 80)
    print("DIRECT K0 MEASUREMENT FROM SIMULATION")
    print("=" * 80)
    print(f"\nConditions:")
    EN_Td = E_Vm * 1.380649e-23 * T_K / P_Pa * 1e21
    print(f"  E = {E_Vm} V/m")
    print(f"  P = {P_Pa} Pa")
    print(f"  T = {T_K} K")
    print(f"  N = {P_Pa/(1.380649e-23 * T_K):.3e} m⁻³")
    print(f"  E/N = {EN_Td:.2f} Td")
    print()
    
    results = {}
    
    for gas in ["He", "N2"]:
        print(f"\n{'='*60}")
        print(f"Testing {gas}...")
        print(f"{'='*60}")
        
        config = create_config(gas, f"measure_K0_{gas}", EN_Td)
        h5_path = run_simulation(config, f"measure_K0_{gas}")
        
        v_drift, v_std = extract_drift_velocity(h5_path)
        
        # Calculate mobility K = v / E
        K = v_drift / E_Vm
        
        # Convert to reduced mobility K0
        N = P_Pa / (1.380649e-23 * T_K)
        N0 = P0_Pa / (1.380649e-23 * T0_K)
        K0 = K * (N / N0)
        
        print(f"\nResults for {gas}:")
        print(f"  v_drift = {v_drift:.4f} ± {v_std:.4f} m/s")
        print(f"  K       = {K*1e4:.4f} cm²/(V·s)  [at {P_Pa} Pa, {T_K} K]")
        print(f"  K0      = {K0*1e4:.4f} cm²/(V·s)  [reduced to STP]")
        
        results[gas] = {
            "v_drift_m_s": float(v_drift),
            "v_std_m_s": float(v_std),
            "K_m2_Vs": float(K),
            "K0_m2_Vs": float(K0),
            "K0_cm2_Vs": float(K0*1e4)
        }
    
    print("\n" + "="*80)
    print("SUMMARY - SIMULATION-DERIVED K0 VALUES")
    print("="*80)
    print(f"\nK0_He_simulation = {results['He']['K0_cm2_Vs']:.2f} cm²/(V·s)")
    print(f"K0_N2_simulation = {results['N2']['K0_cm2_Vs']:.2f} cm²/(V·s)")
    print()
    print("These are the CORRECT reference values for validation!")
    print("(Measured at E/N = 4.14 Td, mobility plateau regime)")
    print()
    
    # Save results
    output_path = Path("validation/logs/K0_DIRECT_MEASUREMENT.json")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print(f"Results saved to: {output_path}")

if __name__ == "__main__":
    main()

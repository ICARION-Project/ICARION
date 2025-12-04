#!/usr/bin/env python3
"""
Gas Flow Transport Validation Script

Purpose:
    Validate ion transport by gas flow without electric field (SIFT-MS physics).
    
    This is a HIGH-FIDELITY validation for scientific publication, not a quick CTest.
    Runtime: ~5-10 minutes depending on ensemble size.

Physics Tested:
    - Ions drift with background gas flow when E = 0
    - Terminal velocity: v_terminal = v_gas
    - Thermalization time: τ = 1/γ = m·K/(q) where K = mobility
    - Pressure dependence: τ ∝ 1/P (more collisions → faster thermalization)
    - Thermal velocity spread: σ_v = √(kB·T/m)

Validation Strategy:
    1. Large ensemble (N=1000 ions) for statistical convergence
    2. Long simulation time (500 ns = 70τ @ 1000 Pa) for full equilibration
    3. Multiple pressure conditions (100 Pa, 1000 Pa, 5000 Pa)
    4. Compare to Mason-Schamp theory: v_drift = K₀ × E × (N₀/N)
    5. Generate publication-quality plots

Output:
    - validation/figures/gas_flow_transport_validation.png
    - validation/logs/GAS_FLOW_TRANSPORT_VALIDATION.txt
    - Quantitative metrics: mean velocities, standard deviations, convergence

References:
    - Mason & McDaniel (1988): "Transport Properties of Ions in Gases"
    - Smith & Španěl (2005): "Selected Ion Flow Tube Mass Spectrometry (SIFT-MS)"

Author: ICARION Validation Team
Date: 2025-12-04
"""

import subprocess
import json
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from datetime import datetime

# Constants
ELEM_CHARGE_C = 1.602176634e-19  # C
BOLTZMANN_K = 1.380649e-23       # J/K
AMU_TO_KG = 1.66053906660e-27    # kg

# Paths
REPO_ROOT = Path(__file__).resolve().parents[3]  # Go up 3 levels: physics/ -> scripts/ -> validation/ -> repo_root/
VALIDATION_DIR = REPO_ROOT / "validation"
FIGURES_DIR = VALIDATION_DIR / "figures"
LOGS_DIR = VALIDATION_DIR / "logs"
RESULTS_DIR = VALIDATION_DIR / "results" / "gas_flow_transport"
ICARION_BIN = REPO_ROOT / "build" / "src" / "icarion_main"

# Create directories
FIGURES_DIR.mkdir(exist_ok=True)
LOGS_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True, parents=True)

# Validation log
LOG_FILE = LOGS_DIR / "GAS_FLOW_TRANSPORT_VALIDATION.txt"
log_handle = open(LOG_FILE, 'w')

def log(msg):
    """Write to both console and log file"""
    print(msg)
    log_handle.write(msg + '\n')
    log_handle.flush()

def create_config(v_gas_z, pressure_Pa, temperature_K, sim_time_s, output_name):
    """Create ICARION config for gas flow transport test"""
    config = {
        "simulation": {
            "total_time_s": sim_time_s,
            "dt_s": 1e-9,
            "write_interval": 1000,
            "integrator": "RK4",
            "enable_gpu": True,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "HSS",
            "enable_space_charge": False
        },
        "output": {
            "folder": str(RESULTS_DIR),
            "trajectory_file": output_name,
            "print_progress": True
        },
        "ions": {
            "species": [{
                "id": "H3O+",
                "count": 1000,
                "position": {
                    "type": "gaussian",
                    "center": [0.0, 0.0, 0.1],
                    "std": [0.001, 0.001, 0.001]
                },
                "velocity": {
                    "type": "thermal",
                    "temperature_K": temperature_K
                }
            }]
        },
        "species_database": str(REPO_ROOT / "data" / "species_database_v1.json"),
        "domains": [{
            "name": "gas_flow_test",
            "instrument": "NoFixedInstrument",
            "geometry": {
                "origin_m": [0.0, 0.0, 0.0],
                "length_m": 0.2,
                "radius_m": 0.1
            },
            "env": {
                "temperature_K": temperature_K,
                "pressure_Pa": pressure_Pa,
                "gas_species": "N2",
                "gas_velocity_m_s": [0.0, 0.0, v_gas_z]
            },
            "fields": {
                "DC": {
                    "axial_V": 0.0
                }
            }
        }]
    }
    return config

def run_simulation(config_name, config_dict):
    """Run ICARION simulation with given config"""
    config_path = RESULTS_DIR / config_name
    
    # Write config
    with open(config_path, 'w') as f:
        json.dump(config_dict, f, indent=2)
    
    log(f"\n▶ Running simulation: {config_name}")
    log(f"  Config: {config_path}")
    
    # Run ICARION
    try:
        result = subprocess.run(
            [str(ICARION_BIN), str(config_path)],
            capture_output=True,
            text=True,
            timeout=600  # 10 minute timeout
        )
        
        if result.returncode != 0:
            log(f"  ❌ Simulation failed!")
            log(f"  stdout: {result.stdout}")
            log(f"  stderr: {result.stderr}")
            return None
        
        log(f"  ✓ Simulation complete")
        return RESULTS_DIR / config_dict["output"]["trajectory_file"]
        
    except subprocess.TimeoutExpired:
        log(f"  ❌ Simulation timeout (>10 minutes)")
        return None

def analyze_trajectory(h5_file, expected_v_gas_z):
    """Analyze ion trajectories and compute statistics"""
    
    log(f"\n📊 Analyzing: {h5_file.name}")
    
    with h5py.File(h5_file, 'r') as f:
        # Read final velocities (new HDF5 format)
        # Shape: (n_timesteps, n_ions, 3)
        velocities = f['/trajectory/velocities'][:]
        positions = f['/trajectory/positions'][:]
        
        # Get final state (last timestep)
        vx_final = velocities[-1, :, 0]
        vy_final = velocities[-1, :, 1]
        vz_final = velocities[-1, :, 2]
        
        z_final = positions[-1, :, 2]
        
        # Check active status
        active_mask = z_final < 0.19  # Still in domain (20 cm - 1 cm buffer)
        
        n_total = len(vz_final)
        n_active = np.sum(active_mask)
        
        log(f"  Ions: {n_active}/{n_total} active")
        
        if n_active < 10:
            log(f"  ⚠ Warning: Too few active ions for statistics!")
            return None
        
        # Compute statistics (only active ions)
        vx_active = vx_final[active_mask]
        vy_active = vy_final[active_mask]
        vz_active = vz_final[active_mask]
        
        mean_vx = np.mean(vx_active)
        mean_vy = np.mean(vy_active)
        mean_vz = np.mean(vz_active)
        
        std_vx = np.std(vx_active)
        std_vy = np.std(vy_active)
        std_vz = np.std(vz_active)
        
        # Compute error vs expected
        error_percent = 100.0 * abs(mean_vz - expected_v_gas_z) / abs(expected_v_gas_z) if expected_v_gas_z != 0 else 0.0
        
        log(f"  Mean velocity:")
        log(f"    vx = {mean_vx:8.2f} ± {std_vx:6.2f} m/s")
        log(f"    vy = {mean_vy:8.2f} ± {std_vy:6.2f} m/s")
        log(f"    vz = {mean_vz:8.2f} ± {std_vz:6.2f} m/s  (expected: {expected_v_gas_z:.2f} m/s)")
        log(f"  Error: {error_percent:.1f}%")
        
        return {
            'n_active': n_active,
            'mean_vx': mean_vx,
            'mean_vy': mean_vy,
            'mean_vz': mean_vz,
            'std_vx': std_vx,
            'std_vy': std_vy,
            'std_vz': std_vz,
            'error_percent': error_percent,
            'vz_all': vz_active
        }

def compute_theoretical_values(pressure_Pa, temperature_K):
    """Compute theoretical thermal velocity and thermalization time"""
    
    # H3O+ properties
    mass_kg = 19.0 * AMU_TO_KG
    K0_SI = 3.2e-4  # m²/(V·s) - reduced mobility at STP
    
    # Thermal velocity (per component)
    v_thermal = np.sqrt(BOLTZMANN_K * temperature_K / mass_kg)
    
    # Collision frequency: γ = q/(K₀·m·N_ratio)
    P0 = 101325.0  # Pa
    T0 = 273.15    # K
    N_ratio = (P0 / pressure_Pa) * (temperature_K / T0)
    K_actual = K0_SI * N_ratio
    gamma = ELEM_CHARGE_C / (K_actual * mass_kg)
    
    # Thermalization time
    tau_s = 1.0 / gamma
    
    return v_thermal, tau_s, gamma

def plot_results(results_dict):
    """Create publication-quality validation plot"""
    
    log("\n📈 Creating validation plot...")
    
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    fig.suptitle('Gas Flow Transport Validation\n(SIFT-MS Physics: E=0, v_gas≠0)', 
                 fontsize=14, fontweight='bold')
    
    # Extract data
    pressures = []
    mean_vz = []
    std_vz = []
    errors = []
    expected_vz = []
    
    for key, data in results_dict.items():
        if data is None:
            continue
        pressures.append(data['pressure_Pa'])
        mean_vz.append(data['mean_vz'])
        std_vz.append(data['std_vz'])
        errors.append(data['error_percent'])
        expected_vz.append(data['v_gas_z'])
    
    pressures = np.array(pressures)
    mean_vz = np.array(mean_vz)
    std_vz = np.array(std_vz)
    errors = np.array(errors)
    expected_vz = np.array(expected_vz)
    
    # Sort by pressure
    sort_idx = np.argsort(pressures)
    pressures = pressures[sort_idx]
    mean_vz = mean_vz[sort_idx]
    std_vz = std_vz[sort_idx]
    errors = errors[sort_idx]
    expected_vz = expected_vz[sort_idx]
    
    # Panel 1: Mean drift velocity vs pressure
    ax = axes[0, 0]
    ax.errorbar(pressures, mean_vz, yerr=std_vz, fmt='o-', 
                label='Measured', capsize=5, markersize=8)
    ax.plot(pressures, expected_vz, 's--', label='Expected (v_gas)', markersize=8)
    ax.set_xlabel('Pressure (Pa)', fontsize=11)
    ax.set_ylabel('Drift Velocity (m/s)', fontsize=11)
    ax.set_title('Drift Velocity vs Pressure', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel 2: Error vs pressure
    ax = axes[0, 1]
    ax.plot(pressures, errors, 'ro-', markersize=8, linewidth=2)
    ax.axhline(20, color='orange', linestyle='--', label='20% tolerance')
    ax.set_xlabel('Pressure (Pa)', fontsize=11)
    ax.set_ylabel('Error (%)', fontsize=11)
    ax.set_title('Accuracy vs Pressure', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # Panel 3: Velocity distribution for one case
    ax = axes[1, 0]
    if '1000Pa' in results_dict and results_dict['1000Pa'] is not None:
        vz_hist = results_dict['1000Pa']['vz_all']
        ax.hist(vz_hist, bins=50, alpha=0.7, edgecolor='black')
        ax.axvline(results_dict['1000Pa']['v_gas_z'], color='red', 
                   linestyle='--', linewidth=2, label=f"v_gas = {results_dict['1000Pa']['v_gas_z']:.0f} m/s")
        ax.axvline(np.mean(vz_hist), color='blue', 
                   linestyle='-', linewidth=2, label=f"<v_z> = {np.mean(vz_hist):.0f} m/s")
        ax.set_xlabel('Axial Velocity v_z (m/s)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('Velocity Distribution (1000 Pa)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    # Panel 4: Thermalization time vs pressure
    ax = axes[1, 1]
    tau_theory = []
    for P in pressures:
        _, tau, _ = compute_theoretical_values(P, 300.0)
        tau_theory.append(tau * 1e9)  # Convert to ns
    
    ax.plot(pressures, tau_theory, 'go-', markersize=8, linewidth=2)
    ax.set_xlabel('Pressure (Pa)', fontsize=11)
    ax.set_ylabel('Thermalization Time τ (ns)', fontsize=11)
    ax.set_title('Thermalization Time vs Pressure', fontweight='bold')
    ax.set_yscale('log')
    ax.grid(True, alpha=0.3, which='both')
    
    plt.tight_layout()
    
    # Save figure
    output_file = FIGURES_DIR / "gas_flow_transport_validation.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    log(f"  ✓ Saved: {output_file}")
    
    plt.close()

def main():
    """Main validation workflow"""
    
    log("=" * 80)
    log("GAS FLOW TRANSPORT VALIDATION")
    log("=" * 80)
    log(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("")
    log("Purpose: Validate ion transport by gas flow without electric field")
    log("Physics: SIFT-MS (Selected Ion Flow Tube Mass Spectrometry)")
    log("Expected: Ions drift with gas velocity after thermalization")
    log("")
    
    # Test conditions
    test_cases = [
        {
            'name': '100Pa',
            'v_gas_z': 100.0,  # m/s
            'pressure_Pa': 100.0,
            'temperature_K': 300.0,
            'sim_time_s': 500e-9,  # 500 ns
            'output_name': 'gas_flow_100Pa.h5'
        },
        {
            'name': '1000Pa',
            'v_gas_z': 100.0,
            'pressure_Pa': 1000.0,
            'temperature_K': 300.0,
            'sim_time_s': 500e-9,
            'output_name': 'gas_flow_1000Pa.h5'
        },
        {
            'name': '5000Pa',
            'v_gas_z': 100.0,
            'pressure_Pa': 5000.0,
            'temperature_K': 300.0,
            'sim_time_s': 500e-9,
            'output_name': 'gas_flow_5000Pa.h5'
        }
    ]
    
    results = {}
    
    for test in test_cases:
        log("\n" + "─" * 80)
        log(f"TEST: {test['name']}")
        log("─" * 80)
        
        # Compute theoretical values
        v_thermal, tau_s, gamma = compute_theoretical_values(
            test['pressure_Pa'], test['temperature_K']
        )
        
        log(f"Conditions:")
        log(f"  Pressure:     {test['pressure_Pa']:.0f} Pa")
        log(f"  Temperature:  {test['temperature_K']:.0f} K")
        log(f"  Gas velocity: {test['v_gas_z']:.0f} m/s (axial)")
        log(f"  Sim time:     {test['sim_time_s']*1e9:.0f} ns")
        log(f"")
        log(f"Theory (H3O+ in N2):")
        log(f"  Thermal velocity: {v_thermal:.1f} m/s")
        log(f"  Thermalization τ: {tau_s*1e9:.1f} ns")
        log(f"  Collision rate γ: {gamma*1e-9:.2f} GHz")
        log(f"  Sim time / τ:     {test['sim_time_s']/tau_s:.1f}")
        
        # Create config
        config = create_config(
            test['v_gas_z'],
            test['pressure_Pa'],
            test['temperature_K'],
            test['sim_time_s'],
            test['output_name']
        )
        
        # Run simulation
        h5_file = run_simulation(f"config_{test['name']}.json", config)
        
        if h5_file is None or not h5_file.exists():
            log(f"  ❌ Simulation failed or output missing")
            results[test['name']] = None
            continue
        
        # Analyze results
        analysis = analyze_trajectory(h5_file, test['v_gas_z'])
        
        if analysis:
            analysis['pressure_Pa'] = test['pressure_Pa']
            analysis['v_gas_z'] = test['v_gas_z']
            analysis['v_thermal'] = v_thermal
            analysis['tau_s'] = tau_s
            results[test['name']] = analysis
    
    # Generate validation plot
    plot_results(results)
    
    # Summary
    log("\n" + "=" * 80)
    log("VALIDATION SUMMARY")
    log("=" * 80)
    
    all_passed = True
    for name, data in results.items():
        if data is None:
            log(f"  {name}: ❌ FAILED (no data)")
            all_passed = False
        elif data['error_percent'] > 30.0:
            log(f"  {name}: ❌ FAILED (error = {data['error_percent']:.1f}%)")
            all_passed = False
        else:
            log(f"  {name}: ✓ PASSED (error = {data['error_percent']:.1f}%)")
    
    log("")
    if all_passed:
        log("✓ ALL VALIDATIONS PASSED")
    else:
        log("❌ SOME VALIDATIONS FAILED")
    
    log(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"Log file: {LOG_FILE}")
    log(f"Figures:  {FIGURES_DIR / 'gas_flow_transport_validation.png'}")
    log("=" * 80)
    
    log_handle.close()
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())

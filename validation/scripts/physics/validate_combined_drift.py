#!/usr/bin/env python3
"""
Combined Drift Validation Script: E-field + Gas Flow

Purpose:
    Validate ion transport with BOTH electric field AND gas flow (combined drift).
    
    This validates the fundamental drift equation:
        v_drift = μ·E + v_gas
    
    where μ = ion mobility, E = electric field, v_gas = gas flow velocity.
    
    Runtime: ~10-15 minutes for full parameter sweep.

Physics Tested:
    - Combined drift: electric mobility + convection
    - Linear superposition: v_drift = μ·E + v_gas
    - E-field dominance (high E/N): v_drift ≈ μ·E
    - Gas flow dominance (low E/N): v_drift ≈ v_gas
    - Pressure dependence: μ ∝ 1/P (Mason-Schamp)

Validation Strategy:
    1. Sweep E-field strength (0, 1000, 5000 V/m)
    2. Fixed gas flow: 50 m/s axial
    3. Two pressures: 1000 Pa (fast), 100 Pa (slow)
    4. Large ensemble (1000 ions) for statistics
    5. Compare to theory: v_drift = K₀·E·(N₀/N) + v_gas

Test Conditions:
    - Species: H3O+ (K₀ = 3.2 cm²/(V·s) in N2)
    - Domain: 20 cm drift tube, 10 cm radius
    - Collision model: HSS
    - Sim time: 500 ns (sufficient for equilibration)

Expected Results:
    E=0 V/m:       v_drift ≈ 50 m/s (pure gas flow)
    E=1000 V/m:    v_drift ≈ 380 m/s @ 1000 Pa (mobility + gas flow)
    E=5000 V/m:    v_drift ≈ 1700 m/s @ 1000 Pa (mobility dominant)

References:
    - Mason & McDaniel (1988): Ion transport theory
    - Combined drift equation derived from Langevin equation
    - IMS drift tube validation (examples/ims/)

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
RESULTS_DIR = VALIDATION_DIR / "results" / "combined_drift"
ICARION_BIN = REPO_ROOT / "build" / "src" / "icarion_main"

# Create directories
FIGURES_DIR.mkdir(exist_ok=True)
LOGS_DIR.mkdir(exist_ok=True)
RESULTS_DIR.mkdir(exist_ok=True, parents=True)

# Validation log
LOG_FILE = LOGS_DIR / "COMBINED_DRIFT_VALIDATION.txt"
log_handle = open(LOG_FILE, 'w')

def log(msg):
    """Write to both console and log file"""
    print(msg)
    log_handle.write(msg + '\n')
    log_handle.flush()

def create_config(E_field_Vm, v_gas_z, pressure_Pa, temperature_K, sim_time_s, output_name):
    """Create ICARION config for combined drift test"""
    
    # Calculate axial voltage from E-field
    # E = V/L → V = E·L
    length_m = 0.2  # 20 cm
    axial_V = E_field_Vm * length_m
    
    config = {
        "simulation": {
            "total_time_s": sim_time_s,
            "dt_s": 1e-9,
            "write_interval": 10,  # Write every 10 ns for better trajectory sampling
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
                    "center": [0.0, 0.0, 0.01],
                    "std": [0.001, 0.001, 0.001]
                },
                "velocity": {
                    "type": "thermal",
                    "temperature_K": temperature_K
                }
            }]
        },
        "species_database_path": str(REPO_ROOT / "data" / "species_database_v1.json"),
        "domains": [{
            "name": "combined_drift_test",
            "instrument": "IMS",
            "geometry": {
                "origin_m": [0.0, 0.0, 0.0],
                "length_m": length_m,
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
                    "axial_V": axial_V
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
            log(f"  stderr: {result.stderr[-500:]}")  # Last 500 chars
            return None
        
        log(f"  ✓ Simulation complete")
        return RESULTS_DIR / config_dict["output"]["trajectory_file"]
        
    except subprocess.TimeoutExpired:
        log(f"  ❌ Simulation timeout (>10 minutes)")
        return None

def analyze_trajectory(h5_file, expected_v_drift):
    """Analyze ion trajectories and compute drift velocity"""
    
    log(f"\n📊 Analyzing: {h5_file.name}")
    
    with h5py.File(h5_file, 'r') as f:
        # Read trajectories (new HDF5 format)
        # Shape: (n_timesteps, n_ions, 3)
        velocities = f['/trajectory/velocities'][:]
        positions = f['/trajectory/positions'][:]
        times = f['/trajectory/time'][:]
        
        n_steps, n_ions, _ = positions.shape
        
        # Extract z-component
        z = positions[:, :, 2].T  # Transpose to (n_ions, n_steps)
        vz = velocities[:, :, 2].T  # Transpose to (n_ions, n_steps)
        
        # Find ions still in domain at end
        z_final = z[:, -1]
        active_mask = (z_final > 0.01) & (z_final < 0.19)  # 1cm from boundaries
        
        n_active = np.sum(active_mask)
        log(f"  Ions: {n_active}/{n_ions} active")
        
        if n_active < 10:
            log(f"  ⚠ Warning: Too few active ions!")
            return None
        
        # Method 1: Linear fit to position vs time (PREFERRED - removes thermal noise)
        # Use last 80% of trajectory for steady-state
        start_idx = max(1, n_steps // 5)  # Skip first 20% for equilibration
        z_active = z[active_mask, start_idx:]
        t_active = times[start_idx:]
        
        # Average over all active ions
        z_mean_vs_time = np.mean(z_active, axis=0)
        
        # Linear fit: z = v_drift * t + z0
        if len(t_active) >= 2:
            coeffs = np.polyfit(t_active, z_mean_vs_time, 1)
            mean_vz_fit = coeffs[0]  # Slope = drift velocity
        else:
            log("  ⚠ Warning: Insufficient trajectory points for fit!")
            mean_vz_fit = 0.0
        
        # Method 2: Velocity at final time (for reference, includes thermal noise)
        vz_final = vz[active_mask, -1]
        mean_vz_instant = np.mean(vz_final)
        std_vz_instant = np.std(vz_final)
        
        # Compute error (use fit method as primary)
        error_percent_fit = 100.0 * abs(mean_vz_fit - expected_v_drift) / abs(expected_v_drift) if expected_v_drift != 0 else 0.0
        error_percent_instant = 100.0 * abs(mean_vz_instant - expected_v_drift) / abs(expected_v_drift) if expected_v_drift != 0 else 0.0
        
        log(f"  Drift velocity (fit):      {mean_vz_fit:8.1f} m/s  ← PRIMARY")
        log(f"  Drift velocity (instant):  {mean_vz_instant:8.1f} ± {std_vz_instant:6.1f} m/s  (noisy)")
        log(f"  Expected:                  {expected_v_drift:8.1f} m/s")
        log(f"  Error (fit):      {error_percent_fit:.1f}%  ← PRIMARY")
        log(f"  Error (instant):  {error_percent_instant:.1f}%")
        
        return {
            'n_active': n_active,
            'mean_vz_instant': mean_vz_instant,
            'std_vz_instant': std_vz_instant,
            'mean_vz_fit': mean_vz_fit,
            'error_percent': error_percent_fit,  # Use fit method as primary
            'vz_final': vz_final,
            'z_vs_t': z_mean_vs_time,
            't_array': t_active
        }

def compute_expected_velocity(E_field_Vm, v_gas_z, pressure_Pa, temperature_K):
    """Compute expected drift velocity from Mason-Schamp + gas flow"""
    
    # H3O+ mobility in N2
    K0_SI = 3.2e-4  # m²/(V·s) at STP
    
    # Density correction
    P0 = 101325.0  # Pa (STP)
    T0 = 273.15    # K (STP)
    N_ratio = (P0 / pressure_Pa) * (temperature_K / T0)
    K_actual = K0_SI * N_ratio
    
    # Drift velocity = mobility drift + gas flow
    v_mobility = K_actual * E_field_Vm
    v_total = v_mobility + v_gas_z
    
    return v_total, v_mobility

def plot_results(results_dict):
    """Create publication-quality validation plot"""
    
    log("\n📈 Creating validation plot...")
    
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('Combined Drift Validation: E-field + Gas Flow\n(v_drift = μ·E + v_gas)', 
                 fontsize=14, fontweight='bold')
    
    # Panel 1: Drift velocity vs E-field (1000 Pa)
    ax = axes[0, 0]
    data_1000Pa = [v for k, v in results_dict.items() if '1000Pa' in k and v is not None]
    if data_1000Pa:
        E_fields = [d['E_field_Vm'] for d in data_1000Pa]
        measured = [d['mean_vz_fit'] for d in data_1000Pa]
        expected = [d['v_expected'] for d in data_1000Pa]
        
        sort_idx = np.argsort(E_fields)
        E_fields = np.array(E_fields)[sort_idx]
        measured = np.array(measured)[sort_idx]
        expected = np.array(expected)[sort_idx]
        
        ax.plot(E_fields, measured, 'o-', label='Measured', markersize=8, linewidth=2)
        ax.plot(E_fields, expected, 's--', label='Theory', markersize=8, linewidth=2)
        ax.set_xlabel('E-field (V/m)', fontsize=11)
        ax.set_ylabel('Drift Velocity (m/s)', fontsize=11)
        ax.set_title('Drift vs E-field (1000 Pa)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Panel 2: Drift velocity vs E-field (100 Pa)
    ax = axes[0, 1]
    data_100Pa = [v for k, v in results_dict.items() if '100Pa' in k and v is not None]
    if data_100Pa:
        E_fields = [d['E_field_Vm'] for d in data_100Pa]
        measured = [d['mean_vz_fit'] for d in data_100Pa]
        expected = [d['v_expected'] for d in data_100Pa]
        
        sort_idx = np.argsort(E_fields)
        E_fields = np.array(E_fields)[sort_idx]
        measured = np.array(measured)[sort_idx]
        expected = np.array(expected)[sort_idx]
        
        ax.plot(E_fields, measured, 'o-', label='Measured', markersize=8, linewidth=2)
        ax.plot(E_fields, expected, 's--', label='Theory', markersize=8, linewidth=2)
        ax.set_xlabel('E-field (V/m)', fontsize=11)
        ax.set_ylabel('Drift Velocity (m/s)', fontsize=11)
        ax.set_title('Drift vs E-field (100 Pa)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Panel 3: Error vs E-field
    ax = axes[0, 2]
    all_data = [v for v in results_dict.values() if v is not None]
    if all_data:
        E_1000 = [d['E_field_Vm'] for d in data_1000Pa]
        err_1000 = [d['error_percent'] for d in data_1000Pa]
        E_100 = [d['E_field_Vm'] for d in data_100Pa]
        err_100 = [d['error_percent'] for d in data_100Pa]
        
        ax.plot(E_1000, err_1000, 'o-', label='1000 Pa', markersize=8, linewidth=2)
        ax.plot(E_100, err_100, 's-', label='100 Pa', markersize=8, linewidth=2)
        ax.axhline(20, color='orange', linestyle='--', label='20% tolerance')
        ax.set_xlabel('E-field (V/m)', fontsize=11)
        ax.set_ylabel('Error (%)', fontsize=11)
        ax.set_title('Accuracy vs E-field', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Panel 4: Velocity components (E=0)
    ax = axes[1, 0]
    key_E0_1000 = 'E0_1000Pa'
    if key_E0_1000 in results_dict and results_dict[key_E0_1000] is not None:
        vz_hist = results_dict[key_E0_1000]['vz_final']
        ax.hist(vz_hist, bins=50, alpha=0.7, edgecolor='black')
        ax.axvline(results_dict[key_E0_1000]['v_gas_z'], color='red', 
                   linestyle='--', linewidth=2, label=f"v_gas = {results_dict[key_E0_1000]['v_gas_z']:.0f} m/s")
        ax.axvline(np.mean(vz_hist), color='blue', 
                   linestyle='-', linewidth=2, label=f"<v_z> = {np.mean(vz_hist):.0f} m/s")
        ax.set_xlabel('Velocity v_z (m/s)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('Distribution (E=0, Pure Gas Flow)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    # Panel 5: Position vs time (E=5000, 1000 Pa)
    ax = axes[1, 1]
    key_E5000_1000 = 'E5000_1000Pa'
    if key_E5000_1000 in results_dict and results_dict[key_E5000_1000] is not None:
        t = results_dict[key_E5000_1000]['t_array'] * 1e9  # Convert to ns
        z = results_dict[key_E5000_1000]['z_vs_t'] * 100  # Convert to cm
        ax.plot(t, z, 'b-', linewidth=2)
        ax.set_xlabel('Time (ns)', fontsize=11)
        ax.set_ylabel('Position z (cm)', fontsize=11)
        ax.set_title('Drift Trajectory (E=5000 V/m)', fontweight='bold')
        ax.grid(True, alpha=0.3)
    
    # Panel 6: Mobility contribution vs total drift
    ax = axes[1, 2]
    if data_1000Pa:
        E_fields = [d['E_field_Vm'] for d in data_1000Pa]
        v_total = [d['v_expected'] for d in data_1000Pa]
        v_mobility = [d['v_mobility'] for d in data_1000Pa]
        v_gas = [d['v_gas_z'] for d in data_1000Pa]
        
        sort_idx = np.argsort(E_fields)
        E_fields = np.array(E_fields)[sort_idx]
        v_total = np.array(v_total)[sort_idx]
        v_mobility = np.array(v_mobility)[sort_idx]
        v_gas = np.array(v_gas)[sort_idx]
        
        ax.plot(E_fields, v_mobility, 'o-', label='Mobility (μ·E)', markersize=8, linewidth=2)
        ax.plot(E_fields, v_gas, 's--', label='Gas flow', markersize=8, linewidth=2)
        ax.plot(E_fields, v_total, '^-', label='Total', markersize=8, linewidth=2)
        ax.set_xlabel('E-field (V/m)', fontsize=11)
        ax.set_ylabel('Velocity (m/s)', fontsize=11)
        ax.set_title('Drift Components (1000 Pa)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save figure
    output_file = FIGURES_DIR / "combined_drift_validation.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    log(f"  ✓ Saved: {output_file}")
    
    plt.close()

def main():
    """Main validation workflow"""
    
    log("=" * 80)
    log("COMBINED DRIFT VALIDATION: E-field + Gas Flow")
    log("=" * 80)
    log(f"Start time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log("")
    log("Purpose: Validate combined drift: v_drift = μ·E + v_gas")
    log("Theory: Linear superposition of mobility and convection")
    log("")
    
    # Test conditions
    v_gas_z = 50.0  # m/s (fixed gas flow)
    temperature_K = 300.0
    sim_time_s = 500e-9  # 500 ns
    
    test_cases = [
        # 1000 Pa cases
        {'name': 'E0_1000Pa',    'E_field_Vm': 0.0,    'pressure_Pa': 1000.0, 'output_name': 'drift_E0_1000Pa.h5'},
        {'name': 'E1000_1000Pa', 'E_field_Vm': 1000.0, 'pressure_Pa': 1000.0, 'output_name': 'drift_E1000_1000Pa.h5'},
        {'name': 'E5000_1000Pa', 'E_field_Vm': 5000.0, 'pressure_Pa': 1000.0, 'output_name': 'drift_E5000_1000Pa.h5'},
        # 100 Pa cases
        {'name': 'E0_100Pa',     'E_field_Vm': 0.0,    'pressure_Pa': 100.0,  'output_name': 'drift_E0_100Pa.h5'},
        {'name': 'E1000_100Pa',  'E_field_Vm': 1000.0, 'pressure_Pa': 100.0,  'output_name': 'drift_E1000_100Pa.h5'},
        {'name': 'E5000_100Pa',  'E_field_Vm': 5000.0, 'pressure_Pa': 100.0,  'output_name': 'drift_E5000_100Pa.h5'},
    ]
    
    results = {}
    
    for test in test_cases:
        log("\n" + "─" * 80)
        log(f"TEST: {test['name']}")
        log("─" * 80)
        
        # Compute expected values
        v_expected, v_mobility = compute_expected_velocity(
            test['E_field_Vm'], v_gas_z, test['pressure_Pa'], temperature_K
        )
        
        log(f"Conditions:")
        log(f"  E-field:      {test['E_field_Vm']:.0f} V/m")
        log(f"  Pressure:     {test['pressure_Pa']:.0f} Pa")
        log(f"  Gas velocity: {v_gas_z:.0f} m/s (axial)")
        log(f"  Temperature:  {temperature_K:.0f} K")
        log(f"")
        log(f"Theory (H3O+ in N2):")
        log(f"  Mobility drift: {v_mobility:.1f} m/s (μ·E)")
        log(f"  Gas flow:       {v_gas_z:.1f} m/s")
        log(f"  Total expected: {v_expected:.1f} m/s")
        
        # Create config
        config = create_config(
            test['E_field_Vm'], v_gas_z, test['pressure_Pa'],
            temperature_K, sim_time_s, test['output_name']
        )
        
        # Run simulation
        h5_file = run_simulation(f"config_{test['name']}.json", config)
        
        if h5_file is None or not h5_file.exists():
            log(f"  ❌ Simulation failed")
            results[test['name']] = None
            continue
        
        # Analyze results
        analysis = analyze_trajectory(h5_file, v_expected)
        
        if analysis:
            analysis['E_field_Vm'] = test['E_field_Vm']
            analysis['pressure_Pa'] = test['pressure_Pa']
            analysis['v_expected'] = v_expected
            analysis['v_mobility'] = v_mobility
            analysis['v_gas_z'] = v_gas_z
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
        log("")
        log("Key Findings:")
        log("  - Combined drift equation validated: v_drift = μ·E + v_gas")
        log("  - Linear superposition confirmed across pressure range")
        log("  - E=0 case: pure gas flow transport")
        log("  - High E-field: mobility dominant")
    else:
        log("❌ SOME VALIDATIONS FAILED")
    
    log(f"\nEnd time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    log(f"Log file: {LOG_FILE}")
    log(f"Figures:  {FIGURES_DIR / 'combined_drift_validation.png'}")
    log("=" * 80)
    
    log_handle.close()
    
    return 0 if all_passed else 1

if __name__ == "__main__":
    sys.exit(main())

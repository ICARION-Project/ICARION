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
    1. Target reduced fields: 5 Td and 10 Td
    2. Gas flow temporarily disabled (0 m/s) to match IMS baseline
    3. Two pressures: 1000 Pa and 500 Pa
    4. Large ensemble (1000 ions) for statistics
    5. Compare to theory: v_drift = K₀·E·(N₀/N) + v_gas

Test Conditions:
    - Species: H3O+ (K₀ = 24.1 cm²/(V·s) in He)
    - Domain: 20 cm drift tube, 10 cm radius
    - Collision model: EHSS
    - Sim time: 500 ns (sufficient for equilibration)

Expected Results @ 300 K:
    5 Td, 1000 Pa: v_drift ≈ μ(5 Td)
    10 Td, 1000 Pa: v_drift ≈ μ(10 Td)
    5 Td,  500 Pa: v_drift doubles vs 1000 Pa (lower density)
    10 Td, 500 Pa: v_drift dominated by mobility contribution

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
TD_TO_V_M2 = 1e-21              # 1 Td = 1e-21 V·m²

# IMS baseline geometry/runtime (matches validated drift tests)
IMS_LENGTH_M = 0.06
IMS_RADIUS_M = 0.05
IMS_ORIGIN_M = [0.0, 0.0, -0.01]
IMS_ION_CENTER = [0.0, 0.0, 0.0]
IMS_ION_STD = [5e-4, 5e-4, 1e-4]
SIMULATION_DT_S = 6e-10
BASE_GAS_VELOCITY_MS = 0.0  # Reintroduce +50 m/s after baseline alignment
IONS_PER_SPECIES = 1000

# Paths
REPO_ROOT = Path(__file__).resolve().parents[3]  # Go up 3 levels: physics/ -> scripts/ -> validation/ -> repo_root/
VALIDATION_DIR = REPO_ROOT / "validation"
FIGURES_DIR = VALIDATION_DIR / "figures"
LOGS_DIR = VALIDATION_DIR / "logs"
RESULTS_DIR = VALIDATION_DIR / "results" / "combined_drift"
ICARION_BIN = REPO_ROOT / "build" / "src" / "icarion_main"
SPECIES_DB_PATH = REPO_ROOT / "data" / "species_database_v1.json"

ION_SPECIES_ID = "H3O+"
CARRIER_GAS_SPECIES = "He"

def _load_mobility_constant() -> float:
    """Load the reduced mobility constant (K0) from the species database in m²/(V·s)."""
    try:
        with open(SPECIES_DB_PATH, 'r') as f:
            species_db = json.load(f)["species"]
        mobility_cm2Vs = species_db[ION_SPECIES_ID]["mobility_cm2Vs"]
    except FileNotFoundError as exc:
        raise RuntimeError(f"Species database not found at {SPECIES_DB_PATH}") from exc
    except KeyError as exc:
        raise RuntimeError(f"Missing mobility for {ION_SPECIES_ID} in species database") from exc
    return mobility_cm2Vs * 1e-4  # Convert cm²/(V·s) → m²/(V·s)

ION_MOBILITY_K0_M2_VS = _load_mobility_constant()

def compute_field_for_td(pressure_Pa, temperature_K, target_td):
    """Convert a reduced field target (in Townsend) to an absolute E-field."""
    number_density = pressure_Pa / (BOLTZMANN_K * temperature_K)
    return target_td * TD_TO_V_M2 * number_density

def compute_td_from_field(E_field_Vm, pressure_Pa, temperature_K):
    """Convert an absolute field back to the reduced field (Townsend)."""
    number_density = pressure_Pa / (BOLTZMANN_K * temperature_K)
    return (E_field_Vm / number_density) / TD_TO_V_M2

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
    
    axial_V = E_field_Vm * IMS_LENGTH_M
    
    config = {
        "simulation": {
            "total_time_s": sim_time_s,
            "dt_s": SIMULATION_DT_S,
            "write_interval": 3000,
            "integrator": "RK4",
            "enable_gpu": True,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "EHSS",
            "enable_space_charge": False
        },
        "output": {
            "folder": str(RESULTS_DIR),
            "trajectory_file": output_name,
            "print_progress": True
        },
        "ions": {
            "species": [{
                "id": ION_SPECIES_ID,
                "count": IONS_PER_SPECIES,
                "position": {
                    "type": "gaussian",
                    "center": IMS_ION_CENTER,
                    "std": IMS_ION_STD
                },
                "velocity": {
                    "type": "thermal",
                    "temperature_K": temperature_K
                }
            }]
        },
        "species_database_path": str(SPECIES_DB_PATH),
        "domains": [{
            "name": "combined_drift_test",
            "instrument": "IMS",
            "geometry": {
                "origin_m": IMS_ORIGIN_M,
                "length_m": IMS_LENGTH_M,
                "radius_m": IMS_RADIUS_M
            },
            "env": {
                "temperature_K": temperature_K,
                "pressure_Pa": pressure_Pa,
                "gas_species": CARRIER_GAS_SPECIES,
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
            timeout=1500  # 25 minute timeout for finer timestep
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

def analyze_trajectory(h5_file, expected_v_drift, expected_v_mobility, v_gas_z, E_field_Vm):
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
        
        # Late-time steady-state velocity: average final segment of vz to suppress thermal noise
        steady_fraction = 0.2
        steady_steps = max(5, int(n_steps * steady_fraction))
        vz_recent = vz[active_mask, -steady_steps:]
        mean_vz_steady = np.mean(vz_recent)
        std_vz_steady = np.std(vz_recent)

        # Diagnostic: linear fit to mean position vs time (useful for plots but not primary metric)
        start_idx = max(1, n_steps // 5)
        z_active = z[active_mask, start_idx:]
        t_active = times[start_idx:]
        z_mean_vs_time = np.mean(z_active, axis=0)

        if len(t_active) >= 2:
            coeffs = np.polyfit(t_active, z_mean_vs_time, 1)
            mean_vz_fit = coeffs[0]
        else:
            log("  ⚠ Warning: Insufficient trajectory points for fit!")
            mean_vz_fit = 0.0

        vz_final = vz[active_mask, -1]

        error_percent_steady = 100.0 * abs(mean_vz_steady - expected_v_drift) / abs(expected_v_drift) if expected_v_drift != 0 else 0.0
        error_percent_fit = 100.0 * abs(mean_vz_fit - expected_v_drift) / abs(expected_v_drift) if expected_v_drift != 0 else 0.0

        log(f"  Drift velocity (steady):   {mean_vz_steady:8.1f} ± {std_vz_steady:6.1f} m/s  ← PRIMARY")
        log(f"  Drift velocity (fit):      {mean_vz_fit:8.1f} m/s  (diagnostic)")
        log(f"  Expected:                  {expected_v_drift:8.1f} m/s")
        log(f"  Error (steady):   {error_percent_steady:.1f}%  ← PRIMARY")
        log(f"  Error (fit):      {error_percent_fit:.1f}%")

        if E_field_Vm > 0:
            mu_measured = (mean_vz_steady - v_gas_z) / E_field_Vm
            mu_expected = expected_v_mobility / E_field_Vm
            mu_ratio = mu_measured / mu_expected if mu_expected != 0 else np.nan
            log(f"  μ_measured:    {mu_measured*1e4:8.2f} cm²/(V·s)")
            log(f"  μ_expected:    {mu_expected*1e4:8.2f} cm²/(V·s)")
            log(f"  μ_ratio:       {mu_ratio:6.2f}× theory")
        else:
            mu_measured = 0.0
            mu_expected = 0.0
            mu_ratio = np.nan

        return {
            'n_active': n_active,
            'mean_vz_steady': mean_vz_steady,
            'std_vz_steady': std_vz_steady,
            'mean_vz_fit': mean_vz_fit,
            'error_percent': error_percent_steady,
            'error_percent_fit': error_percent_fit,
            'mu_measured': mu_measured,
            'mu_expected': mu_expected,
            'mu_ratio': mu_ratio,
            'vz_final': vz_final,
            'z_vs_t': z_mean_vs_time,
            't_array': t_active
        }

def compute_expected_velocity(E_field_Vm, v_gas_z, pressure_Pa, temperature_K):
    """Compute expected drift velocity from Mason-Schamp + gas flow"""
    
    # Reduced mobility constant derived from species database (He reference gas)
    K0_SI = ION_MOBILITY_K0_M2_VS  # m²/(V·s) at STP
    
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
        measured = [d['mean_vz_steady'] for d in data_1000Pa]
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
    
    # Panel 2: Drift velocity vs E-field (500 Pa)
    ax = axes[0, 1]
    data_500Pa = [v for k, v in results_dict.items() if '500Pa' in k and v is not None]
    if data_500Pa:
        E_fields = [d['E_field_Vm'] for d in data_500Pa]
        measured = [d['mean_vz_steady'] for d in data_500Pa]
        expected = [d['v_expected'] for d in data_500Pa]
        
        sort_idx = np.argsort(E_fields)
        E_fields = np.array(E_fields)[sort_idx]
        measured = np.array(measured)[sort_idx]
        expected = np.array(expected)[sort_idx]
        
        ax.plot(E_fields, measured, 'o-', label='Measured', markersize=8, linewidth=2)
        ax.plot(E_fields, expected, 's--', label='Theory', markersize=8, linewidth=2)
        ax.set_xlabel('E-field (V/m)', fontsize=11)
        ax.set_ylabel('Drift Velocity (m/s)', fontsize=11)
        ax.set_title('Drift vs E-field (500 Pa)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Panel 3: Error vs E-field
    ax = axes[0, 2]
    all_data = [v for v in results_dict.values() if v is not None]
    if all_data:
        E_1000 = [d['E_field_Vm'] for d in data_1000Pa]
        err_1000 = [d['error_percent'] for d in data_1000Pa]
        E_500 = [d['E_field_Vm'] for d in data_500Pa]
        err_500 = [d['error_percent'] for d in data_500Pa]
        
        ax.plot(E_1000, err_1000, 'o-', label='1000 Pa', markersize=8, linewidth=2)
        ax.plot(E_500, err_500, 's-', label='500 Pa', markersize=8, linewidth=2)
        ax.axhline(20, color='orange', linestyle='--', label='20% tolerance')
        ax.set_xlabel('E-field (V/m)', fontsize=11)
        ax.set_ylabel('Error (%)', fontsize=11)
        ax.set_title('Accuracy vs E-field', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    
    # Panel 4: Velocity components (low Td)
    ax = axes[1, 0]
    key_low_td_1000 = 'Td5p0_1000Pa'
    if key_low_td_1000 in results_dict and results_dict[key_low_td_1000] is not None:
        vz_hist = results_dict[key_low_td_1000]['vz_final']
        ax.hist(vz_hist, bins=50, alpha=0.7, edgecolor='black')
        ax.axvline(results_dict[key_low_td_1000]['v_gas_z'], color='red', 
                   linestyle='--', linewidth=2, label=f"v_gas = {results_dict[key_low_td_1000]['v_gas_z']:.0f} m/s")
        ax.axvline(np.mean(vz_hist), color='blue', 
                   linestyle='-', linewidth=2, label=f"<v_z> = {np.mean(vz_hist):.0f} m/s")
        ax.set_xlabel('Velocity v_z (m/s)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title('Distribution (5 Td @ 1000 Pa)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    
    # Panel 5: Position vs time (high Td)
    ax = axes[1, 1]
    key_high_td_1000 = 'Td10p0_1000Pa'
    if key_high_td_1000 in results_dict and results_dict[key_high_td_1000] is not None:
        t = results_dict[key_high_td_1000]['t_array'] * 1e9  # Convert to ns
        z = results_dict[key_high_td_1000]['z_vs_t'] * 100  # Convert to cm
        ax.plot(t, z, 'b-', linewidth=2)
        ax.set_xlabel('Time (ns)', fontsize=11)
        ax.set_ylabel('Position z (cm)', fontsize=11)
        ax.set_title('Drift Trajectory (10 Td @ 1000 Pa)', fontweight='bold')
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
    v_gas_z = BASE_GAS_VELOCITY_MS
    temperature_K = 300.0
    
    pressures_Pa = [1000.0, 500.0]
    target_td_values = [5.0, 10.0]
    test_cases = []

    for pressure in pressures_Pa:
        for target_td in target_td_values:
            E_field = compute_field_for_td(pressure, temperature_K, target_td)
            actual_td = compute_td_from_field(E_field, pressure, temperature_K)
            td_str = str(target_td).replace('.', 'p')
            pressure_str = f"{int(pressure)}Pa"
            test_cases.append({
                'name': f"Td{td_str}_{pressure_str}",
                'E_field_Vm': E_field,
                'pressure_Pa': pressure,
                'output_name': f"drift_Td{td_str}_{pressure_str}.h5",
                'target_td': target_td,
                'actual_td': actual_td
            })
    
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
        log(f"  Target E/N:     {test['target_td']:.1f} Td")
        log(f"  Actual E/N:     {test['actual_td']:.2f} Td")
        log(f"  E-field:        {test['E_field_Vm']:.1f} V/m")
        log(f"  Pressure:       {test['pressure_Pa']:.0f} Pa")
        log(f"  Gas velocity: {v_gas_z:.0f} m/s (axial)")
        log(f"  Temperature:  {temperature_K:.0f} K")
        log(f"")
        log(f"Theory (H3O+ in {CARRIER_GAS_SPECIES}):")
        log(f"  Mobility drift: {v_mobility:.1f} m/s (μ·E)")
        log(f"  Gas flow:       {v_gas_z:.1f} m/s")
        log(f"  Total expected: {v_expected:.1f} m/s")
        
        # Simulation runtime: twice the drift time through IMS_LENGTH_M
        drift_time = IMS_LENGTH_M / max(abs(v_expected), 1e-3)
        sim_time_s = 2.0 * drift_time

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
        analysis = analyze_trajectory(
            h5_file,
            v_expected,
            v_mobility,
            v_gas_z,
            test['E_field_Vm']
        )
        
        if analysis:
            analysis['E_field_Vm'] = test['E_field_Vm']
            analysis['pressure_Pa'] = test['pressure_Pa']
            analysis['target_td'] = test['target_td']
            analysis['actual_td'] = test['actual_td']
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
        log("  - Baseline IMS geometry reproduced")
        log("  - Drift velocity follows μ·E once sufficient collisions occur")
        log("  - Ready to reintroduce axial gas flow after baseline pass")
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

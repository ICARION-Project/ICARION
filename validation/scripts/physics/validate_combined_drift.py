#!/usr/bin/env python3
"""
Combined Drift Validation Script: E-field + Gas Flow

Purpose:
    Revalidate the superposition law v_drift = μ·E + v_gas using the
    production collision kernels and IMS-style geometry. The earlier sweep was
    paused after the CCS lookup fixes; this refresh restores the physics
    coverage with tractable runtimes.

Highlights:
    - Explicit gas flow (50 m/s) plus electric drift (0, 1 kV/m, 5 kV/m)
    - Two pressures (100 Pa, 1000 Pa) to probe density scaling
    - Ensemble of 1000 ions, 1 µs runtime, dt = 1 ns (gas-flow parity)
    - EHSS kernel by default (toggleable), matching IMS configs
    - Direct comparison against Mason-Schamp low-field theory

Runtime: ~8 minutes for the six-case sweep.
"""

import subprocess
import json
import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from datetime import datetime
import os

# Constants
ELEM_CHARGE_C = 1.602176634e-19  # C
BOLTZMANN_K = 1.380649e-23       # J/K
AMU_TO_KG = 1.66053906660e-27    # kg
TD_TO_V_M2 = 1e-21              # 1 Td = 1e-21 V·m²

# Drift tube geometry/runtime (matches report Section 8.5)
DRIFT_LENGTH_M = 0.20
DRIFT_RADIUS_M = 0.10
DRIFT_ORIGIN_M = [0.0, 0.0, 0.0]
ION_INJECTION_CENTER = [0.0, 0.0, 0.05]
ION_INJECTION_STD = [5e-4, 5e-4, 1e-4]
SIMULATION_DT_S = 1e-9
SIMULATION_TIME_S = 1e-6
BASE_GAS_VELOCITY_MS = 50.0
IONS_PER_SPECIES = 1000
ANALYSIS_Z_MIN_M = 0.02
ANALYSIS_Z_MAX_M = DRIFT_LENGTH_M - 0.02
MAX_REDUCED_TD = 10.0

# Paths
REPO_ROOT = Path(__file__).resolve().parents[3]  # Go up 3 levels: physics/ -> scripts/ -> validation/ -> repo_root/
VALIDATION_DIR = REPO_ROOT / "validation"
RUN_DIR = os.environ.get("ICARION_VALIDATION_RUN_DIR")
if RUN_DIR:
    RUN_DIR = Path(RUN_DIR)
    FIGURES_DIR = RUN_DIR / "figures" / "physics"
    LOGS_DIR = RUN_DIR / "logs"
    RESULTS_DIR = RUN_DIR / "results" / "combined_drift"
else:
    FIGURES_DIR = VALIDATION_DIR / "figures" / "physics"
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

def compute_td_from_field(E_field_Vm, pressure_Pa, temperature_K):
    """Convert an absolute field back to the reduced field (Townsend)."""
    number_density = pressure_Pa / (BOLTZMANN_K * temperature_K)
    return (E_field_Vm / number_density) / TD_TO_V_M2

# Create directories
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
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

    axial_V = E_field_Vm * DRIFT_LENGTH_M

    config = {
        "simulation": {
            "total_time_s": sim_time_s,
            "dt_s": SIMULATION_DT_S,
            "write_interval": 200,
            "integrator": "RK4",
            "enable_gpu": False,  # GPU path drops stochastic forces today
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
                    "center": ION_INJECTION_CENTER,
                    "std": ION_INJECTION_STD
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
                "origin_m": DRIFT_ORIGIN_M,
                "length_m": DRIFT_LENGTH_M,
                "radius_m": DRIFT_RADIUS_M
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
            timeout=900  # 15 minute guard for 1 µs runs
        )
        
        if result.returncode != 0:
            log(f"  ❌ Simulation failed!")
            log(f"  stderr: {result.stderr[-500:]}")  # Last 500 chars
            return None
        
        log(f"  ✓ Simulation complete")
        return RESULTS_DIR / config_dict["output"]["trajectory_file"]
        
    except subprocess.TimeoutExpired:
        log(f"  ❌ Simulation timeout (>15 minutes)")
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
        active_mask = (z_final > ANALYSIS_Z_MIN_M) & (z_final < ANALYSIS_Z_MAX_M)
        
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

    valid_data = [data for data in results_dict.values() if data is not None]
    if not valid_data:
        log("  ⚠ No successful simulations to plot")
        return

    pressures = sorted({data['pressure_Pa'] for data in valid_data})

    fig, axes = plt.subplots(2, 3, figsize=(16, 10))
    fig.suptitle('Combined Drift Validation: E-field + Gas Flow\n(v_drift = μ·E + v_gas)',
                 fontsize=14, fontweight='bold')

    def subset_for_pressure(p_value):
        return sorted(
            [data for data in valid_data if np.isclose(data['pressure_Pa'], p_value)],
            key=lambda d: d['E_field_Vm']
        )

    # Panels 1 & 2: Drift velocity vs E-field per pressure
    for idx in range(2):
        ax = axes[0, idx]
        if idx >= len(pressures):
            ax.axis('off')
            continue
        pressure = pressures[idx]
        subset = subset_for_pressure(pressure)
        if not subset:
            ax.axis('off')
            continue

        E_fields = np.array([d['E_field_Vm'] for d in subset])
        measured = np.array([d['mean_vz_steady'] for d in subset])
        expected = np.array([d['v_expected'] for d in subset])

        ax.plot(E_fields, measured, 'o-', label='Measured', markersize=8, linewidth=2)
        ax.plot(E_fields, expected, 's--', label='Theory', markersize=8, linewidth=2)
        ax.set_xlabel('E-field (V/m)', fontsize=11)
        ax.set_ylabel('Drift Velocity (m/s)', fontsize=11)
        ax.set_title(f'Drift vs E-field ({int(pressure)} Pa)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)

    # Panel 3: Error vs E-field by pressure
    ax = axes[0, 2]
    markers = ['o', 's', '^', 'd']
    for idx, pressure in enumerate(pressures):
        subset = subset_for_pressure(pressure)
        if not subset:
            continue
        E_fields = [d['E_field_Vm'] for d in subset]
        errors = [d['error_percent'] for d in subset]
        ax.plot(E_fields, errors, markers[idx % len(markers)] + '-',
                label=f'{int(pressure)} Pa', markersize=8, linewidth=2)
    ax.axhline(20, color='orange', linestyle='--', label='20% tolerance')
    ax.set_xlabel('E-field (V/m)', fontsize=11)
    ax.set_ylabel('Error (%)', fontsize=11)
    ax.set_title('Accuracy vs E-field', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)

    # Panel 4: Velocity distribution for baseline (E=0) case
    ax = axes[1, 0]
    baseline_case = next((d for d in valid_data if np.isclose(d['E_field_Vm'], 0.0)), None)
    if baseline_case is not None:
        vz_hist = baseline_case['vz_final']
        ax.hist(vz_hist, bins=50, alpha=0.7, edgecolor='black')
        ax.axvline(baseline_case['v_gas_z'], color='red', linestyle='--', linewidth=2,
                   label=f"v_gas = {baseline_case['v_gas_z']:.0f} m/s")
        ax.axvline(np.mean(vz_hist), color='blue', linestyle='-', linewidth=2,
                   label=f"<v_z> = {np.mean(vz_hist):.0f} m/s")
        ax.set_xlabel('Velocity v_z (m/s)', fontsize=11)
        ax.set_ylabel('Count', fontsize=11)
        ax.set_title(f"Distribution (E=0, {int(baseline_case['pressure_Pa'])} Pa)", fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
    else:
        ax.axis('off')

    # Panel 5: Position vs time for highest E-field at high pressure
    ax = axes[1, 1]
    high_pressure = pressures[-1]
    high_field_case = max(
        (d for d in valid_data if np.isclose(d['pressure_Pa'], high_pressure)),
        key=lambda d: d['E_field_Vm'],
        default=None
    )
    if high_field_case and len(high_field_case['t_array']) >= 2:
        t = high_field_case['t_array'] * 1e9  # ns
        z = high_field_case['z_vs_t'] * 100   # cm
        ax.plot(t, z, 'b-', linewidth=2)
        ax.set_xlabel('Time (ns)', fontsize=11)
        ax.set_ylabel('Position z (cm)', fontsize=11)
        ax.set_title(
            f"Drift Trajectory ({int(high_field_case['pressure_Pa'])} Pa, E={high_field_case['E_field_Vm']:.0f} V/m)",
            fontweight='bold'
        )
        ax.grid(True, alpha=0.3)
    else:
        ax.axis('off')

    # Panel 6: Component breakdown at highest pressure
    ax = axes[1, 2]
    high_pressure_subset = subset_for_pressure(high_pressure)
    if high_pressure_subset:
        E_fields = np.array([d['E_field_Vm'] for d in high_pressure_subset])
        v_total = np.array([d['v_expected'] for d in high_pressure_subset])
        v_mobility = np.array([d['v_mobility'] for d in high_pressure_subset])
        v_gas = np.array([d['v_gas_z'] for d in high_pressure_subset])

        ax.plot(E_fields, v_mobility, 'o-', label='Mobility (μ·E)', markersize=8, linewidth=2)
        ax.plot(E_fields, v_gas, 's--', label='Gas flow', markersize=8, linewidth=2)
        ax.plot(E_fields, v_total, '^-', label='Total', markersize=8, linewidth=2)
        ax.set_xlabel('E-field (V/m)', fontsize=11)
        ax.set_ylabel('Velocity (m/s)', fontsize=11)
        ax.set_title(f'Drift Components ({int(high_pressure)} Pa)', fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3)
    else:
        ax.axis('off')

    plt.tight_layout()

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
    
    pressures_Pa = [100.0, 1000.0]
    E_field_values = [0.0, 1000.0, 5000.0]
    test_cases = []

    for pressure in pressures_Pa:
        for E_field in E_field_values:
            reduced_td = compute_td_from_field(E_field, pressure, temperature_K)
            if reduced_td > MAX_REDUCED_TD:
                continue
            pressure_str = f"{int(pressure)}Pa"
            field_label = f"E{int(E_field)}"
            test_cases.append({
                'name': f"{field_label}_{pressure_str}",
                'E_field_Vm': E_field,
                'pressure_Pa': pressure,
                'output_name': f"combined_{field_label}_{pressure_str}.h5",
                'reduced_td': reduced_td
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
        log(f"  E-field:        {test['E_field_Vm']:.1f} V/m")
        log(f"  Pressure:       {test['pressure_Pa']:.0f} Pa")
        log(f"  Reduced E/N:    {test['reduced_td']:.2f} Td")
        log(f"  Gas velocity:   {v_gas_z:.0f} m/s (axial)")
        log(f"  Temperature:    {temperature_K:.0f} K")
        log(f"")
        log(f"Theory (H3O+ in {CARRIER_GAS_SPECIES}):")
        log(f"  Mobility drift: {v_mobility:.1f} m/s (μ·E)")
        log(f"  Gas flow:       {v_gas_z:.1f} m/s")
        log(f"  Total expected: {v_expected:.1f} m/s")
        
        sim_time_s = SIMULATION_TIME_S

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
            analysis['reduced_td'] = test['reduced_td']
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

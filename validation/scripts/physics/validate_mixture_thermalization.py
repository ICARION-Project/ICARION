#!/usr/bin/env python3
"""
ICARION Physics Validation: Gas Mixture Thermalization

Simple test: Verify ions thermalize correctly in gas mixtures (no E-field).
Start ions at 0.1 K, let them thermalize via collisions to ambient T.

Test Species: H3O+ in He/N2 mixtures
Conditions: P=1000 Pa, T=300 K, E=0 V/m
Duration: ~50 collision times
"""

import subprocess
import json
import h5py
import numpy as np
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
import sys

# Shared HDF5 helpers (species IDs)
COMMON_DIR = Path(__file__).resolve().parents[2] / "common"
if str(COMMON_DIR) not in sys.path:
    sys.path.append(str(COMMON_DIR))
from hdf5_utils import load_species_ids  # noqa: E402

# ============================================================================
# CONFIGURATION
# ============================================================================

# Test conditions
P_Pa = 1000.0
T_K = 300.0
N_ions = 500

# Ion species
ion_species = "H3O+"
ion_mass_amu = 19.02

# Mole fractions to test: He/N2 ratios
test_mixtures = [
    {"name": "100He_0N2",   "He": 1.00, "N2": 0.00},
    {"name": "75He_25N2",   "He": 0.75, "N2": 0.25},
    {"name": "50He_50N2",   "He": 0.50, "N2": 0.50},
    {"name": "25He_75N2",   "He": 0.25, "N2": 0.75},
    {"name": "0He_100N2",   "He": 0.00, "N2": 1.00},
]

# Physical constants
k_B = 1.380649e-23  # J/K
amu_to_kg = 1.66053906660e-27

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
ICARION_BIN = PROJECT_ROOT / "build" / "src" / "icarion_main"
RESULTS_DIR = PROJECT_ROOT / "validation" / "results" / "mixture_thermalization"
FIGURES_DIR = PROJECT_ROOT / "validation" / "figures"
LOGS_DIR = PROJECT_ROOT / "validation" / "logs"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOGS_DIR / "MIXTURE_THERMALIZATION_VALIDATION.txt"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def log(msg, file=None):
    """Print and log message"""
    print(msg)
    if file:
        file.write(msg + "\n")
        file.flush()

def calculate_collision_frequency(m_ion_amu, x_He, x_N2, T_K, P_Pa):
    """Calculate weighted collision frequency for gas mixture"""
    # CCS values for H3O+ (Å²)
    CCS_He = 25.56
    CCS_N2 = 104.02
    
    # Gas masses
    m_He_amu = 4.002602
    m_N2_amu = 28.014
    
    # Convert to kg
    m_ion = m_ion_amu * amu_to_kg
    m_He = m_He_amu * amu_to_kg
    m_N2 = m_N2_amu * amu_to_kg
    
    # Number density (total)
    N_total = P_Pa / (k_B * T_K)
    
    # Collision frequencies for each gas
    nu_total = 0.0
    
    if x_He > 0:
        mu_He = (m_ion * m_He) / (m_ion + m_He)
        v_rel_He = np.sqrt(8 * k_B * T_K / (np.pi * mu_He))
        sigma_He = CCS_He * 1e-20  # m²
        nu_He = x_He * N_total * sigma_He * v_rel_He
        nu_total += nu_He
    
    if x_N2 > 0:
        mu_N2 = (m_ion * m_N2) / (m_ion + m_N2)
        v_rel_N2 = np.sqrt(8 * k_B * T_K / (np.pi * mu_N2))
        sigma_N2 = CCS_N2 * 1e-20  # m²
        nu_N2 = x_N2 * N_total * sigma_N2 * v_rel_N2
        nu_total += nu_N2
    
    tau_collision = 1.0 / nu_total if nu_total > 0 else 1e-6
    
    return nu_total, tau_collision

def create_config(mixture_name, x_He, x_N2, output_file):
    """Create simulation config for mixture thermalization test"""
    
    # Calculate collision time for this mixture
    nu, tau_coll = calculate_collision_frequency(ion_mass_amu, x_He, x_N2, T_K, P_Pa)
    
    # Simulation parameters
    n_collision_times = 2000  # Very long time to check convergence
    total_time_s = n_collision_times * tau_coll
    dt_s = tau_coll / 50.0  # 50 steps per collision time
    write_interval = max(1, int(n_collision_times // 10))  # ~10 writes total
    
    # Domain size
    length_m = 0.2
    radius_m = 0.1
    
    # Build gas_mixture array
    gas_mixture = []
    if x_He > 0:
        gas_mixture.append({"species": "He", "mole_fraction": x_He})
    if x_N2 > 0:
        gas_mixture.append({"species": "N2", "mole_fraction": x_N2})
    
    config = {
        "simulation": {
            "total_time_s": total_time_s,
            "dt_s": dt_s,
            "write_interval": write_interval,
            "integrator": "RK4",
            "enable_gpu": False,  # Multi-gas only on CPU
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "HSS",
            "enable_space_charge": False
        },
        "output": {
            "folder": str(RESULTS_DIR),
            "trajectory_file": f"therm_{mixture_name}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "id": ion_species,
                    "count": N_ions,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.1],
                        "std": [0.001, 0.001, 0.001]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": 0.1  # Cold start
                    }
                }
            ]
        },
        "species_database": str(PROJECT_ROOT / "data" / "species_database_v1.json"),
        "domains": [
            {
                "name": f"therm_test_{mixture_name}",
                "instrument": "IMS",
                "geometry": {
                    "origin_m": [0.0, 0.0, 0.0],
                    "length_m": length_m,
                    "radius_m": radius_m
                },
                "env": {
                    "temperature_K": T_K,
                    "pressure_Pa": P_Pa,
                    "gas_mixture": gas_mixture
                },
                "fields": {
                    "DC": {
                        "axial_V": 0.0  # No E-field
                    }
                }
            }
        ]
    }
    
    with open(output_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    return output_file, tau_coll, total_time_s

def analyze_thermalization(h5_file, mixture_name, T_expected, logf):
    """Analyze thermalization from simulation output"""
    
    log(f"\n📊 Analyzing: {h5_file.name}", logf)
    
    with h5py.File(h5_file, 'r') as f:
        velocities = f['trajectory/velocities'][:]  # (n_frames, n_ions, 3)
        times = f['trajectory/time'][:]  # (n_frames,)
        species_ids = load_species_ids(f)  # (n_frames, n_ions) or (n_ions,)
        if species_ids.ndim == 1:
            species_ids = species_ids[np.newaxis, :]
        
        n_frames, n_ions, _ = velocities.shape
        
        # Calculate temperature at each frame
        T_ions = np.zeros(n_frames)
        n_active = np.zeros(n_frames)
        
        for i in range(n_frames):
            active_mask = np.array([s != '' for s in species_ids[i, :]])
            n_active[i] = np.sum(active_mask)
            
            if n_active[i] > 0:
                # KE = 0.5 * m * v²  =>  T = (m * <v²>) / (3 * k_B)
                v_squared = np.sum(velocities[i, active_mask, :]**2, axis=1)
                mean_v_squared = np.mean(v_squared)
                T_ions[i] = (ion_mass_amu * amu_to_kg * mean_v_squared) / (3.0 * k_B)
        
        # Final temperature (average of last 10% of trajectory)
        final_idx = int(0.9 * n_frames)
        T_final = np.mean(T_ions[final_idx:])
        T_final_std = np.std(T_ions[final_idx:])
        
        error_pct = abs(T_final - T_expected) / T_expected * 100
        status = "✅" if error_pct < 10 else "⚠️" if error_pct < 20 else "❌"
        
        log(f"  Ions active (final): {int(n_active[-1])}/{n_ions}", logf)
        log(f"  Temperature (final): {T_final:.1f} ± {T_final_std:.1f} K", logf)
        log(f"  Expected:            {T_expected:.1f} K", logf)
        log(f"  Error:               {error_pct:.1f}%  {status}", logf)
        
        return {
            'T_final': T_final,
            'T_final_std': T_final_std,
            'T_expected': T_expected,
            'error_pct': error_pct,
            'status': status,
            'times': times,
            'T_ions': T_ions,
            'n_active': n_active,
            'mixture_name': mixture_name
        }

def run_simulation(config_file, logf):
    """Run ICARION simulation"""
    cmd = [str(ICARION_BIN), str(config_file)]
    
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=PROJECT_ROOT)
    
    if result.returncode != 0:
        log(f"  ❌ Simulation failed!", logf)
        log(f"  STDERR: {result.stderr}", logf)
        return False
    
    log(f"  ✓ Simulation complete", logf)
    return True

def create_validation_plot(results, output_file):
    """Create validation plot"""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Plot 1: Temperature evolution
    ax = axes[0, 0]
    for r in results:
        ax.plot(r['times'] * 1e3, r['T_ions'], label=r['mixture_name'], linewidth=2)
    ax.axhline(y=T_K, color='black', linestyle='--', linewidth=2, label='Expected (300 K)')
    ax.set_xlabel("Time [ms]", fontsize=12)
    ax.set_ylabel("Temperature [K]", fontsize=12)
    ax.set_title("Thermalization in Gas Mixtures", fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Error bar chart
    ax = axes[0, 1]
    names = [r['mixture_name'] for r in results]
    errors = [r['error_pct'] for r in results]
    colors = ['green' if abs(e) < 10 else 'orange' if abs(e) < 20 else 'red' for e in errors]
    ax.bar(range(len(errors)), errors, color=colors, alpha=0.7, edgecolor='black')
    ax.axhline(y=10, color='orange', linestyle='--', label='10% threshold')
    ax.axhline(y=-10, color='orange', linestyle='--')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xticks(range(len(errors)))
    ax.set_xticklabels([n.replace('_', '\n') for n in names], fontsize=10)
    ax.set_ylabel("Error [%]", fontsize=12)
    ax.set_title("Temperature Error", fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 3: Final temperature vs mixture
    ax = axes[1, 0]
    T_final = [r['T_final'] for r in results]
    T_std = [r['T_final_std'] for r in results]
    x_He = [float(n.split('He')[0]) / 100 for n in names]
    ax.errorbar(x_He, T_final, yerr=T_std, fmt='o', markersize=10, linewidth=2, capsize=5)
    ax.axhline(y=T_K, color='red', linestyle='--', linewidth=2, label='Expected (300 K)')
    ax.set_xlabel("He Mole Fraction", fontsize=12)
    ax.set_ylabel("Final Temperature [K]", fontsize=12)
    ax.set_title("Final Temperature vs Gas Composition", fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Plot 4: Number of active ions
    ax = axes[1, 1]
    for r in results:
        ax.plot(r['times'] * 1e3, r['n_active'], label=r['mixture_name'], linewidth=2)
    ax.set_xlabel("Time [ms]", fontsize=12)
    ax.set_ylabel("Active Ions", fontsize=12)
    ax.set_title("Ion Survival", fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"  ✓ Saved: {output_file}")

# ============================================================================
# MAIN VALIDATION
# ============================================================================

def main():
    start_time = datetime.now()
    
    with open(LOG_FILE, 'w') as logf:
        log("="*80, logf)
        log("GAS MIXTURE THERMALIZATION VALIDATION", logf)
        log("="*80, logf)
        log(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}", logf)
        log("", logf)
        log("Purpose: Validate HSS collision model thermalization in gas mixtures", logf)
        log("Method: Cold start (0.1 K) → thermalize to ambient (300 K)", logf)
        log("", logf)
        log(f"Test Species: {ion_species}", logf)
        log(f"Gas Mixture: He/N2 (5 ratios)", logf)
        log(f"Conditions: P={P_Pa} Pa, T={T_K} K, E=0 V/m", logf)
        log(f"Duration: ~50 collision times each", logf)
        log("", logf)
        
        results = []
        
        for mixture in test_mixtures:
            name = mixture['name']
            x_He = mixture['He']
            x_N2 = mixture['N2']
            
            log("─"*80, logf)
            log(f"TEST: {name}", logf)
            log("─"*80, logf)
            log(f"Composition: He={x_He*100:.0f}%, N2={x_N2*100:.0f}%", logf)
            
            # Calculate collision parameters
            nu, tau_coll = calculate_collision_frequency(ion_mass_amu, x_He, x_N2, T_K, P_Pa)
            log(f"Collision frequency: {nu:.2e} Hz", logf)
            log(f"Collision time:      {tau_coll*1e6:.2f} µs", logf)
            log("", logf)
            
            # Create config
            config_file, tau, total_time = create_config(name, x_He, x_N2, RESULTS_DIR / f"config_{name}.json")
            log(f"▶ Running simulation: {config_file.name}", logf)
            log(f"  Total time: {total_time*1e3:.3f} ms ({total_time/tau:.0f} collision times)", logf)
            
            # Run simulation
            if not run_simulation(config_file, logf):
                log(f"  ❌ FAILED - skipping analysis", logf)
                continue
            
            # Analyze results
            h5_file = RESULTS_DIR / f"therm_{name}.h5"
            analysis = analyze_thermalization(h5_file, name, T_K, logf)
            
            if analysis is None:
                log(f"  ❌ FAILED - analysis error", logf)
                continue
            
            results.append(analysis)
        
        # Summary
        log("", logf)
        log("="*80, logf)
        log("VALIDATION SUMMARY", logf)
        log("="*80, logf)
        
        for r in results:
            log(f"  {r['mixture_name']:15s}: {r['status']} {r['error_pct']:6.1f}% error", logf)
        
        log("", logf)
        
        n_pass = sum(1 for r in results if r['error_pct'] < 10)
        n_total = len(results)
        
        if n_pass == n_total:
            log(f"✅ ALL VALIDATIONS PASSED ({n_pass}/{n_total})", logf)
        elif n_pass >= 0.8 * n_total:
            log(f"⚠️ MOST VALIDATIONS PASSED ({n_pass}/{n_total})", logf)
        else:
            log(f"❌ SOME VALIDATIONS FAILED ({n_pass}/{n_total})", logf)
        
        # Create plots
        if results:
            log("", logf)
            log("📈 Creating validation plot...", logf)
            plot_file = FIGURES_DIR / "mixture_thermalization_validation.png"
            create_validation_plot(results, plot_file)
            log(f"  ✓ Saved: {plot_file}", logf)
        
        end_time = datetime.now()
        log("", logf)
        log("="*80, logf)
        log(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}", logf)
        log(f"Duration: {(end_time - start_time).total_seconds():.1f} seconds", logf)
        log(f"Log file: {LOG_FILE}", logf)
        log("="*80, logf)
    
    # Exit code
    n_pass = sum(1 for r in results if r['error_pct'] < 10)
    return 0 if n_pass == len(results) else 1

if __name__ == "__main__":
    exit(main())

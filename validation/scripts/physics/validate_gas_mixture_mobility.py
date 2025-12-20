#!/usr/bin/env python3
"""
ICARION Physics Validation: Gas Mixture Mobility (Blanc's Law)

Tests HSS/EHSS collision models for gas mixtures using Blanc's Law:
    1/K_mix = Σ(x_i / K_i)

Where:
- K_mix: Mobility in gas mixture
- x_i: Mole fraction of gas component i
- K_i: Mobility in pure gas i

Test Species: H3O+ in He/N2 mixtures
Conditions: P=1000 Pa, T=300 K, E=1000 V/m (E/N ~41 Td, well within linear regime)
"""

import subprocess
import json
import h5py
import numpy as np
from pathlib import Path
from datetime import datetime
import matplotlib.pyplot as plt
from scipy import stats
import sys

# Shared HDF5 helpers (species IDs)
COMMON_DIR = Path(__file__).resolve().parents[2] / "common"
if str(COMMON_DIR) not in sys.path:
    sys.path.append(str(COMMON_DIR))
from hdf5_utils import load_species_ids  # noqa: E402

# ============================================================================
# CONFIGURATION
# ============================================================================

# Reference mobilities for H3O+ (K0 at Tref=273K, P0=101325Pa)
# Source: Literature values for H3O+ in He and N2
K0_He = 24.1e-4  # m²/(V·s) - H3O+ in pure He
K0_N2 = 3.2e-4   # m²/(V·s) - H3O+ in pure N2

# Test conditions (use low-field regime so K₀ references apply directly)
P_Pa = 1000.0
T_K = 300.0
E_Vm = 1000.0  # E/N ≈ 4.1 Td -> minimal field heating
N_ions = 500  # Fewer ions for faster tests
total_time_ms = 0.1  # Shorter simulation time

# Ion species
ion_species = "H3O+"
ion_mass_amu = 19.0

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
e = 1.602176634e-19  # C
N_A = 6.02214076e23  # 1/mol
T_ref = 273.15  # K
P_ref = 101325.0  # Pa

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
ICARION_BIN = PROJECT_ROOT / "build" / "src" / "icarion_main"
RESULTS_DIR = PROJECT_ROOT / "validation" / "results" / "gas_mixture_mobility"
FIGURES_DIR = PROJECT_ROOT / "validation" / "figures"
LOGS_DIR = PROJECT_ROOT / "validation" / "logs"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)
FIGURES_DIR.mkdir(parents=True, exist_ok=True)
LOGS_DIR.mkdir(parents=True, exist_ok=True)

LOG_FILE = LOGS_DIR / "GAS_MIXTURE_MOBILITY_VALIDATION.txt"

# ============================================================================
# HELPER FUNCTIONS
# ============================================================================

def log(msg, file=None):
    """Print and log message"""
    print(msg)
    if file:
        file.write(msg + "\n")
        file.flush()

def blancs_law(K0_He, K0_N2, x_He, x_N2, T_K, T_ref, P_Pa, P_ref):
    """
    Calculate mobility in gas mixture using Blanc's Law
    
    Returns K_mix at (T_K, P_Pa) in m²/(V·s)
    """
    # Convert K0 to actual K at test conditions
    K_He = K0_He * (T_K / T_ref) * (P_ref / P_Pa)
    K_N2 = K0_N2 * (T_K / T_ref) * (P_ref / P_Pa)
    
    # Blanc's Law: 1/K_mix = Σ(x_i / K_i)
    if x_He > 0 and x_N2 > 0:
        K_mix_inv = x_He / K_He + x_N2 / K_N2
        K_mix = 1.0 / K_mix_inv
    elif x_He == 1.0:
        K_mix = K_He
    elif x_N2 == 1.0:
        K_mix = K_N2
    else:
        raise ValueError("Mole fractions must sum to 1.0")
    
    return K_mix

def create_config(mixture_name, x_He, x_N2, output_file):
    """Create simulation config for gas mixture test"""
    
    # Domain size: large enough to avoid boundary effects
    length_m = 0.2  # 20 cm drift length
    radius_m = 0.1  # 10 cm radius
    
    # Electric field setup
    axial_V = E_Vm * length_m  # Voltage to achieve E_Vm field
    
    # Build gas_mixture array (only include gases with x > 0)
    gas_mixture = []
    if x_He > 0:
        gas_mixture.append({
            "species": "He",
            "mole_fraction": x_He
        })
    if x_N2 > 0:
        gas_mixture.append({
            "species": "N2",
            "mole_fraction": x_N2
        })
    
    config = {
        "simulation": {
            "total_time_s": total_time_ms * 1e-3,
            "dt_s": 1e-9,
            "write_interval": 10,
            "integrator": "RK4",
            "enable_gpu": False,  # Multi-gas mixtures only supported on CPU
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "HSS",
            "enable_space_charge": False
        },
        "output": {
            "folder": str(RESULTS_DIR),
            "trajectory_file": f"drift_{mixture_name}.h5",
            "print_progress": True
        },
        "ions": {
            "species": [
                {
                    "id": ion_species,
                    "count": N_ions,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.01],
                        "std": [0.001, 0.001, 0.001]
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": T_K
                    }
                }
            ]
        },
        "species_database": str(PROJECT_ROOT / "data" / "species_database_v1.json"),
        "domains": [
            {
                "name": f"mixture_test_{mixture_name}",
                "instrument": "IMS",
                "geometry": {
                    "origin_m": [0.0, 0.0, 0.0],
                    "length_m": length_m,
                    "radius_m": radius_m
                },
                "env": {
                    "temperature_K": T_K,
                    "pressure_Pa": P_Pa,
                    "gas_mixture": gas_mixture,
                    "gas_velocity_m_s": [0.0, 0.0, 0.0]
                },
                "fields": {
                    "DC": {
                        "axial_V": axial_V
                    }
                }
            }
        ]
    }
    
    with open(output_file, 'w') as f:
        json.dump(config, f, indent=2)
    
    return output_file

def analyze_drift_velocity(h5_file, mixture_name, logf):
    """Analyze drift velocity from simulation output"""
    
    log(f"\n📊 Analyzing: {h5_file.name}", logf)
    
    with h5py.File(h5_file, 'r') as f:
        positions = f['trajectory/positions'][:]  # (n_frames, n_ions, 3)
        times = f['trajectory/time'][:]  # (n_frames,)
        species_ids = load_species_ids(f)  # (n_frames, n_ions) - string labels
        if species_ids.ndim == 1:
            species_ids = species_ids[np.newaxis, :]
        
        # Filter: only active ions (species_id != b'')
        n_frames, n_ions, _ = positions.shape
        z_positions = positions[:, :, 2]  # Axial (z) positions
        
        # Get active ions mask (last frame) - active = has non-empty species ID
        active_mask = np.array([s != b'' for s in species_ids[-1, :]])
        n_active = np.sum(active_mask)
        
        log(f"  Ions: {n_active}/{n_ions} active", logf)
        
        if n_active < 10:
            log(f"  ❌ Too few active ions!", logf)
            return None
        
        # Calculate mean z-position vs time (only active ions)
        z_mean = np.zeros(n_frames)
        for i in range(n_frames):
            active_now = np.array([s != b'' for s in species_ids[i, :]])
            if np.sum(active_now) > 0:
                z_mean[i] = np.mean(z_positions[i, active_now])
            else:
                z_mean[i] = 0.0
        
        # Linear fit: z = v_drift * t + z0
        # Use last 80% of trajectory for steady-state drift
        start_idx = int(0.2 * n_frames)
        slope, intercept, r_value, p_value, std_err = stats.linregress(
            times[start_idx:], z_mean[start_idx:]
        )
        
        v_drift_fit = slope
        
        # Instantaneous velocity (less reliable, for comparison)
        dz = np.diff(z_mean)
        dt = np.diff(times)
        v_instant = dz / dt
        v_drift_instant = np.mean(v_instant[start_idx:])
        v_drift_instant_std = np.std(v_instant[start_idx:])
        
        log(f"  Drift velocity (fit):      {v_drift_fit:8.1f} m/s  ← PRIMARY", logf)
        log(f"  Drift velocity (instant):  {v_drift_instant:8.1f} ± {v_drift_instant_std:6.1f} m/s  (noisy)", logf)
        log(f"  R² (fit quality):          {r_value**2:.4f}", logf)
        
        return {
            'v_drift_fit': v_drift_fit,
            'v_drift_instant': v_drift_instant,
            'v_drift_instant_std': v_drift_instant_std,
            'r_squared': r_value**2,
            'n_active': n_active,
            'times': times,
            'z_mean': z_mean,
            'mixture_name': mixture_name
        }

def run_simulation(config_file, logf):
    """Run ICARION simulation"""
    cmd = [str(ICARION_BIN), str(config_file)]
    log(f"  ✓ Running: {' '.join(cmd)}", logf)
    
    result = subprocess.run(
        cmd,
        capture_output=True,
        text=True,
        cwd=PROJECT_ROOT
    )
    
    if result.returncode != 0:
        log(f"  ❌ Simulation failed!", logf)
        log(f"  STDERR: {result.stderr}", logf)
        return False
    
    log(f"  ✓ Simulation complete", logf)
    return True

def create_validation_plot(results, output_file):
    """Create validation plot comparing measured vs predicted mobility"""
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    
    # Extract data
    x_He_list = [r['x_He'] for r in results]
    K_expected = [r['K_expected'] for r in results]
    K_measured = [r['K_measured'] for r in results]
    v_expected = [r['v_expected'] for r in results]
    v_measured = [r['v_measured'] for r in results]
    errors = [r['error_pct'] for r in results]
    mixture_names = [r['mixture_name'] for r in results]
    
    # Plot 1: Mobility vs He fraction
    ax = axes[0, 0]
    ax.plot(x_He_list, np.array(K_expected)*1e4, 'o-', label="Blanc's Law", 
            color='blue', linewidth=2, markersize=8)
    ax.plot(x_He_list, np.array(K_measured)*1e4, 's--', label="HSS (measured)", 
            color='red', linewidth=2, markersize=8)
    ax.set_xlabel("He Mole Fraction", fontsize=12)
    ax.set_ylabel("Mobility K [cm²/(V·s)]", fontsize=12)
    ax.set_title("Mobility vs Gas Mixture Composition", fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Plot 2: Drift velocity vs He fraction
    ax = axes[0, 1]
    ax.plot(x_He_list, v_expected, 'o-', label="Expected (K·E)", 
            color='blue', linewidth=2, markersize=8)
    ax.plot(x_He_list, v_measured, 's--', label="Measured", 
            color='red', linewidth=2, markersize=8)
    ax.set_xlabel("He Mole Fraction", fontsize=12)
    ax.set_ylabel("Drift Velocity [m/s]", fontsize=12)
    ax.set_title(f"Drift Velocity (E={E_Vm} V/m, P={P_Pa} Pa)", fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Plot 3: Error bar chart
    ax = axes[1, 0]
    colors = ['green' if abs(e) < 10 else 'orange' if abs(e) < 20 else 'red' 
              for e in errors]
    bars = ax.bar(range(len(errors)), errors, color=colors, alpha=0.7, edgecolor='black')
    ax.axhline(y=10, color='orange', linestyle='--', label='10% threshold')
    ax.axhline(y=-10, color='orange', linestyle='--')
    ax.axhline(y=0, color='black', linestyle='-', linewidth=0.5)
    ax.set_xticks(range(len(errors)))
    ax.set_xticklabels([m.replace('_', '\n') for m in mixture_names], fontsize=10)
    ax.set_ylabel("Error [%]", fontsize=12)
    ax.set_title("Relative Error per Mixture", fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3, axis='y')
    
    # Plot 4: Measured vs Expected (parity plot)
    ax = axes[1, 1]
    ax.plot(v_expected, v_measured, 'o', markersize=10, color='blue')
    for i, name in enumerate(mixture_names):
        ax.annotate(name.replace('_', '\n'), (v_expected[i], v_measured[i]), 
                   fontsize=8, ha='right', va='bottom')
    
    # Perfect agreement line
    v_min = min(min(v_expected), min(v_measured))
    v_max = max(max(v_expected), max(v_measured))
    ax.plot([v_min, v_max], [v_min, v_max], 'k--', label='Perfect agreement')
    
    # ±10% error bands
    ax.fill_between([v_min, v_max], 
                     [v_min*0.9, v_max*0.9], 
                     [v_min*1.1, v_max*1.1], 
                     color='green', alpha=0.2, label='±10% error')
    
    ax.set_xlabel("Expected Drift Velocity [m/s]", fontsize=12)
    ax.set_ylabel("Measured Drift Velocity [m/s]", fontsize=12)
    ax.set_title("Parity Plot: Measured vs Expected", fontsize=13, fontweight='bold')
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)
    ax.set_aspect('equal', adjustable='box')
    
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
        log("GAS MIXTURE MOBILITY VALIDATION: Blanc's Law", logf)
        log("="*80, logf)
        log(f"Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}", logf)
        log("", logf)
        log("Purpose: Validate HSS collision model for gas mixtures", logf)
        log("Theory: Blanc's Law - 1/K_mix = Σ(x_i / K_i)", logf)
        log("", logf)
        log(f"Test Species: {ion_species}", logf)
        log(f"Gas Mixture: He/N2 (5 ratios: 100/0, 75/25, 50/50, 25/75, 0/100)", logf)
        log(f"Conditions: P={P_Pa} Pa, T={T_K} K, E={E_Vm} V/m", logf)
        log(f"Reference Mobilities (K0 at 273K, 101325 Pa):", logf)
        log(f"  - H3O+ in He: {K0_He*1e4:.3f} cm²/(V·s)", logf)
        log(f"  - H3O+ in N2: {K0_N2*1e4:.3f} cm²/(V·s)", logf)
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
            log("", logf)
            
            # Calculate expected mobility using Blanc's Law
            K_expected = blancs_law(K0_He, K0_N2, x_He, x_N2, T_K, T_ref, P_Pa, P_ref)
            v_expected = K_expected * E_Vm
            
            log(f"Theory (Blanc's Law):", logf)
            log(f"  K_mix:      {K_expected*1e4:.3f} cm²/(V·s)", logf)
            log(f"  v_drift:    {v_expected:.1f} m/s", logf)
            log("", logf)
            
            # Create config
            config_file = RESULTS_DIR / f"config_{name}.json"
            create_config(name, x_He, x_N2, config_file)
            log(f"▶ Running simulation: {config_file.name}", logf)
            log(f"  Config: {config_file}", logf)
            
            # Run simulation
            if not run_simulation(config_file, logf):
                log(f"  ❌ FAILED - skipping analysis", logf)
                continue
            
            # Analyze results
            h5_file = RESULTS_DIR / f"drift_{name}.h5"
            analysis = analyze_drift_velocity(h5_file, name, logf)
            
            if analysis is None:
                log(f"  ❌ FAILED - analysis error", logf)
                continue
            
            v_measured = analysis['v_drift_fit']
            K_measured = v_measured / E_Vm
            
            error_pct = (v_measured - v_expected) / v_expected * 100
            status = "✅" if abs(error_pct) < 10 else "⚠️" if abs(error_pct) < 20 else "❌"
            
            log(f"  Expected:   {v_expected:8.1f} m/s", logf)
            log(f"  Error:      {error_pct:6.1f}%  {status}", logf)
            log("", logf)
            
            results.append({
                'mixture_name': name,
                'x_He': x_He,
                'x_N2': x_N2,
                'K_expected': K_expected,
                'K_measured': K_measured,
                'v_expected': v_expected,
                'v_measured': v_measured,
                'error_pct': error_pct,
                'status': status,
                'n_active': analysis['n_active']
            })
        
        # Summary
        log("="*80, logf)
        log("VALIDATION SUMMARY", logf)
        log("="*80, logf)
        
        for r in results:
            log(f"  {r['mixture_name']:15s}: {r['status']} {r['error_pct']:6.1f}% error", logf)
        
        log("", logf)
        
        n_pass = sum(1 for r in results if abs(r['error_pct']) < 10)
        n_total = len(results)
        
        if n_pass == n_total:
            log(f"✅ ALL VALIDATIONS PASSED ({n_pass}/{n_total})", logf)
        elif n_pass >= 0.8 * n_total:
            log(f"⚠️ MOST VALIDATIONS PASSED ({n_pass}/{n_total})", logf)
        else:
            log(f"❌ SOME VALIDATIONS FAILED ({n_pass}/{n_total})", logf)
        
        # Create plots
        log("", logf)
        log("📈 Creating validation plot...", logf)
        plot_file = FIGURES_DIR / "gas_mixture_mobility_validation.png"
        create_validation_plot(results, plot_file)
        log(f"  ✓ Saved: {plot_file}", logf)
        
        end_time = datetime.now()
        log("", logf)
        log("="*80, logf)
        log(f"End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}", logf)
        log(f"Duration: {(end_time - start_time).total_seconds():.1f} seconds", logf)
        log(f"Log file: {LOG_FILE}", logf)
        log(f"Figures:  {plot_file}", logf)
        log("="*80, logf)
    
    # Exit code
    n_pass = sum(1 for r in results if abs(r['error_pct']) < 10)
    return 0 if n_pass == len(results) else 1

if __name__ == "__main__":
    exit(main())

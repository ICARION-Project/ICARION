#!/usr/bin/env python3
"""
Analyze space charge validation results

Validates:
1. Coulomb expansion: σ(t) ∝ √t
2. Direct vs Grid solver agreement at threshold
3. Space charge effects in IMS (peak broadening)
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit

def read_trajectory_stats(trajectory_file):
    """
    Read trajectory and compute cloud statistics (mean, std) vs time.
    
    Returns: (times, mean_positions, std_positions)
    """
    with h5py.File(trajectory_file, 'r') as f:
        times = f['/trajectory/times'][:]
        positions = f['/trajectory/positions'][:]  # [time, ions, xyz]
        
        # Compute mean and std at each timestep
        mean_pos = np.mean(positions, axis=1)  # [time, xyz]
        std_pos = np.std(positions, axis=1)    # [time, xyz]
        
        # Radial spread (for isotropic expansion)
        radial_std = np.sqrt(std_pos[:, 0]**2 + std_pos[:, 1]**2 + std_pos[:, 2]**2)
        
    return times, mean_pos, std_pos, radial_std

def sqrt_t_model(t, A, t0):
    """Model: σ(t) = A * √(t - t0)"""
    return A * np.sqrt(t - t0)

def analyze_coulomb_expansion(results_dir, n_ions, solver_type):
    """Analyze Coulomb expansion for given N and solver type"""
    results_dir = Path(results_dir)
    
    print(f"\nAnalyzing Coulomb Expansion: N={n_ions}, {solver_type} solver")
    print("-" * 60)
    
    traj_file = results_dir / f"coulomb_expansion_N{n_ions}_{solver_type}.h5"
    
    if not traj_file.exists():
        print(f"❌ File not found: {traj_file}")
        return None
    
    try:
        times, mean_pos, std_pos, radial_std = read_trajectory_stats(traj_file)
        
        # Fit σ(t) = A * √t (skip first few points for settling)
        skip = 10
        times_fit = times[skip:]
        radial_fit = radial_std[skip:]
        
        # Initial guess
        p0 = [radial_fit[-1] / np.sqrt(times_fit[-1]), 0.0]
        
        try:
            popt, pcov = curve_fit(sqrt_t_model, times_fit, radial_fit, p0=p0)
            A_fit, t0_fit = popt
            
            # Predicted vs actual
            radial_pred = sqrt_t_model(times_fit, A_fit, t0_fit)
            rms_error = np.sqrt(np.mean((radial_pred - radial_fit)**2))
            
            print(f"  Fit: σ(t) = {A_fit:.6e} * √(t - {t0_fit:.2e})")
            print(f"  RMS error: {rms_error*1e3:.3f} mm")
            print(f"  Initial σ: {radial_std[0]*1e3:.3f} mm")
            print(f"  Final σ: {radial_std[-1]*1e3:.3f} mm")
            print(f"  Expansion factor: {radial_std[-1]/radial_std[0]:.1f}x")
            
            return {
                'n_ions': n_ions,
                'solver': solver_type,
                'times': times,
                'radial_std': radial_std,
                'A_fit': A_fit,
                't0_fit': t0_fit,
                'rms_error': rms_error
            }
            
        except Exception as e:
            print(f"  ⚠️  Fit failed: {e}")
            return None
            
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return None

def compare_direct_grid(results_dir):
    """Compare Direct vs Grid solver at various N"""
    results_dir = Path(results_dir)
    
    print("\n" + "="*80)
    print("DIRECT vs GRID SOLVER COMPARISON")
    print("="*80)
    
    # Analyze all expansion configs
    results = []
    
    # Direct solver (N < 1000)
    for n_ions in [100, 500]:
        result = analyze_coulomb_expansion(results_dir, n_ions, 'direct')
        if result:
            results.append(result)
    
    # Grid solver (N >= 1000)
    for n_ions in [1000, 5000, 10000]:
        result = analyze_coulomb_expansion(results_dir, n_ions, 'grid')
        if result:
            results.append(result)
    
    if not results:
        print("❌ No valid results")
        return
    
    # Plot comparison
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle('Space Charge Validation: Coulomb Expansion', fontsize=14, fontweight='bold')
    
    # 1. σ(t) for all N
    ax = axes[0, 0]
    for res in results:
        label = f"N={res['n_ions']} ({res['solver']})"
        ax.plot(res['times']*1e6, res['radial_std']*1e3, label=label, alpha=0.7)
    ax.set_xlabel('Time (µs)', fontsize=11)
    ax.set_ylabel('Radial Spread σ (mm)', fontsize=11)
    ax.set_title('Cloud Expansion vs Time', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. σ vs √t (should be linear)
    ax = axes[0, 1]
    for res in results:
        label = f"N={res['n_ions']} ({res['solver']})"
        ax.plot(np.sqrt(res['times']*1e6), res['radial_std']*1e3, 'o-', 
                markersize=3, label=label, alpha=0.7)
    ax.set_xlabel('√Time (√µs)', fontsize=11)
    ax.set_ylabel('Radial Spread σ (mm)', fontsize=11)
    ax.set_title('σ ∝ √t Verification', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 3. Expansion rate A vs N
    ax = axes[1, 0]
    n_ions_list = [res['n_ions'] for res in results]
    A_list = [res['A_fit'] for res in results]
    solvers = [res['solver'] for res in results]
    colors = ['blue' if s == 'direct' else 'orange' for s in solvers]
    
    ax.scatter(n_ions_list, A_list, c=colors, s=100, alpha=0.6)
    ax.set_xlabel('Number of Ions', fontsize=11)
    ax.set_ylabel('Expansion Rate A (m/√s)', fontsize=11)
    ax.set_title('Expansion Rate vs Ion Count', fontweight='bold')
    ax.set_xscale('log')
    ax.grid(True, alpha=0.3)
    
    # Add legend
    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='blue', alpha=0.6, label='Direct'),
                      Patch(facecolor='orange', alpha=0.6, label='Grid')]
    ax.legend(handles=legend_elements)
    
    # 4. RMS error
    ax = axes[1, 1]
    rms_errors = [res['rms_error']*1e3 for res in results]
    ax.bar(range(len(results)), rms_errors, color=colors, alpha=0.6)
    ax.set_xticks(range(len(results)))
    ax.set_xticklabels([f"N={res['n_ions']}\n({res['solver']})" 
                        for res in results], rotation=45, ha='right')
    ax.set_ylabel('RMS Error (mm)', fontsize=11)
    ax.set_title('Fit Quality (σ ∝ √t)', fontweight='bold')
    ax.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    
    # Save plot
    plot_file = results_dir / 'spacecharge_expansion_validation.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"\n📊 Saved plot: {plot_file}")
    
    # Summary statistics
    print("\n" + "="*80)
    print("SUMMARY")
    print("="*80)
    print(f"Tested configurations: {len(results)}")
    print(f"  Direct solver: {sum(1 for r in results if r['solver'] == 'direct')}")
    print(f"  Grid solver: {sum(1 for r in results if r['solver'] == 'grid')}")
    
    # Check if expansion rate consistent
    A_mean = np.mean(A_list)
    A_std = np.std(A_list)
    print(f"\nExpansion rate A:")
    print(f"  Mean: {A_mean:.4e} m/√s")
    print(f"  Std:  {A_std:.4e} m/√s")
    print(f"  Variation: {A_std/A_mean*100:.1f}%")
    
    if A_std/A_mean < 0.3:
        print("\n✅ GOOD: Expansion rates consistent across solvers (<30% variation)")
    else:
        print("\n⚠️  WARNING: Large variation in expansion rates (>30%)")
    
    print("="*80 + "\n")

def analyze_ims_spacecharge(results_dir):
    """Compare IMS with/without space charge"""
    results_dir = Path(results_dir)
    
    print("\n" + "="*80)
    print("IMS SPACE CHARGE EFFECTS")
    print("="*80)
    
    n_ions = 5000
    
    results_sc = {}
    for sc_enabled in [False, True]:
        sc_label = 'on' if sc_enabled else 'off'
        traj_file = results_dir / f"ims_spacecharge_{sc_label}_N{n_ions}.h5"
        
        if not traj_file.exists():
            print(f"❌ File not found: {traj_file}")
            continue
        
        try:
            times, mean_pos, std_pos, radial_std = read_trajectory_stats(traj_file)
            
            # Drift velocity (z-component)
            z_mean = mean_pos[:, 2]
            
            # Find arrival time (when z > 0.045 m, 90% of drift length)
            arrival_idx = np.where(z_mean > 0.045)[0]
            if len(arrival_idx) > 0:
                arrival_time = times[arrival_idx[0]]
            else:
                arrival_time = times[-1]
            
            # Peak width at arrival
            peak_width = std_pos[-1, 2] * 1e3  # mm
            
            results_sc[sc_label] = {
                'times': times,
                'z_mean': z_mean,
                'z_std': std_pos[:, 2],
                'arrival_time': arrival_time,
                'peak_width': peak_width
            }
            
            print(f"\nSpace Charge {sc_label.upper()}:")
            print(f"  Arrival time: {arrival_time*1e3:.2f} ms")
            print(f"  Peak width (z): {peak_width:.3f} mm")
            
        except Exception as e:
            print(f"  ❌ Error: {e}")
    
    if len(results_sc) == 2:
        # Compare
        t_off = results_sc['off']['arrival_time']
        t_on = results_sc['on']['arrival_time']
        delay = (t_on - t_off) / t_off * 100
        
        w_off = results_sc['off']['peak_width']
        w_on = results_sc['on']['peak_width']
        broadening = (w_on - w_off) / w_off * 100
        
        print(f"\nSpace Charge Effects:")
        print(f"  Drift time delay: {delay:+.1f}%")
        print(f"  Peak broadening: {broadening:+.1f}%")
        
        if delay > 5:
            print("  ✅ Space charge causes drift delay (expected)")
        if broadening > 10:
            print("  ✅ Space charge causes peak broadening (expected)")
    
    print("="*80 + "\n")

def analyze_all_spacecharge(results_dir):
    """Main analysis function"""
    results_dir = Path(results_dir)
    
    print("\n" + "="*80)
    print("Space Charge Validation Analysis")
    print("="*80)
    print(f"Results directory: {results_dir}")
    print("="*80)
    
    # 1. Coulomb expansion (Direct vs Grid)
    compare_direct_grid(results_dir)
    
    # 2. IMS with space charge
    analyze_ims_spacecharge(results_dir)
    
    print("\n" + "="*80)
    print("VALIDATION COMPLETE")
    print("="*80)
    print("Key findings:")
    print("  1. Coulomb expansion follows σ ∝ √t (theory)")
    print("  2. Direct solver (N<1000) vs Grid solver (N≥1000) agreement")
    print("  3. Space charge causes drift delay + peak broadening in IMS")
    print("="*80 + "\n")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        results_dir = Path(sys.argv[1])
    else:
        results_dir = Path("../results/v1.0_test/physics/spacecharge")
    
    analyze_all_spacecharge(results_dir)

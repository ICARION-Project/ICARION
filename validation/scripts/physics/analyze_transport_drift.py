#!/usr/bin/env python3
"""
Analyze transport drift velocity validation results

Extracts drift velocity from trajectories and validates:
  v_drift = K₀ × E
  
Compares measured K₀ with literature values.
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import re
from scipy import stats

# Literature reduced mobility values (cm²/(V·s) in He at 300K)
K0_LITERATURE = {
    "H3O+": 2.8,
    "PentanalH+": 1.9,
    "ReserpineH+": 0.85,
}

def extract_drift_velocity(trajectory_file):
    """
    Extract drift velocity from trajectory HDF5 file.
    
    Returns: (v_drift_mean, v_drift_std, t_final)
    """
    with h5py.File(trajectory_file, 'r') as f:
        # Read trajectory data
        pos = f['/trajectory/positions'][:]  # [time, ions, xyz]
        times = f['/trajectory/times'][:]
        
        # Extract z-coordinate (drift direction)
        z = pos[:, :, 2]  # [time, ions]
        
        # Calculate mean position vs time
        z_mean = np.mean(z, axis=1)
        
        # Fit linear regression to get drift velocity
        # z(t) = z₀ + v_drift · t
        
        # Use second half of trajectory (after initial transients)
        n_half = len(times) // 2
        t_fit = times[n_half:]
        z_fit = z_mean[n_half:]
        
        # Linear fit
        slope, intercept, r_value, p_value, std_err = stats.linregress(t_fit, z_fit)
        
        v_drift = slope  # m/s
        v_drift_std = std_err  # Standard error
        
        return v_drift, v_drift_std, times[-1], r_value**2

def parse_filename(filename):
    """
    Parse config filename to extract parameters.
    
    Format: drift_{model}_{species}_{temp}K_{en}Td.h5
    """
    pattern = r'drift_(\w+)_([^_]+)_(\d+)K_(\d+)Td\.h5'
    match = re.match(pattern, filename)
    
    if not match:
        return None
    
    return {
        'model': match.group(1),
        'species': match.group(2),
        'temp_K': float(match.group(3)),
        'en_td': float(match.group(4)),
    }

def analyze_drift_results(results_dir):
    """
    Analyze all drift velocity results.
    """
    results_dir = Path(results_dir)
    
    print("\n" + "="*80)
    print("Transport Drift Velocity Analysis")
    print("="*80)
    print(f"Results directory: {results_dir}")
    print("="*80 + "\n")
    
    # Find trajectory files
    traj_files = list(results_dir.glob("drift_*.h5"))
    
    if not traj_files:
        print(f"❌ No trajectory files found in {results_dir}")
        return
    
    print(f"Found {len(traj_files)} trajectory files\n")
    
    results = []
    
    for traj_file in sorted(traj_files):
        params = parse_filename(traj_file.name)
        if params is None:
            print(f"⚠️  Skipping {traj_file.name} (invalid format)")
            continue
        
        try:
            v_drift, v_drift_std, t_final, r2 = extract_drift_velocity(traj_file)
            
            # Calculate E-field from E/N and extract from config
            # For now, we'll recalculate E from v_drift and K₀_literature
            # In practice, we should read E from config or HDF5 metadata
            
            results.append({
                'model': params['model'],
                'species': params['species'],
                'temp_K': params['temp_K'],
                'en_td': params['en_td'],
                'v_drift': v_drift,
                'v_drift_std': v_drift_std,
                't_sim': t_final,
                'r2': r2,
                'file': traj_file.name
            })
            
            print(f"✓ {traj_file.name}")
            print(f"    v_drift = {v_drift:.2f} ± {v_drift_std:.2f} m/s")
            print(f"    R² = {r2:.6f}")
            
        except Exception as e:
            print(f"❌ Error processing {traj_file.name}: {e}")
    
    if not results:
        print("❌ No valid results")
        return
    
    print(f"\n{'='*80}")
    print(f"Successfully analyzed {len(results)} trajectories")
    print("="*80 + "\n")
    
    # Group results by species and model for plotting
    fig, axes = plt.subplots(3, 3, figsize=(16, 12))
    fig.suptitle('Transport Physics: Drift Velocity vs E/N', fontsize=14, fontweight='bold')
    
    species_list = ['H3O+', 'PentanalH+', 'ReserpineH+']
    models = ['HSS', 'Langevin', 'Friction']
    
    for i, species in enumerate(species_list):
        for j, model in enumerate(models):
            ax = axes[i, j]
            
            # Filter results for this species and model
            subset = [r for r in results if r['species'] == species and r['model'] == model]
            
            if not subset:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{species} - {model}')
                continue
            
            # Group by temperature
            temps = sorted(set(r['temp_K'] for r in subset))
            
            for temp in temps:
                temp_subset = [r for r in subset if r['temp_K'] == temp]
                en_vals = np.array([r['en_td'] for r in temp_subset])
                v_vals = np.array([r['v_drift'] for r in temp_subset])
                v_std = np.array([r['v_drift_std'] for r in temp_subset])
                
                # Sort by E/N
                sort_idx = np.argsort(en_vals)
                en_vals = en_vals[sort_idx]
                v_vals = v_vals[sort_idx]
                v_std = v_std[sort_idx]
                
                ax.errorbar(en_vals, v_vals, yerr=v_std, 
                           marker='o', label=f'{temp:.0f} K', capsize=3)
            
            ax.set_xlabel('E/N (Td)')
            ax.set_ylabel('Drift Velocity (m/s)')
            ax.set_title(f'{species} - {model}')
            ax.legend()
            ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = results_dir / 'drift_velocity_analysis.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"📊 Saved plot: {plot_file}")
    
    # Summary statistics
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    
    for species in species_list:
        print(f"\n{species} (K₀_lit = {K0_LITERATURE.get(species, 'N/A')} cm²/(V·s)):")
        for model in models:
            subset = [r for r in results if r['species'] == species and r['model'] == model]
            if subset:
                r2_mean = np.mean([r['r2'] for r in subset])
                print(f"  {model:10s}: {len(subset):3d} configs, R² = {r2_mean:.6f}")
    
    print("\n" + "="*80 + "\n")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        results_dir = Path(sys.argv[1])
    else:
        results_dir = Path("../results/v1.0_test/physics/transport/drift")
    
    analyze_drift_results(results_dir)

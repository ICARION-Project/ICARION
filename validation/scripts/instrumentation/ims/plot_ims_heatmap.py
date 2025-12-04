#!/usr/bin/env python3
"""
Create 2D heatmap of relative error vs E/N and pressure for IMS validation
"""

import csv
import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from pathlib import Path
import sys

# Physical constants
kB = 1.380649e-23  # J/K
amu = 1.66053906660e-27  # kg
K0_H3Op_He = 24.1e-4  # m²/(V·s)
M_He = 4.003 * amu  # kg
LOSCHMIDT = 2.6867774e25  # m^-3
P0 = 101325.0  # Pa
T0 = 273.15  # K

def analyze_single_file(h5_file):
    """Analyze a single HDF5 file"""
    with h5py.File(h5_file, 'r') as f:
        axial_V = f['domains/domain_0/fields/dc/axial_V'][()]
        length_m = f['domains/domain_0/geometry/length_m'][()]
        E_Vm = axial_V / length_m
        P_Pa = f['domains/domain_0/environment/pressure_Pa'][()]
        T_K = f['domains/domain_0/environment/temperature_K'][()]
        N_m3 = f['domains/domain_0/environment/particle_density_m3'][()]
        positions = f['/trajectory/positions'][:]
        times = f['/trajectory/time'][:]
    
    z_pos = positions[:, :, 2]
    drift_distance = z_pos[-1, :] - z_pos[0, :]
    drifted_mask = drift_distance > 1e-3
    n_drifted = int(np.sum(drifted_mask))
    
    if n_drifted == 0:
        return None
    
    t_total = float(times[-1] - times[0])
    mean_drift_distance = float(np.mean(drift_distance[drifted_mask]))
    v_drift_mean = mean_drift_distance / t_total
    
    E_N_Td = (E_Vm / N_m3) * 1e21
    N_ratio = LOSCHMIDT / N_m3
    v_drift_field = K0_H3Op_He * E_Vm * N_ratio
    T_heating = (M_He / (3 * kB)) * v_drift_field**2
    T_eff = T_K + T_heating
    K_eff = K0_H3Op_He * np.sqrt(T0 / T_eff)
    v_expected = K_eff * E_Vm * N_ratio
    
    if v_expected == 0 or E_Vm == 0:
        return None
    
    error_pct = 100 * (v_drift_mean - v_expected) / v_expected
    
    mean_drift_time_us = t_total * 1e6
    K_m2Vs = v_drift_mean / E_Vm
    K0_m2Vs = K_m2Vs * (P_Pa / P0) * (T0 / T_K)
    K0_cm2Vs = K0_m2Vs * 1e4
    
    return {
        'E_N_Td': float(E_N_Td),
        'P_Pa': float(P_Pa),
        'T_K': float(T_K),
        'N_m3': float(N_m3),
        'E_Vm': float(E_Vm),
        'num_crossings': n_drifted,
        'mean_drift_time_us': float(mean_drift_time_us),
        'drift_velocity_ms': float(v_drift_mean),
        'K_m2Vs': float(K_m2Vs),
        'K0_m2Vs': float(K0_m2Vs),
        'K0_cm2Vs': float(K0_cm2Vs),
        'error_pct': float(error_pct)
    }

def write_csv(rows, output_path):
    """Write IMS analysis results to CSV in benchmark-compatible format"""
    if not rows:
        print("No results to write to CSV.")
        return
    
    fieldnames = [
        'filename', 'collision_model', 'EN_Td', 'pressure_Pa', 'temperature_K', 'N_(1/m3)',
        'E_(V/m)', 'num_crossings', 'mean_drift_time_us',
        'drift_velocity_(m/s)', 'K_(m2/Vs)', 'K0_(m2/Vs)', 'K0_(cm2/Vs)'
    ]
    
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    rows_sorted = sorted(rows, key=lambda r: (r['EN_Td'], r['pressure_Pa'], r['filename']))
    
    with output_path.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows_sorted:
            writer.writerow(row)
    
    print(f"CSV results saved to: {output_path}")

def create_heatmap(results_dir, output_file, csv_output):
    """Create 2D heatmap of error vs E/N and pressure"""
    
    results_path = Path(results_dir)
    h5_files = list(results_path.glob('**/*.h5'))
    
    # Collect results by model
    results = {'HSS': [], 'EHSS': [], 'Friction': []}
    csv_rows = []
    
    for h5_file in h5_files:
        basename = h5_file.stem
        
        if 'ehss' in basename.lower():
            model = 'EHSS'
        elif 'hss' in basename.lower():
            model = 'HSS'
        elif 'friction' in basename.lower():
            model = 'Friction'
        else:
            continue
        
        try:
            data = analyze_single_file(str(h5_file))
            if data and 10 <= data['E_N_Td'] <= 50:  # Only 10-50 Td
                csv_rows.append({
                    'filename': h5_file.name,
                    'collision_model': model,
                    'EN_Td': data['E_N_Td'],
                    'pressure_Pa': data['P_Pa'],
                    'temperature_K': data['T_K'],
                    'N_(1/m3)': data['N_m3'],
                    'E_(V/m)': data['E_Vm'],
                    'num_crossings': data['num_crossings'],
                    'mean_drift_time_us': data['mean_drift_time_us'],
                    'drift_velocity_(m/s)': data['drift_velocity_ms'],
                    'K_(m2/Vs)': data['K_m2Vs'],
                    'K0_(m2/Vs)': data['K0_m2Vs'],
                    'K0_(cm2/Vs)': data['K0_cm2Vs']
                })
                
                # Round E/N to nearest 10 Td (10, 20, 30, 40, 50)
                en_rounded = round(data['E_N_Td'] / 10) * 10
                p_rounded = round(data['P_Pa'])
                
                # Update data with rounded values
                data['E_N_Td'] = en_rounded
                data['P_Pa'] = p_rounded
                
                # Check if we already have this combination
                if not any(d['E_N_Td'] == en_rounded and d['P_Pa'] == p_rounded 
                          for d in results[model]):
                    results[model].append(data)
        except:
            pass
    
    write_csv(csv_rows, csv_output)
    
    # Create figure with 3 subplots
    fig, axes = plt.subplots(1, 3, figsize=(18, 5))
    
    for idx, model in enumerate(['HSS', 'EHSS', 'Friction']):
        ax = axes[idx]
        
        if not results[model]:
            ax.text(0.5, 0.5, 'No data', ha='center', va='center')
            ax.set_title(model)
            continue
        
        # Get unique E/N and pressure values
        en_values = sorted(set(d['E_N_Td'] for d in results[model]))
        p_values = sorted(set(d['P_Pa'] for d in results[model]))
        
        # Create grid - use absolute error for coloring
        error_grid = np.full((len(en_values), len(p_values)), np.nan)
        error_grid_raw = np.full((len(en_values), len(p_values)), np.nan)
        
        for d in results[model]:
            en_idx = en_values.index(d['E_N_Td'])
            p_idx = p_values.index(d['P_Pa'])
            error_grid[en_idx, p_idx] = abs(d['error_pct'])  # Absolute value for color
            error_grid_raw[en_idx, p_idx] = d['error_pct']  # Keep raw for annotation
        
        # Plot heatmap with sequential colormap (green=good, red=bad)
        vmax = min(np.nanmax(error_grid), 100)  # Cap at 100%
        
        im = ax.imshow(error_grid, aspect='auto', cmap='RdYlGn_r',
                      vmin=0, vmax=vmax,
                      origin='lower', interpolation='nearest')
        
        # Set ticks
        ax.set_xticks(range(len(p_values)))
        ax.set_xticklabels([f'{int(p)}' for p in p_values])
        ax.set_yticks(range(len(en_values)))
        ax.set_yticklabels([f'{int(en)}' for en in en_values])
        
        # Add text annotations with signed values
        for i in range(len(en_values)):
            for j in range(len(p_values)):
                val_abs = error_grid[i, j]
                val_raw = error_grid_raw[i, j]
                if not np.isnan(val_abs):
                    color = 'white' if val_abs > 40 else 'black'
                    ax.text(j, i, f'{val_raw:+.0f}%', ha='center', va='center',
                           color=color, fontsize=9, weight='bold')
        
        ax.set_xlabel('Pressure (Pa)', fontsize=12)
        ax.set_ylabel('E/N (Td)', fontsize=12)
        ax.set_title(f'{model} Model', fontsize=14, weight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Relative Error (%)', fontsize=10)
    
    plt.suptitle('IMS Drift Velocity Validation: Relative Error vs E/N and Pressure\n' +
                 'H3O⁺ in He, 300K, with T_eff correction',
                 fontsize=16, weight='bold', y=1.02)
    
    plt.tight_layout()
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to: {output_file}")
    
    return fig

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: plot_ims_heatmap.py <results_dir>")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    output_file = 'results/v1.0_test/instruments/ims/ims_EN_heatmap.png'
    csv_output = Path(__file__).with_name('ims_results_ICARION.csv')
    
    create_heatmap(results_dir, output_file, csv_output)
    print("\nHeatmap Summary:")
    print("  Green: Good agreement (error near 0%)")
    print("  Yellow: Moderate error")
    print("  Red: Large error")

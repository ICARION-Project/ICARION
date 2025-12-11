#!/usr/bin/env python3
"""
Create 2D heatmap of relative error vs E/N and pressure for IMS validation
"""

import argparse
import csv
import h5py
import numpy as np
import matplotlib.pyplot as plt
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
        'T_eff_K': float(T_eff),
        'v_expected_ms': float(v_expected),
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
        'E_(V/m)', 'T_eff_K', 'drift_velocity_expected_(m/s)', 'num_crossings', 'mean_drift_time_us',
        'drift_velocity_(m/s)', 'K_(m2/Vs)', 'K0_(m2/Vs)', 'K0_(cm2/Vs)', 'error_pct'
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

def write_summary_csv(results_by_model, output_path):
    """Write per-model error summary for doc ingestion"""

    rows = []
    for model, data_list in results_by_model.items():
        for data in data_list:
            rows.append({
                'collision_model': model,
                'EN_Td': int(data['E_N_Td']),
                'pressure_Pa': int(data['P_Pa']),
                'T_eff_K': round(data['T_eff_K'], 2),
                'drift_velocity_(m/s)': round(data['drift_velocity_ms'], 2),
                'expected_velocity_(m/s)': round(data['v_expected_ms'], 2),
                'error_pct': round(data['error_pct'], 2)
            })

    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with output_path.open('w', newline='') as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                'collision_model', 'EN_Td', 'pressure_Pa', 'T_eff_K',
                'drift_velocity_(m/s)', 'expected_velocity_(m/s)', 'error_pct'
            ]
        )
        writer.writeheader()
        for row in sorted(rows, key=lambda r: (r['collision_model'], r['EN_Td'], r['pressure_Pa'])):
            writer.writerow(row)

    print(f"Summary CSV saved to: {output_path}")


def create_heatmap(results_dir, output_file, csv_output, summary_csv):
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
            if data and 0 <= data['E_N_Td'] <= 10:
                csv_rows.append({
                    'filename': h5_file.name,
                    'collision_model': model,
                    'EN_Td': data['E_N_Td'],
                    'pressure_Pa': data['P_Pa'],
                    'temperature_K': data['T_K'],
                    'N_(1/m3)': data['N_m3'],
                    'E_(V/m)': data['E_Vm'],
                    'T_eff_K': data['T_eff_K'],
                    'drift_velocity_expected_(m/s)': data['v_expected_ms'],
                    'num_crossings': data['num_crossings'],
                    'mean_drift_time_us': data['mean_drift_time_us'],
                    'drift_velocity_(m/s)': data['drift_velocity_ms'],
                    'K_(m2/Vs)': data['K_m2Vs'],
                    'K0_(m2/Vs)': data['K0_m2Vs'],
                    'K0_(cm2/Vs)': data['K0_cm2Vs'],
                    'error_pct': data['error_pct']
                })
                
                # Round E/N to nearest 1 Td
                en_rounded = round(data['E_N_Td']) 
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
    write_summary_csv(results, summary_csv)
    
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
        # Instead of plotting error, we want to plot mobilities now
        # Create grid - use mobility values
        mobility_grid = np.full((len(en_values), len(p_values)), np.nan)
        mobility_grid_raw = np.full((len(en_values), len(p_values)), np.nan)
        
        for d in results[model]:
            en_idx = en_values.index(d['E_N_Td'])
            p_idx = p_values.index(d['P_Pa'])
            mobility_grid[en_idx, p_idx] = abs(d['K0_cm2Vs'])  # Absolute value for color
            mobility_grid_raw[en_idx, p_idx] = d['K0_cm2Vs']  # Keep raw for annotation
        
        # Plot heatmap with sequential colormap (green=good, red=bad)
        vmax = np.nanmax(mobility_grid)  
        vmin = np.nanmin(mobility_grid)
        
        im = ax.imshow(mobility_grid, aspect='auto', cmap='RdYlGn_r',
                      vmin=vmin, vmax=vmax,
                      origin='lower', interpolation='nearest')
        
        # Set ticks
        ax.set_xticks(range(len(p_values)))
        ax.set_xticklabels([f'{int(p)}' for p in p_values])
        ax.set_yticks(range(len(en_values)))
        ax.set_yticklabels([f'{int(en)}' for en in en_values])
        
        # Add text annotations with signed values
        for i in range(len(en_values)):
            for j in range(len(p_values)):
                val_abs = mobility_grid[i, j]
                val_raw = mobility_grid_raw[i, j]
                if not np.isnan(val_abs):
                    color = 'white' if val_abs > 40 else 'black'
                    # plot the raw value without unit or sign and with 1 decimal place
                    ax.text(j, i, f'{val_raw:.1f}', ha='center', va='center',
                           color=color, fontsize=9, weight='bold')
        
        ax.set_xlabel('Pressure (Pa)', fontsize=12)
        ax.set_ylabel('E/N (Td)', fontsize=12)
        ax.set_title(f'{model} Model', fontsize=14, weight='bold')
        
        # Add colorbar
        cbar = plt.colorbar(im, ax=ax)
        cbar.set_label('Reduced Mobility (cm²/Vs)', fontsize=10)
        # scale colorbar axis limits to match data range from 10 to 25
        cbar_ticks = np.linspace(10, 25, num=6)
        cbar.set_ticks(cbar_ticks)
        cbar.set_ticklabels([f'{tick:.1f}' for tick in cbar_ticks])
    
    plt.suptitle('IMS Drift Velocity Validation: Reduced Mobility vs E/N and Pressure\n' +
                 'H3O⁺ in He, 300K',
                 fontsize=16, weight='bold', y=1.02)
    
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    plt.tight_layout()
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Heatmap saved to: {output_file}")
    
    return fig

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Generate IMS heatmap and rollups")
    parser.add_argument("results_dir", help="Directory containing IMS HDF5 outputs")
    parser.add_argument(
        "--figure-path",
        default="validation/figures/ims/ims_EN_heatmap.png",
        help="Destination for the heatmap figure"
    )
    parser.add_argument(
        "--csv",
        default="validation/results/v1.0_test/instruments/ims/ims_measurements.csv",
        help="Path for detailed measurement CSV"
    )
    parser.add_argument(
        "--summary-csv",
        default="validation/results/v1.0_test/instruments/ims/ims_error_summary.csv",
        help="Path for per-model error summary CSV"
    )

    args = parser.parse_args()

    create_heatmap(args.results_dir, args.figure_path, args.csv, args.summary_csv)
    print("\nHeatmap Summary:")
    print("  Green: Good agreement (error near 0%)")
    print("  Yellow: Moderate error")
    print("  Red: Large error")

#!/usr/bin/env python3
"""
Analyze IMS E/N mapping study and generate heatmap of relative errors
"""

import h5py
import sys
import os
from pathlib import Path

# Physical constants
kB = 1.380649e-23  # J/K - Boltzmann constant
amu = 1.66053906660e-27  # kg - atomic mass unit

# Mason-Schamp parameters
K0_H3Op_He = 24.1e-4  # m²/(V·s) - H3O+ in He reduced mobility
M_He = 4.003 * amu  # kg - He mass
LOSCHMIDT = 2.6867774e25  # m^-3 - number density at STP
T0 = 273.15  # K

def analyze_single_file(h5_file):
    """Analyze a single HDF5 file and return drift data"""
    
    with h5py.File(h5_file, 'r') as f:
        # Read simulation parameters
        axial_V = f['domains/domain_0/fields/dc/axial_V'][()]
        length_m = f['domains/domain_0/geometry/length_m'][()]
        E_Vm = axial_V / length_m
        P_Pa = f['domains/domain_0/environment/pressure_Pa'][()]
        T_K = f['domains/domain_0/environment/temperature_K'][()]
        N_m3 = f['domains/domain_0/environment/particle_density_m3'][()]
        
        # Read trajectory data
        positions = f['/trajectory/positions'][:]
        times = f['/trajectory/time'][:]
    
    # Get z-positions (drift direction)
    z_pos = positions[:, :, 2]
    z_start = z_pos[0, :]
    z_end = z_pos[-1, :]
    drift_distance = z_end - z_start
    
    # Filter ions that drifted > 1mm
    drifted_mask = drift_distance > 1e-3
    n_drifted = sum(drifted_mask)
    
    if n_drifted == 0:
        return None
    
    # Calculate measured drift velocity
    t_total = times[-1] - times[0]
    drift_velocities = drift_distance[drifted_mask] / t_total
    v_drift_mean = sum(drift_velocities) / len(drift_velocities)
    
    # Calculate E/N in Td
    E_N_Td = (E_Vm / N_m3) * 1e21
    
    # Calculate T_eff with He mass
    N_ratio = LOSCHMIDT / N_m3
    v_drift_field = K0_H3Op_He * E_Vm * N_ratio
    T_heating = (M_He / (3 * kB)) * v_drift_field**2
    T_eff = T_K + T_heating
    
    # Mason-Schamp with T_eff correction
    K_eff = K0_H3Op_He * (T0 / T_eff)**0.5
    v_expected = K_eff * E_Vm * N_ratio
    
    error_pct = 100 * (v_drift_mean - v_expected) / v_expected
    
    return {
        'E_N_Td': E_N_Td,
        'P_Pa': P_Pa,
        'T_K': T_K,
        'T_eff': T_eff,
        'v_meas': v_drift_mean,
        'v_expected': v_expected,
        'error_pct': error_pct,
        'n_ions': n_drifted
    }

def main(results_dir):
    """Analyze all HDF5 files in results directory"""
    
    results_path = Path(results_dir)
    
    # Find all HDF5 files
    h5_files = list(results_path.glob('**/*.h5'))
    
    if not h5_files:
        print(f"No HDF5 files found in {results_dir}")
        return
    
    print("="*80)
    print("IMS E/N Mapping Analysis")
    print("="*80)
    print(f"Found {len(h5_files)} HDF5 files")
    print()
    
    # Organize results by model, E/N, and pressure
    results = {}
    
    for h5_file in sorted(h5_files):
        basename = h5_file.stem
        
        # Extract model from filename
        if 'hss' in basename.lower():
            if 'ehss' in basename.lower():
                model = 'EHSS'
            else:
                model = 'HSS'
        elif 'friction' in basename.lower():
            model = 'Friction'
        else:
            continue
        
        try:
            data = analyze_single_file(str(h5_file))
            if data:
                if model not in results:
                    results[model] = []
                results[model].append(data)
        except Exception as e:
            print(f"Error analyzing {basename}: {e}")
    
    # Print results organized by model
    for model in sorted(results.keys()):
        print(f"\n{'='*80}")
        print(f"{model} Model")
        print('='*80)
        print(f"{'E/N (Td)':<10} {'P (Pa)':<10} {'T_eff (K)':<12} {'v_meas':<10} {'v_exp':<10} {'Error':<10}")
        print('-'*80)
        
        # Sort by E/N then pressure
        for data in sorted(results[model], key=lambda x: (x['E_N_Td'], x['P_Pa'])):
            print(f"{data['E_N_Td']:>8.1f}   "
                  f"{data['P_Pa']:>8.0f}   "
                  f"{data['T_eff']:>10.1f}   "
                  f"{data['v_meas']:>8.1f}   "
                  f"{data['v_expected']:>8.1f}   "
                  f"{data['error_pct']:>+7.1f}%")
    
    # Generate summary table by E/N and pressure
    print("\n" + "="*80)
    print("Summary: Error by E/N and Pressure")
    print("="*80)
    
    for model in sorted(results.keys()):
        print(f"\n{model}:")
        
        # Get unique E/N and pressure values
        en_values = sorted(set(d['E_N_Td'] for d in results[model]))
        p_values = sorted(set(d['P_Pa'] for d in results[model]))
        
        # Create lookup dict
        error_map = {}
        for d in results[model]:
            error_map[(d['E_N_Td'], d['P_Pa'])] = d['error_pct']
        
        # Print table
        print(f"{'E/N (Td)':<12}", end='')
        for p in p_values:
            print(f"{int(p):>10} Pa", end='')
        print()
        print('-'*70)
        
        for en in en_values:
            print(f"{en:>8.0f} Td  ", end='')
            for p in p_values:
                err = error_map.get((en, p), None)
                if err is not None:
                    status = '✅' if abs(err) < 10 else '⚠️'
                    print(f"{status} {err:>+6.1f}%", end=' ')
                else:
                    print("      ---", end=' ')
            print()
    
    print("\n" + "="*80)
    print("Analysis Complete")
    print("="*80)
    print("\nNote: ✅ = error < 10%, ⚠️ = error >= 10%")
    print(f"All results use T_eff correction with He mass")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: analyze_ims_EN_mapping.py <results_dir>")
        sys.exit(1)
    
    main(sys.argv[1])

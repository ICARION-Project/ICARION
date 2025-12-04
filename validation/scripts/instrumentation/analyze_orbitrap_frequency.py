#!/usr/bin/env python3
"""
Analyze Orbitrap axial oscillation frequency from trajectory data.

Theory:
  ω_z = sqrt(k * q / m)
  where k = field curvature constant (V/m²)
  
For hyperlogarithmic field: V(r,z) = (k/2)(z² - r²/2)
  k = 2*V_rad / (r_char² * ln(r_out/r_in))
  
Frequency:
  f_z = ω_z / (2π) = (1/2π) * sqrt(k * q / m)
"""

import h5py
import numpy as np
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import sys
from pathlib import Path

def calculate_theoretical_frequency(mass_u, charge_e, V_rad, r_char, r_in, r_out):
    """Calculate theoretical axial frequency for Orbitrap.
    
    Field curvature constant from hyperlogarithmic potential:
    k = 2*V_rad / (r_char² * ln(r_out/r_in) - 0.5*(r_out² - r_in²))
    
    This matches ElectricFieldForce.cpp line 337-339.
    """
    denominator = r_char**2 * np.log(r_out / r_in) - 0.5 * (r_out**2 - r_in**2)
    k = 2 * V_rad / denominator  # V/m²
    q_C = charge_e * 1.602176634e-19  # C
    m_kg = mass_u * 1.66053906660e-27  # kg
    omega_z = np.sqrt(abs(k) * q_C / m_kg)  # rad/s
    f_z = omega_z / (2 * np.pi)  # Hz
    return f_z, k

def analyze_trajectory(h5_path, config_path=None):
    """Extract axial oscillation frequency from HDF5 trajectory."""
    with h5py.File(h5_path, 'r') as f:
        # Read trajectory
        positions = f['trajectory/positions'][:]  # (n_steps, n_ions, 3)
        time = f['trajectory/time'][:]  # (n_steps,)
        
        # Take first ion if multiple
        if positions.ndim == 3:
            positions = positions[:, 0, :]
    
    # Read metadata from config file (HDF5 may not have complete metadata)
    if config_path is None:
        # Infer config path from h5 path
        config_path = str(h5_path).replace('.h5', '.json').replace('results/v1.0_test/instruments/', 'configs/instruments/')
    
    import json
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # Extract parameters from config
    species_id = config['ions']['species'][0].get('species_id', 'Unknown')
    
    # Read species database if species_id is given
    if 'species_database_path' in config:
        import json
        with open(config['species_database_path'], 'r') as f:
            species_db = json.load(f)
        # Database structure: {"species": {"H3O+": {...}, ...}}
        if 'species' in species_db and species_id in species_db['species']:
            species_data = species_db['species'][species_id]
            mass_u = species_data.get('mass_amu', species_data.get('mass_Da', 0))
            charge = species_data.get('charge', 1)
        else:
            mass_u = config['ions']['species'][0].get('mass_Da', config['ions']['species'][0].get('mass_amu', 0))
            charge = config['ions']['species'][0].get('charge', 1)
    else:
        mass_u = config['ions']['species'][0].get('mass_Da', config['ions']['species'][0].get('mass_amu', 0))
        charge = config['ions']['species'][0].get('charge', 1)
    
    # Read field parameters
    V_rad = config['domains'][0]['fields']['DC']['radial_V']
    r_char = config['domains'][0]['geometry']['radius_char_m']
    r_in = config['domains'][0]['geometry']['radius_in_m']
    r_out = config['domains'][0]['geometry']['radius_out_m']
    
    # Extract z-coordinate (axial oscillation)
    z = positions[:, 2]
    dt = time[1] - time[0]
    
    # Remove linear drift (ion may have initial z-velocity)
    z_detrended = z - np.linspace(z[0], z[-1], len(z))
    
    # FFT to find dominant frequency
    N = len(z_detrended)
    fft_vals = fft(z_detrended - np.mean(z_detrended))
    freqs = fftfreq(N, dt)
    
    # Only positive frequencies
    pos_mask = freqs > 0
    freqs_pos = freqs[pos_mask]
    fft_mag = np.abs(fft_vals[pos_mask])
    
    # Find peak frequency
    peaks, properties = find_peaks(fft_mag, height=np.max(fft_mag) * 0.1)
    if len(peaks) > 0:
        peak_idx = peaks[np.argmax(properties['peak_heights'])]
        f_measured = freqs_pos[peak_idx]
    else:
        # Fallback: just take max
        f_measured = freqs_pos[np.argmax(fft_mag)]
    
    # Calculate theoretical frequency
    f_theory, k = calculate_theoretical_frequency(mass_u, charge, V_rad, r_char, r_in, r_out)
    
    # Calculate error
    error_pct = 100 * (f_measured - f_theory) / f_theory
    
    # Calculate period
    T_measured = 1.0 / f_measured if f_measured > 0 else 0
    T_theory = 1.0 / f_theory if f_theory > 0 else 0
    
    return {
        'species': species_id,
        'mass_u': mass_u,
        'charge': charge,
        'V_rad': V_rad,
        'k': k,
        'f_theory_Hz': f_theory,
        'f_measured_Hz': f_measured,
        'T_theory_us': T_theory * 1e6,
        'T_measured_us': T_measured * 1e6,
        'error_pct': error_pct,
        'n_points': N,
        'duration_ms': time[-1] * 1e3
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_orbitrap_frequency.py <trajectory.h5> [<trajectory2.h5> ...]")
        sys.exit(1)
    
    results = []
    for h5_path in sys.argv[1:]:
        if not Path(h5_path).exists():
            print(f"⚠ File not found: {h5_path}")
            continue
        
        try:
            result = analyze_trajectory(h5_path)
            results.append(result)
        except Exception as e:
            print(f"✗ Error analyzing {h5_path}: {e}")
            continue
    
    if not results:
        print("No valid results to display")
        return
    
    # Print results table
    print("\n" + "="*100)
    print("ORBITRAP AXIAL FREQUENCY VALIDATION")
    print("="*100)
    print(f"{'Species':<15} {'Mass [u]':<10} {'V_rad [V]':<12} {'f_theory [Hz]':<15} {'f_meas [Hz]':<15} {'Error [%]':<10} {'Status':<10}")
    print("-"*100)
    
    for r in results:
        status = "✅ PASS" if abs(r['error_pct']) < 1.0 else "⚠ WARN" if abs(r['error_pct']) < 5.0 else "✗ FAIL"
        print(f"{r['species']:<15} {r['mass_u']:<10.2f} {r['V_rad']:<12.1f} {r['f_theory_Hz']:<15.2f} {r['f_measured_Hz']:<15.2f} {r['error_pct']:>+9.2f} {status:<10}")
    
    print("-"*100)
    print(f"\nField curvature constant k = {results[0]['k']:.2e} V/m²")
    print(f"Theory: f ∝ 1/√m  →  Frequency ratio test:")
    
    # Verify mass scaling: f1/f2 = sqrt(m2/m1)
    if len(results) >= 2:
        for i in range(len(results)-1):
            r1, r2 = results[i], results[i+1]
            ratio_measured = r1['f_measured_Hz'] / r2['f_measured_Hz']
            ratio_theory = np.sqrt(r2['mass_u'] / r1['mass_u'])
            ratio_error = 100 * (ratio_measured - ratio_theory) / ratio_theory
            print(f"  {r1['species']}/{r2['species']}: measured={ratio_measured:.4f}, theory={ratio_theory:.4f}, Δ={ratio_error:+.2f}%")
    
    print("="*100)

if __name__ == '__main__':
    main()

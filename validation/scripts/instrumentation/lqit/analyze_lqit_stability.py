#!/usr/bin/env python3
"""
Analyze LQIT stability and secular motion

Checks ion confinement, measures secular frequency, and validates Mathieu stability.
"""

import h5py
import numpy as np
import sys
from pathlib import Path

def analyze_lqit_stability(h5_file):
    """Analyze LQIT ion motion and stability"""
    
    with h5py.File(h5_file, 'r') as f:
        # Read simulation parameters
        r0_m = f['domains/domain_0/geometry/radius_m'][()]
        length_m = f['domains/domain_0/geometry/length_m'][()]
        origin_m = f['domains/domain_0/geometry/origin_m'][:]
        rf_freq_Hz = f['domains/domain_0/fields/rf/frequency_Hz'][()]
        rf_V = f['domains/domain_0/fields/rf/voltage_V'][()]
        
        # Check if DC is set
        if 'quad_V' in f['domains/domain_0/fields/dc']:
            dc_V = f['domains/domain_0/fields/dc/quad_V'][()]
        else:
            dc_V = 0.0
        
        # Read trajectory data
        positions = f['/trajectory/positions'][:]  # (time, ion, xyz)
        times = f['/trajectory/time'][:]
        
        # Get ion parameters from metadata
        species = f['metadata/species/names'][0].decode('utf-8')
        mass_kg = f['metadata/species/mass_kg'][0]
        mass_amu = mass_kg / 1.66054e-27
    
    # Calculate Mathieu parameters
    m_kg = mass_amu * 1.66054e-27
    q_C = 1.602176634e-19
    omega = 2 * np.pi * rf_freq_Hz
    
    q_param = (4 * q_C * rf_V) / (m_kg * omega**2 * r0_m**2)
    a_param = (-8 * q_C * dc_V) / (m_kg * omega**2 * r0_m**2)
    
    # Expected secular frequency
    omega_sec_theory = (q_param * omega) / (2 * np.sqrt(2))
    f_sec_theory = omega_sec_theory / (2 * np.pi)
    
    # Analyze ion positions
    n_times, n_ions, _ = positions.shape
    
    # Extract radial positions (sqrt(x^2 + y^2))
    x = positions[:, :, 0]
    y = positions[:, :, 1]
    z = positions[:, :, 2]
    r = np.sqrt(x**2 + y**2)
    
    # Check confinement; origin defines one end of the trap, not the midpoint
    z_start = origin_m[2]
    z_end = z_start + length_m
    z_min_trap = min(z_start, z_end)
    z_max_trap = max(z_start, z_end)
    
    r_max = r.max(axis=0)  # Max radius for each ion
    z_min = z.min(axis=0)  # Min z for each ion
    z_max = z.max(axis=0)  # Max z for each ion
    
    confined_radial = r_max < r0_m
    confined_axial = (z_min >= z_min_trap) & (z_max <= z_max_trap)
    confined = confined_radial & confined_axial
    
    n_confined = np.sum(confined)
    confinement_pct = 100 * n_confined / n_ions
    
    # Measure secular frequency from confined ions
    if n_confined > 0:
        # Take first confined ion
        ion_idx = np.where(confined)[0][0]
        
        # Use x-component (NOT r!) to avoid frequency doubling
        # r² = x² + y² oscillates at 2*f_sec for harmonic motion
        x_trajectory = x[:, ion_idx]
        
        # FFT to find dominant frequency
        dt = times[1] - times[0]
        freqs = np.fft.rfftfreq(len(x_trajectory), dt)
        fft = np.abs(np.fft.rfft(x_trajectory))
        
        # Find peak (ignore DC component at freq=0)
        peak_idx = np.argmax(fft[1:]) + 1
        f_sec_measured = freqs[peak_idx]
        
        freq_error = ((f_sec_measured - f_sec_theory) / f_sec_theory) * 100
    else:
        f_sec_measured = None
        freq_error = None
    
    # Stability classification
    stable_theory = (q_param < 0.908) and (abs(a_param) < 0.2)
    status_theory = "STABLE" if stable_theory else "UNSTABLE"
    
    stable_measured = confinement_pct > 90
    status_measured = "CONFINED" if stable_measured else "EJECTED"
    
    # Print results
    print(f"  Trap: r₀={r0_m*1e3:.1f} mm, L={length_m*1e3:.1f} mm")
    print(f"  RF: f={rf_freq_Hz/1e6:.2f} MHz, V={rf_V:.2f} V")
    if dc_V != 0:
        print(f"  DC: V={dc_V:.2f} V")
    print(f"  Ion: {species} ({mass_amu:.2f} u)")
    print()
    print(f"  Mathieu parameters:")
    print(f"    q = {q_param:.3f}")
    print(f"    a = {a_param:.3f}")
    print(f"    Theory: {status_theory}")
    print()
    print(f"  Confinement:")
    print(f"    Ions confined: {n_confined}/{n_ions} ({confinement_pct:.1f}%)")
    print(f"    Status: {status_measured}")
    
    if f_sec_measured is not None:
        print()
        print(f"  Secular frequency:")
        print(f"    Theory:   {f_sec_theory/1e3:.1f} kHz")
        print(f"    Measured: {f_sec_measured/1e3:.1f} kHz")
        print(f"    Error:    {freq_error:+.1f}%")
    
    # Return summary
    return {
        'q': q_param,
        'a': a_param,
        'stable_theory': stable_theory,
        'confinement_pct': confinement_pct,
        'f_sec_theory': f_sec_theory,
        'f_sec_measured': f_sec_measured,
        'freq_error': freq_error
    }

def main():
    if len(sys.argv) < 2:
        print("Usage: analyze_lqit_stability.py <h5_file1> [h5_file2] ...")
        sys.exit(1)
    
    print("=" * 70)
    print("LQIT Stability Analysis")
    print("=" * 70)
    
    results = []
    for h5_file in sys.argv[1:]:
        # Extract model/config from filename
        basename = Path(h5_file).stem
        parts = basename.split('_')
        model = parts[1].upper() if len(parts) > 1 else "UNKNOWN"
        
        print(f"\n{basename}:")
        print("-" * 70)
        
        try:
            result = analyze_lqit_stability(h5_file)
            result['model'] = model
            result['file'] = basename
            results.append(result)
        except Exception as e:
            print(f"  ❌ Error: {e}")
            import traceback
            traceback.print_exc()
    
    # Summary table
    if results:
        print("\n" + "=" * 70)
        print("SUMMARY")
        print("=" * 70)
        print(f"{'Model':<10} {'q':>6} {'a':>6} {'Theory':>8} {'Confined':>10} {'f_sec_err':>10}")
        print("-" * 70)
        
        for r in results:
            theory_status = "STABLE" if r['stable_theory'] else "UNSTABLE"
            confined_str = f"{r['confinement_pct']:.0f}%"
            
            if r['freq_error'] is not None:
                freq_str = f"{r['freq_error']:+.1f}%"
            else:
                freq_str = "N/A"
            
            # Status emoji
            if r['stable_theory'] and r['confinement_pct'] > 90:
                status = "✅"
            elif not r['stable_theory'] and r['confinement_pct'] < 10:
                status = "✅"
            else:
                status = "⚠️"
            
            print(f"{status} {r['model']:<8} {r['q']:6.3f} {r['a']:6.3f} {theory_status:>8} {confined_str:>10} {freq_str:>10}")

if __name__ == "__main__":
    main()

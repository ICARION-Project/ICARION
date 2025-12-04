#!/usr/bin/env python3
"""
Analyze FT-ICR cyclotron frequencies

Validates:
  f_c = q*B / (2π*m)
  
Measures cyclotron frequency from XY trajectory using FFT.
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import fft
from scipy.signal import find_peaks

# Constants
Q_E = 1.602176634e-19  # Elementary charge
AMU_TO_KG = 1.66054e-27

# Species mass lookup
SPECIES_MASSES = {
    "H3O+": 19.0,
    "PentanalH+": 87.0,
    "CaffeineH+": 195.0,
    "ReserpineH+": 609.0,
}

def measure_cyclotron_frequency(trajectory_file, species_id=None):
    """
    Measure cyclotron frequency from XY plane trajectory.
    
    For multi-species, species_id specifies which species to analyze.
    Returns: (freq_Hz, radius_m, period_s)
    """
    with h5py.File(trajectory_file, 'r') as f:
        # Read trajectory data
        pos = f['/trajectory/positions'][:]  # [time, ions, xyz]
        times = f['/trajectory/time'][:]
        
        # Get species IDs if multi-species
        if '/trajectory/species_ids' in f:
            species_ids = f['/trajectory/species_ids'][:]
            if species_id is not None:
                # Filter for specific species
                mask = species_ids == species_id.encode()
                if not np.any(mask):
                    raise ValueError(f"Species {species_id} not found")
                pos = pos[:, mask, :]
        
        # Extract X and Y coordinates
        x = pos[:, :, 0]  # [time, ions]
        y = pos[:, :, 1]
        
        # Average over ions
        x_avg = np.mean(x, axis=1)
        y_avg = np.mean(y, axis=1)
        
        # Compute radial distance (cyclotron radius)
        r = np.sqrt(x_avg**2 + y_avg**2)
        r_mean = np.mean(r[len(r)//2:])  # Mean radius (after settling)
        
        # FFT on X coordinate to find cyclotron frequency
        dt = times[1] - times[0]
        n = len(x_avg)
        freqs = fft.fftfreq(n, dt)
        fft_x = fft.fft(x_avg - np.mean(x_avg))  # Remove DC offset
        psd = np.abs(fft_x)**2
        
        # Only positive frequencies
        pos_mask = freqs > 0
        freqs_pos = freqs[pos_mask]
        psd_pos = psd[pos_mask]
        
        # Find peak frequency (cyclotron frequency)
        peak_idx = np.argmax(psd_pos)
        f_c = freqs_pos[peak_idx]
        
        # Estimate period
        T_c = 1.0 / f_c if f_c > 0 else 0.0
        
        return f_c, r_mean, T_c, freqs_pos, psd_pos

def analyze_fticr_results(results_dir, b_field_T):
    """
    Analyze FT-ICR results: measure cyclotron frequencies and validate scaling.
    """
    results_dir = Path(results_dir)
    
    print("\n" + "="*80)
    print("FT-ICR Cyclotron Frequency Analysis")
    print("="*80)
    print(f"Magnetic field: {b_field_T:.1f} T")
    print(f"Results directory: {results_dir}")
    print("="*80 + "\n")
    
    # Find trajectory files
    pattern = f"fticr_*_B{b_field_T:.1f}T.h5"
    traj_files = list(results_dir.glob(pattern))
    
    if not traj_files:
        print(f"❌ No trajectory files matching {pattern}")
        return
    
    print(f"Found {len(traj_files)} trajectory files\n")
    
    results = []
    
    for traj_file in sorted(traj_files):
        print(f"Analyzing: {traj_file.name}")
        
        # Extract species from filename
        if "multi_species" in traj_file.name:
            # Analyze each species separately
            for species_id, mass_amu in SPECIES_MASSES.items():
                print(f"  Species: {species_id} (m={mass_amu:.1f} u)")
                try:
                    f_c, r_c, T_c, freqs, psd = measure_cyclotron_frequency(
                        traj_file, species_id=species_id
                    )
                    
                    # Theoretical frequency
                    m_kg = mass_amu * AMU_TO_KG
                    f_c_theory = (Q_E * b_field_T) / (2 * np.pi * m_kg)
                    
                    # Error
                    error_pct = 100 * (f_c - f_c_theory) / f_c_theory
                    
                    print(f"    Measured f_c: {f_c/1e6:.4f} MHz")
                    print(f"    Theory f_c:   {f_c_theory/1e6:.4f} MHz")
                    print(f"    Error:        {error_pct:+.2f}%")
                    print(f"    Radius:       {r_c*1e3:.3f} mm")
                    print(f"    Period:       {T_c*1e6:.3f} µs")
                    
                    results.append({
                        'species': species_id,
                        'mass_amu': mass_amu,
                        'f_c_measured': f_c,
                        'f_c_theory': f_c_theory,
                        'error_pct': error_pct,
                        'radius_m': r_c,
                        'period_s': T_c,
                        'freqs': freqs,
                        'psd': psd
                    })
                    
                except Exception as e:
                    print(f"    ⚠️  Error: {e}")
        else:
            # Single species file
            for species_id, mass_amu in SPECIES_MASSES.items():
                if species_id in traj_file.name:
                    print(f"  Species: {species_id} (m={mass_amu:.1f} u)")
                    try:
                        f_c, r_c, T_c, freqs, psd = measure_cyclotron_frequency(traj_file)
                        
                        # Theoretical frequency
                        m_kg = mass_amu * AMU_TO_KG
                        f_c_theory = (Q_E * b_field_T) / (2 * np.pi * m_kg)
                        
                        # Error
                        error_pct = 100 * (f_c - f_c_theory) / f_c_theory
                        
                        print(f"    Measured f_c: {f_c/1e6:.4f} MHz")
                        print(f"    Theory f_c:   {f_c_theory/1e6:.4f} MHz")
                        print(f"    Error:        {error_pct:+.2f}%")
                        print(f"    Radius:       {r_c*1e3:.3f} mm")
                        print(f"    Period:       {T_c*1e6:.3f} µs")
                        
                        results.append({
                            'species': species_id,
                            'mass_amu': mass_amu,
                            'f_c_measured': f_c,
                            'f_c_theory': f_c_theory,
                            'error_pct': error_pct,
                            'radius_m': r_c,
                            'period_s': T_c,
                            'freqs': freqs,
                            'psd': psd
                        })
                        
                    except Exception as e:
                        print(f"    ⚠️  Error: {e}")
                    break
        print()
    
    if not results:
        print("❌ No valid results")
        return
    
    # Plot results
    fig, axes = plt.subplots(2, 2, figsize=(14, 10))
    fig.suptitle(f'FT-ICR Cyclotron Frequency Validation (B = {b_field_T:.1f} T)', 
                 fontsize=14, fontweight='bold')
    
    # Sort by mass
    results.sort(key=lambda r: r['mass_amu'])
    
    masses = np.array([r['mass_amu'] for r in results])
    f_c_meas = np.array([r['f_c_measured'] for r in results])
    f_c_theo = np.array([r['f_c_theory'] for r in results])
    errors = np.array([r['error_pct'] for r in results])
    species_names = [r['species'] for r in results]
    
    # 1. f_c vs m (should be 1/m)
    ax = axes[0, 0]
    ax.plot(masses, f_c_theo / 1e6, 'k--', label='Theory: f_c = qB/(2πm)', linewidth=2)
    ax.plot(masses, f_c_meas / 1e6, 'o-', markersize=8, label='Measured')
    ax.set_xlabel('Mass (u)', fontsize=11)
    ax.set_ylabel('Cyclotron Frequency (MHz)', fontsize=11)
    ax.set_title('Cyclotron Frequency vs Mass', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    # 2. Log-log plot (should be linear with slope -1)
    ax = axes[0, 1]
    ax.loglog(masses, f_c_theo / 1e6, 'k--', label='Theory: f_c ∝ 1/m', linewidth=2)
    ax.loglog(masses, f_c_meas / 1e6, 'o-', markersize=8, label='Measured')
    ax.set_xlabel('Mass (u)', fontsize=11)
    ax.set_ylabel('Cyclotron Frequency (MHz)', fontsize=11)
    ax.set_title('Log-Log: f_c ∝ m⁻¹ Scaling', fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3, which='both')
    
    # Fit slope
    log_m = np.log10(masses)
    log_fc = np.log10(f_c_meas / 1e6)
    fit = np.polyfit(log_m, log_fc, 1)
    slope = fit[0]
    ax.text(0.05, 0.95, f'Slope: {slope:.3f}\n(Expected: -1.000)',
            transform=ax.transAxes, verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    # 3. Error distribution
    ax = axes[1, 0]
    x_pos = np.arange(len(species_names))
    colors = ['green' if abs(e) < 1 else 'orange' if abs(e) < 5 else 'red' 
              for e in errors]
    ax.bar(x_pos, errors, color=colors, alpha=0.6)
    ax.axhline(0, color='k', linestyle='-', linewidth=0.8)
    ax.axhline(1, color='gray', linestyle='--', linewidth=0.8)
    ax.axhline(-1, color='gray', linestyle='--', linewidth=0.8)
    ax.set_xticks(x_pos)
    ax.set_xticklabels(species_names, rotation=45, ha='right')
    ax.set_ylabel('Error (%)', fontsize=11)
    ax.set_title('Frequency Measurement Error', fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    # Statistics text
    mean_err = np.mean(np.abs(errors))
    max_err = np.max(np.abs(errors))
    ax.text(0.95, 0.95, f'|Error|:\nMean: {mean_err:.2f}%\nMax: {max_err:.2f}%',
            transform=ax.transAxes, verticalalignment='top',
            horizontalalignment='right',
            bbox=dict(boxstyle='round', facecolor='lightblue', alpha=0.5))
    
    # 4. Example FFT spectrum (first species)
    ax = axes[1, 1]
    example = results[0]
    freqs_mhz = example['freqs'] / 1e6
    # Plot only up to 10 MHz
    mask = freqs_mhz < 10
    ax.semilogy(freqs_mhz[mask], example['psd'][mask])
    ax.axvline(example['f_c_measured'] / 1e6, color='r', linestyle='--', 
               label=f"f_c = {example['f_c_measured']/1e6:.3f} MHz")
    ax.set_xlabel('Frequency (MHz)', fontsize=11)
    ax.set_ylabel('Power Spectral Density', fontsize=11)
    ax.set_title(f"Example FFT: {example['species']}", fontweight='bold')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = results_dir / f'fticr_validation_B{b_field_T:.1f}T.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"📊 Saved plot: {plot_file}")
    
    # Summary
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    print(f"Tested species: {len(results)}")
    print(f"Mean |error|:   {mean_err:.2f}%")
    print(f"Max |error|:    {max_err:.2f}%")
    print(f"Scaling slope:  {slope:.3f} (expected: -1.000)")
    
    if mean_err < 1.0:
        print("\n✅ EXCELLENT: Mean error < 1%")
    elif mean_err < 5.0:
        print("\n✅ GOOD: Mean error < 5%")
    else:
        print("\n⚠️  WARNING: Mean error > 5%")
    
    if abs(slope + 1.0) < 0.05:
        print("✅ EXCELLENT: f_c ∝ m⁻¹ scaling verified (slope ≈ -1)")
    elif abs(slope + 1.0) < 0.1:
        print("✅ GOOD: f_c ∝ m⁻¹ scaling reasonable")
    else:
        print("⚠️  WARNING: Scaling deviates from theory")
    
    print("="*80 + "\n")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        results_dir = Path(sys.argv[1])
    else:
        results_dir = Path("../results/v1.0_test/instruments/fticr")
    
    # Magnetic field from configs
    b_field_T = 7.0
    
    analyze_fticr_results(results_dir, b_field_T)

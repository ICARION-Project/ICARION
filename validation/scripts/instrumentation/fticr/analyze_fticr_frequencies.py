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
import sys
from scipy import fft
from scipy.signal import find_peaks

# Allow importing shared helpers (species IDs, embedded inputs)
COMMON_DIR = Path(__file__).resolve().parents[2] / "common"
if str(COMMON_DIR) not in sys.path:
    sys.path.append(str(COMMON_DIR))
from hdf5_utils import load_species_ids  # noqa: E402

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

_IRREGULAR_SAMPLING_WARNED = set()


def _enforce_uniform_sampling(times, positions, trajectory_path):
    """Detect irregular timestamps and decimate to a uniform grid for FFTs."""
    if times.size < 3:
        return times, positions

    diffs = np.diff(times)
    median_dt = np.median(diffs)
    if median_dt <= 0:
        return times, positions

    jitter = np.median(np.abs(diffs - median_dt))
    if jitter <= 0.05 * median_dt:
        return times, positions

    threshold = 0.5 * median_dt
    keep = [0]
    last_t = times[0]
    for idx in range(1, len(times)):
        if times[idx] - last_t >= threshold:
            keep.append(idx)
            last_t = times[idx]

    if len(keep) < 3:
        return times, positions

    keep_idx = np.array(keep, dtype=int)
    key = str(trajectory_path)
    if key not in _IRREGULAR_SAMPLING_WARNED:
        print(
            f"    ⚠️  Detected irregular sampling in {Path(trajectory_path).name}: "
            f"{len(times)} snapshots → {keep_idx.size}, median Δt={median_dt*1e9:.3f} ns"
        )
        _IRREGULAR_SAMPLING_WARNED.add(key)

    return times[keep_idx], positions[keep_idx]

def _load_species_labels(h5_file):
    """Return 1D array of species labels for each ion in the trajectory."""
    species = load_species_ids(h5_file)
    # Use the first time slice if 2D (time, ion)
    if species.ndim == 2:
        return species[0, :]
    return species


def measure_cyclotron_frequency(trajectory_file, species_id=None, expected_freq=None, search_window=0.3):
    """
    Measure cyclotron frequency from XY plane trajectory.
    
    For multi-species, species_id specifies which species to analyze.
    Returns: (freq_Hz, radius_m, period_s)
    """
    with h5py.File(trajectory_file, 'r') as f:
        # Read trajectory data and enforce uniform sampling for FFT stability
        pos = f['/trajectory/positions'][:]  # [time, ions, xyz]
        times = f['/trajectory/time'][:]
        times, pos = _enforce_uniform_sampling(times, pos, trajectory_file)
        
        # Get species IDs if multi-species
        multi_species = species_id is not None
        if multi_species:
            species_labels = _load_species_labels(f)
            mask = species_labels == species_id
            if not np.any(mask):
                raise ValueError(f"Species {species_id} not found")
            pos = pos[:, mask, :]
        
        # Extract X and Y coordinates
        x = pos[:, :, 0]  # [time, ions]
        y = pos[:, :, 1]

        if times.size < 2:
            raise ValueError("Trajectory file has fewer than two time samples")

        diffs = np.diff(times)
        if np.any(diffs <= 0):
            raise ValueError("Trajectory time axis is not strictly increasing")
        dt = np.median(diffs)

        if multi_species:
            # Use whichever axial projection carries the strongest signal to avoid
            # phase-cancellation when combining X/Y directly.
            x_centered = x - np.mean(x, axis=0)
            y_centered = y - np.mean(y, axis=0)
            rms_x = np.sqrt(np.mean(x_centered**2, axis=0))
            rms_y = np.sqrt(np.mean(y_centered**2, axis=0))
            if np.max(rms_x) >= np.max(rms_y):
                best_idx = int(np.argmax(rms_x))
                samples = x_centered[:, best_idx]
            else:
                best_idx = int(np.argmax(rms_y))
                samples = y_centered[:, best_idx]
        else:
            x_avg = np.mean(x, axis=1)
            samples = x_avg - np.mean(x_avg)

        n = samples.shape[0]
        freqs = fft.fftfreq(n, dt)
        fft_vals = fft.fft(samples)
        psd = np.abs(fft_vals)**2

        pos_mask = freqs > 0
        freqs_pos = freqs[pos_mask]
        psd_pos = psd[pos_mask]

        if psd_pos.size == 0:
            raise ValueError("Unable to compute PSD peak for cyclotron frequency")

        search_mask = np.ones_like(freqs_pos, dtype=bool)
        if expected_freq is not None:
            lower = max(expected_freq * (1 - search_window), freqs_pos[0])
            upper = min(expected_freq * (1 + search_window), freqs_pos[-1])
            mask = (freqs_pos >= lower) & (freqs_pos <= upper)
            if np.any(mask):
                search_mask = mask
            else:
                # Fallback to alias-aware window if expected frequency exceeds Nyquist
                fs = 1.0 / dt
                alias_freq = abs(expected_freq - np.round(expected_freq / fs) * fs)
                lower = max(alias_freq * (1 - search_window), freqs_pos[0])
                upper = min(alias_freq * (1 + search_window), freqs_pos[-1])
                mask = (freqs_pos >= lower) & (freqs_pos <= upper)
                if np.any(mask):
                    search_mask = mask

        search_freqs = freqs_pos[search_mask]
        search_psd = psd_pos[search_mask]

        peak_idx = np.argmax(search_psd)
        f_c = search_freqs[peak_idx]
        # Compute radial distance (mean over latter half of trajectory, averaged across ions)
        r = np.sqrt(x**2 + y**2)
        r_mean = np.mean(r[len(r)//2:, :])
        
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
                    # Theoretical frequency used to guide peak search
                    m_kg = mass_amu * AMU_TO_KG
                    f_c_theory = (Q_E * b_field_T) / (2 * np.pi * m_kg)
                    f_c, r_c, T_c, freqs, psd = measure_cyclotron_frequency(
                        traj_file,
                        species_id=species_id,
                        expected_freq=f_c_theory
                    )
                    
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
                        m_kg = mass_amu * AMU_TO_KG
                        f_c_theory = (Q_E * b_field_T) / (2 * np.pi * m_kg)
                        f_c, r_c, T_c, freqs, psd = measure_cyclotron_frequency(
                            traj_file,
                            expected_freq=f_c_theory
                        )
                        
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

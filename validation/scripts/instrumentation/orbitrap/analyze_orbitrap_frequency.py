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

import argparse
import json
import h5py
import numpy as np
from scipy.fft import fft, fftfreq
from pathlib import Path
import matplotlib.pyplot as plt

AXIS_INDICES = {'x': 0, 'y': 1, 'z': 2}
DEFAULT_AXIAL_AXIS = 'z'
BAND_PROMINENCE_RATIO = 0.1

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


def _prepare_signal(signal):
    """Remove linear drift and DC offset, then apply a Hann window."""
    n = len(signal)
    if n < 2:
        raise ValueError("Signal must contain at least two samples")
    trend = np.linspace(signal[0], signal[-1], n)
    detrended = signal - trend
    detrended -= np.mean(detrended)
    window = np.hanning(n)
    if np.all(window == 0):
        return detrended
    return detrended * window


def _measure_axis_frequency(signal, time, f_theory):
    """Measure dominant frequency for a detrended axis signal."""
    cleaned = _prepare_signal(signal)
    dt = time[1] - time[0]
    fft_vals = fft(cleaned)
    freqs = fftfreq(len(cleaned), dt)
    pos_mask = freqs > 0
    freqs_pos = freqs[pos_mask]
    fft_mag = np.abs(fft_vals[pos_mask])
    if len(freqs_pos) == 0:
        raise ValueError("Unable to compute positive frequency components")

    if f_theory > 0:
        band_mask = (freqs_pos >= f_theory * 0.25) & (freqs_pos <= f_theory * 4.0)
        if np.any(band_mask):
            freqs_band = freqs_pos[band_mask]
            fft_band = fft_mag[band_mask]
            peak_idx = np.argmax(fft_band)
            if fft_band[peak_idx] >= BAND_PROMINENCE_RATIO * np.max(fft_mag):
                return freqs_band[peak_idx], fft_band[peak_idx]

    peak_idx = np.argmax(fft_mag)
    return freqs_pos[peak_idx], fft_mag[peak_idx]
def _normalize_species_label(label):
    if isinstance(label, bytes):
        return label.decode('utf-8')
    if isinstance(label, np.bytes_):
        return label.decode('utf-8')
    return str(label)


def _resolve_species_metadata(spec_cfg, species_db):
    species_id = spec_cfg.get('species_id') or spec_cfg.get('id') or spec_cfg.get('name')
    if not species_id:
        raise ValueError("Species entry missing 'species_id'")
    db_entry = species_db.get(species_id, {})
    mass_u = db_entry.get('mass_amu', db_entry.get('mass_Da'))
    if mass_u is None:
        mass_u = spec_cfg.get('mass_Da', spec_cfg.get('mass_amu', 0))
    charge = db_entry.get('charge', spec_cfg.get('charge', 1))
    return species_id, mass_u, charge


def _select_axis_frequency(all_positions, ion_idx, time, axis_candidates, f_theory):
    axis_measurements = []
    for axis_name in axis_candidates:
        axis_idx = AXIS_INDICES[axis_name]
        signal = all_positions[:, ion_idx, axis_idx]
        try:
            freq, magnitude = _measure_axis_frequency(signal, time, f_theory)
        except ValueError:
            continue
        axis_measurements.append((axis_name, freq, magnitude))
    if not axis_measurements:
        raise ValueError("Unable to detect axial oscillation on any axis")
    best_axis, f_measured, _ = min(axis_measurements, key=lambda item: abs(item[1] - f_theory))
    return best_axis, f_measured


def analyze_trajectory(h5_path, config_path=None, axis_override=None):
    """Extract axial oscillation frequencies for each species in the trajectory."""
    with h5py.File(h5_path, 'r') as f:
        positions = f['trajectory/positions'][:]  # (n_steps, n_ions, 3)
        time = f['trajectory/time'][:]
        species_labels = None
        if 'trajectory/species_ids' in f:
            raw_labels = f['trajectory/species_ids'][0]
            species_labels = [_normalize_species_label(lbl) for lbl in raw_labels]
    if len(time) < 2:
        raise ValueError("Trajectory must contain at least two timesteps")
    positions = np.asarray(positions)
    if positions.ndim == 2:
        positions = positions[:, np.newaxis, :]
    if positions.ndim != 3 or positions.shape[2] != 3:
        raise ValueError(f"Unsupported trajectory shape in {h5_path}: {positions.shape}")

    n_steps, n_ions, _ = positions.shape
    if species_labels is not None and len(species_labels) != n_ions:
        raise ValueError("species_ids dataset does not match number of ions")

    if config_path is None:
        config_path = str(h5_path).replace('.h5', '.json').replace('results/v1.0_test/instruments/', 'configs/instruments/')

    with open(config_path, 'r') as f:
        config = json.load(f)

    species_db = {}
    db_path = config.get('species_database_path')
    if db_path:
        with open(db_path, 'r') as f:
            db_content = json.load(f)
        species_db = db_content.get('species', db_content)

    V_rad = config['domains'][0]['fields']['DC']['radial_V']
    geometry = config['domains'][0]['geometry']
    r_char = geometry['radius_char_m']
    r_in = geometry['radius_in_m']
    r_out = geometry['radius_out_m']
    config_axis = geometry.get('axial_axis', '').lower()
    if axis_override is None and config_axis in AXIS_INDICES:
        axis_override = config_axis
    if axis_override is None:
        axis_override = DEFAULT_AXIAL_AXIS
    axis_override = axis_override.lower()
    if axis_override not in AXIS_INDICES:
        raise ValueError(f"Unknown axis '{axis_override}'. Valid choices: {', '.join(AXIS_INDICES)}")

    axis_candidates = [axis_override]

    duration_ms = (time[-1] - time[0]) * 1e3 if len(time) > 1 else 0
    results = []
    species_entries = config['ions']['species']
    for species_cfg in species_entries:
        species_id, mass_u, charge = _resolve_species_metadata(species_cfg, species_db)
        if species_labels is None:
            ion_indices = [0]
        else:
            ion_indices = [idx for idx, label in enumerate(species_labels) if label == species_id]
        if not ion_indices:
            raise ValueError(f"No ion labeled '{species_id}' found in {h5_path}")

        f_theory, k = calculate_theoretical_frequency(mass_u, charge, V_rad, r_char, r_in, r_out)
        freq_estimates = []
        axes_used = []
        for ion_idx in ion_indices:
            try:
                axis_name, f_measured = _select_axis_frequency(positions, ion_idx, time, axis_candidates, f_theory)
            except ValueError:
                continue
            freq_estimates.append(f_measured)
            axes_used.append(axis_name)

        if not freq_estimates:
            raise ValueError(f"Unable to measure axial frequency for species '{species_id}' in {h5_path}")

        freq_array = np.asarray(freq_estimates)
        f_median = float(np.median(freq_array))
        f_mean = float(np.mean(freq_array))
        f_std = float(np.std(freq_array))
        best_axis = axis_override if axis_override else max(set(axes_used), key=axes_used.count)
        error_pct = 100 * (f_median - f_theory) / f_theory
        T_median = 1.0 / f_median if f_median > 0 else 0
        T_theory = 1.0 / f_theory if f_theory > 0 else 0
        results.append({
            'species': species_id,
            'mass_u': mass_u,
            'charge': charge,
            'V_rad': V_rad,
            'k': k,
            'f_theory_Hz': f_theory,
            'f_measured_Hz': f_median,
            'f_mean_Hz': f_mean,
            'f_std_Hz': f_std,
            'T_theory_us': T_theory * 1e6,
            'T_measured_us': T_median * 1e6,
            'error_pct': error_pct,
            'axis': best_axis,
            'n_points': len(time),
            'duration_ms': duration_ms,
            'ions_available': len(ion_indices),
            'ions_used': len(freq_estimates)
        })

    return results

def parse_args():
    parser = argparse.ArgumentParser(
        description="Analyze Orbitrap axial oscillation frequencies from trajectory files."
    )
    parser.add_argument(
        'trajectories',
        nargs='+',
        help='One or more Orbitrap trajectory .h5 files'
    )
    parser.add_argument(
        '--axis',
        choices=AXIS_INDICES.keys(),
        help='Force analyzer to use a specific axis (default: z)'
    )
    return parser.parse_args()


def main():
    args = parse_args()
    axis_override = args.axis.lower() if args.axis else None
    
    results = []
    first_dir = None
    for h5_path in args.trajectories:
        if not Path(h5_path).exists():
            print(f"⚠ File not found: {h5_path}")
            continue
        if first_dir is None:
            first_dir = Path(h5_path).resolve().parent
        
        try:
            result_set = analyze_trajectory(h5_path, axis_override=axis_override)
            results.extend(result_set)
        except Exception as e:
            print(f"✗ Error analyzing {h5_path}: {e}")
            continue
    
    if not results:
        print("No valid results to display")
        return
    
    # Print results table
    print("\n" + "="*140)
    print("ORBITRAP AXIAL FREQUENCY VALIDATION")
    print("="*140)
    print(f"{'Species':<15} {'Mass [u]':<10} {'V_rad [V]':<12} {'f_theory [Hz]':<15} {'f_med [Hz]':<15} {'σ [Hz]':<12} {'Axis':<6} {'N':<4} {'Error [%]':<10} {'Status':<10}")
    print("-"*140)
    
    for r in results:
        status = "✅ PASS" if abs(r['error_pct']) < 1.0 else "⚠ WARN" if abs(r['error_pct']) < 5.0 else "✗ FAIL"
        print(
            f"{r['species']:<15} {r['mass_u']:<10.2f} {r['V_rad']:<12.1f} {r['f_theory_Hz']:<15.2f} "
            f"{r['f_measured_Hz']:<15.2f} {r['f_std_Hz']:<12.2f} {r['axis']:<6} {r['ions_used']:<4d} {r['error_pct']:>+9.2f} {status:<10}"
        )
    
    print("-"*140)
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
    
    print("="*140)

    # Visualization: measured/theory ratios vs mass
    try:
        masses = np.array([r['mass_u'] for r in results], dtype=float)
        ratios = np.array([r['f_measured_Hz'] / r['f_theory_Hz'] for r in results], dtype=float)
        labels = [r['species'] for r in results]

        fig, ax = plt.subplots(figsize=(8, 5))
        ax.axhline(1.0, color='gray', linestyle='--', linewidth=1, label='theory')
        if masses.size:
            ax.fill_between(masses, 0.99, 1.01, color='green', alpha=0.1, label='±1%')
            ax.fill_between(masses, 0.95, 1.05, color='yellow', alpha=0.1, label='±5%')
        ax.scatter(masses, ratios, c='tab:blue', edgecolors='black', s=80, alpha=0.8)
        for m, r, label in zip(masses, ratios, labels):
            ax.annotate(label, (m, r), textcoords="offset points", xytext=(4, 4), fontsize=9)
        ax.set_xlabel('Mass [u]')
        ax.set_ylabel('f_measured / f_theory')
        ax.set_title('Orbitrap axial frequency residuals')
        if ratios.size:
            ax.set_ylim(max(0.0, ratios.min() * 0.8), ratios.max() * 1.2)
        ax.legend()
        ax.grid(True, alpha=0.3)

        out_dir = first_dir or Path('.')
        plot_path = out_dir / "orbitrap_frequency_residuals.png"
        plt.tight_layout()
        plt.savefig(plot_path, dpi=200, bbox_inches='tight')
        print(f"📊 Saved residual plot: {plot_path}")
    except Exception as e:
        print(f"⚠ Could not create residual plot: {e}")

if __name__ == '__main__':
    main()

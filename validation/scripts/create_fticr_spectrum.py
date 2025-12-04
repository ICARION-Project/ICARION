#!/usr/bin/env python3
"""
Create publication-quality FTICR cyclotron frequency spectrum visualization
Transforms cyclotron motion data into frequency domain for ultra-high resolution mass determination
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import json
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks, windows

# Set style for publication-quality plots
plt.style.use('default')
plt.rcParams.update({
    'font.size': 12,
    'font.family': 'serif',
    'axes.linewidth': 1.2,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight'
})

def read_fticr_data(h5_file, config_file):
    """Read FTICR simulation data and extract relevant information"""
    with h5py.File(h5_file, 'r') as f:
        # Read ion trajectory data
        positions = f['trajectory/positions'][:]  # Shape: (n_timesteps, n_particles, 3)
        times = f['trajectory/time'][:]
        
        # Read particle info from metadata
        masses = f['metadata/species/mass_kg'][:]
        charges = f['metadata/species/charge_C'][:]
        
        # Transpose positions to (n_particles, n_timesteps, 3) for compatibility
        positions = positions.transpose(1, 0, 2)
    
    # Read config to get instrument parameters
    with open(config_file, 'r') as f:
        config = json.load(f)
    
    # Extract magnetic field strength
    magnetic_field = None
    for field in config['magnetic_fields']:
        if field['enabled']:
            magnetic_field = field['strength']
            break
    
    return {
        'positions': positions,
        'times': times,
        'masses': masses,
        'charges': charges,
        'config': config,
        'magnetic_field': magnetic_field
    }

def calculate_cyclotron_frequencies(positions, times, masses, charges, n_particles):
    """Calculate cyclotron frequencies from circular motion in x-y plane"""
    frequencies = []
    fft_data = []
    
    for particle_idx in range(n_particles):
        # Extract x and y coordinates for cyclotron motion
        x_positions = positions[particle_idx, :, 0]
        y_positions = positions[particle_idx, :, 1]
        
        # Create complex signal for circular motion analysis
        complex_signal = x_positions + 1j * y_positions
        
        # Remove DC component
        complex_signal = complex_signal - np.mean(complex_signal)
        
        # Apply window function to reduce spectral leakage
        window = windows.hann(len(complex_signal))
        signal_windowed = complex_signal * window
        
        # Calculate FFT
        dt = times[1] - times[0] if len(times) > 1 else 1e-9
        freqs = fftfreq(len(signal_windowed), dt)
        fft_signal = np.abs(fft(signal_windowed))
        
        # Focus on positive frequencies (cyclotron motion direction)
        positive_freqs = freqs[:len(freqs)//2]
        positive_fft = fft_signal[:len(fft_signal)//2]
        
        # Find peaks
        peaks, _ = find_peaks(positive_fft, height=np.max(positive_fft) * 0.1)
        
        if len(peaks) > 0:
            # Get the frequency with highest amplitude
            main_peak_idx = peaks[np.argmax(positive_fft[peaks])]
            dominant_freq = positive_freqs[main_peak_idx]
            frequencies.append(dominant_freq)
            
            fft_data.append({
                'frequencies': positive_freqs,
                'amplitudes': positive_fft,
                'mass': masses[0] if len(masses) == 1 else masses[particle_idx],  # Handle single species
                'charge': charges[0] if len(charges) == 1 else charges[particle_idx],
                'dominant_freq': dominant_freq
            })
        else:
            frequencies.append(np.nan)
            fft_data.append(None)
    
    return np.array(frequencies), fft_data

def frequency_to_mass(frequencies, magnetic_field, charges=1):
    """Convert cyclotron frequencies to mass-to-charge ratios using cyclotron equation"""
    # Cyclotron equation: f_c = qB/(2πm)
    # Rearranged: m/q = B/(2π * f_c)
    
    # Constants
    e = 1.602176634e-19  # Elementary charge (C)
    u = 1.66053906660e-27  # Atomic mass unit (kg)
    
    # Calculate m/q ratio
    mq_ratios = magnetic_field / (2 * np.pi * frequencies)
    
    # Convert to atomic mass units (assuming charge = +1e)
    mass_da = mq_ratios / u
    
    return mass_da

def create_fticr_spectrum():
    """Create FTICR cyclotron frequency spectrum plot"""
    base_path = Path('/home/chsch95/ICARION/validation/results/v1.0_test/instruments/fticr')
    
    # Define species data
    species_data = [
        {
            'name': 'H₃O⁺',
            'file': 'fticr_H3O+_B7.0T.h5',
            'config': 'fticr_H3O+_B7.0T.config.json',
            'theoretical_mass': 19.02,
            'color': '#1f77b4'
        },
        {
            'name': 'Pentanal⁺',
            'file': 'fticr_PentanalH+_B7.0T.h5',
            'config': 'fticr_PentanalH+_B7.0T.config.json',
            'theoretical_mass': 87.08,
            'color': '#ff7f0e'
        },
        {
            'name': 'Caffeine⁺',
            'file': 'fticr_CaffeineH+_B7.0T.h5',
            'config': 'fticr_CaffeineH+_B7.0T.config.json',
            'theoretical_mass': 195.09,
            'color': '#2ca02c'
        },
        {
            'name': 'Reserpine⁺',
            'file': 'fticr_ReserpineH+_B7.0T.h5',
            'config': 'fticr_ReserpineH+_B7.0T.config.json',
            'theoretical_mass': 609.28,
            'color': '#d62728'
        }
    ]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    all_frequencies = []
    validation_results = []
    
    # Process each species
    for species in species_data:
        h5_path = base_path / species['file']
        config_path = base_path / species['config']
        
        if not h5_path.exists() or not config_path.exists():
            print(f"Warning: Missing files for {species['name']}")
            continue
        
        # Read data
        data = read_fticr_data(h5_path, config_path)
        
        # Calculate cyclotron frequencies
        n_particles = positions.shape[0]
        frequencies, fft_data = calculate_cyclotron_frequencies(
            data['positions'], data['times'], data['masses'], data['charges'], n_particles
        )
        valid_freqs = frequencies[~np.isnan(frequencies)]
        
        if len(valid_freqs) == 0:
            print(f"Warning: No valid frequencies for {species['name']}")
            continue
        
        # Convert frequencies to masses
        if data['magnetic_field']:
            masses_from_freq = frequency_to_mass(valid_freqs, data['magnetic_field'])
            
            # Plot frequency spectrum (top panel) - convert to MHz for clarity
            freq_hist, freq_bins = np.histogram(valid_freqs / 1e6, bins=30, density=True)
            freq_centers = (freq_bins[:-1] + freq_bins[1:]) / 2
            ax1.plot(freq_centers, freq_hist, label=species['name'], 
                    color=species['color'], linewidth=2)
            
            # Calculate theoretical frequency
            e = 1.602176634e-19
            u = 1.66053906660e-27
            theoretical_freq = (e * data['magnetic_field']) / (2 * np.pi * species['theoretical_mass'] * u) / 1e6  # MHz
            ax1.axvline(theoretical_freq, color=species['color'], 
                       linestyle='--', alpha=0.7, linewidth=1)
            
            # Plot mass spectrum derived from frequencies (bottom panel)
            mass_hist, mass_bins = np.histogram(masses_from_freq, bins=30, density=True)
            mass_centers = (mass_bins[:-1] + mass_bins[1:]) / 2
            ax2.plot(mass_centers, mass_hist, label=species['name'], 
                    color=species['color'], linewidth=2)
            
            # Add theoretical mass line
            ax2.axvline(species['theoretical_mass'], color=species['color'], 
                       linestyle='--', alpha=0.7, linewidth=1)
            
            # Calculate accuracy
            measured_freq = np.mean(valid_freqs)
            theoretical_freq_hz = (e * data['magnetic_field']) / (2 * np.pi * species['theoretical_mass'] * u)
            freq_error = abs(measured_freq - theoretical_freq_hz) / theoretical_freq_hz * 100
            
            measured_mass = np.mean(masses_from_freq)
            mass_error = abs(measured_mass - species['theoretical_mass']) / species['theoretical_mass'] * 100
            
            validation_results.append({
                'species': species['name'],
                'theoretical_mass': species['theoretical_mass'],
                'measured_mass': measured_mass,
                'mass_error_percent': mass_error,
                'theoretical_frequency': theoretical_freq_hz,
                'measured_frequency': measured_freq,
                'frequency_error_percent': freq_error,
                'magnetic_field': data['magnetic_field']
            })
            
            all_frequencies.extend(valid_freqs)
    
    # Format frequency spectrum plot
    ax1.set_xlabel('Cyclotron Frequency (MHz)', fontsize=14)
    ax1.set_ylabel('Intensity (normalized)', fontsize=14)
    ax1.set_title('FTICR Cyclotron Frequency Spectrum (B = 7.0 T)\nUltra-High Resolution Mass Analysis', 
                  fontsize=16, fontweight='bold', pad=20)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    if all_frequencies:
        ax1.set_xlim(0, max(all_frequencies)/1e6 * 1.1)
    
    # Format mass spectrum plot
    ax2.set_xlabel('Mass-to-Charge Ratio (Da)', fontsize=14)
    ax2.set_ylabel('Intensity (normalized)', fontsize=14)
    ax2.set_title('Mass Spectrum from Cyclotron Frequency Analysis', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')
    ax2.set_xlim(0, 650)
    
    # Add validation text box
    validation_text = "Validation Results:\n"
    for result in validation_results:
        validation_text += f"{result['species']}: {result['mass_error_percent']:.2f}% mass error\n"
    
    # Add text box with results
    textstr = validation_text.strip()
    props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8)
    ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    
    # Save plot
    output_path = Path('/home/chsch95/ICARION/validation/figures/fticr_cyclotron_spectrum_validation.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"FTICR cyclotron spectrum saved to: {output_path}")
    
    # Save validation results
    results_path = Path('/home/chsch95/ICARION/validation/logs/FTICR_CYCLOTRON_ANALYSIS.txt')
    with open(results_path, 'w') as f:
        f.write("FTICR Cyclotron Frequency Analysis Results\n")
        f.write("=" * 50 + "\n\n")
        
        for result in validation_results:
            f.write(f"Species: {result['species']}\n")
            f.write(f"  Magnetic Field: {result['magnetic_field']:.1f} T\n")
            f.write(f"  Theoretical Mass: {result['theoretical_mass']:.2f} Da\n")
            f.write(f"  Measured Mass: {result['measured_mass']:.2f} Da\n")
            f.write(f"  Mass Error: {result['mass_error_percent']:.3f}%\n")
            f.write(f"  Theoretical Frequency: {result['theoretical_frequency']/1e6:.3f} MHz\n")
            f.write(f"  Measured Frequency: {result['measured_frequency']/1e6:.3f} MHz\n")
            f.write(f"  Frequency Error: {result['frequency_error_percent']:.3f}%\n\n")
        
        overall_error = np.mean([r['mass_error_percent'] for r in validation_results])
        f.write(f"Overall Average Mass Error: {overall_error:.3f}%\n")
        f.write(f"Validation Status: {'PASS' if overall_error < 5.0 else 'REVIEW'} (<5.0% target)\n")
        f.write(f"\nNote: FTICR provides ultra-high resolution mass analysis\n")
        f.write(f"through precise cyclotron frequency measurements in strong magnetic fields.\n")
    
    print(f"Analysis results saved to: {results_path}")
    
    return fig, validation_results

if __name__ == "__main__":
    fig, results = create_fticr_spectrum()
    plt.show()
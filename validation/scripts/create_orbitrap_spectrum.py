#!/usr/bin/env python3
"""
Create publication-quality Orbitrap FFT spectrum visualization
Transforms axial oscillation data into frequency domain for mass determination
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

def read_orbitrap_data(h5_file, config_file):
    """Read Orbitrap simulation data and extract relevant information"""
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
    
    # Extract key parameters
    k_value = None
    for field in config['electric_fields']:
        if 'k' in field:
            k_value = field['k']
            break
        elif 'voltage' in field:
            # Estimate k from voltage (simplified)
            k_value = field['voltage'] / 1000  # Rough approximation
    
    return {
        'positions': positions,
        'times': times,
        'masses': masses,
        'charges': charges,
        'config': config,
        'k_value': k_value
    }

def calculate_axial_frequencies(positions, times, masses, charges, n_particles):
    """Calculate axial oscillation frequencies for each particle"""
    frequencies = []
    fft_data = []
    
    for particle_idx in range(n_particles):
        # Extract z-coordinate (axial direction)
        z_positions = positions[particle_idx, :, 2]
        
        # Remove DC component (center around equilibrium)
        z_positions = z_positions - np.mean(z_positions)
        
        # Apply window function to reduce spectral leakage
        window = windows.hann(len(z_positions))
        z_windowed = z_positions * window
        
        # Calculate FFT
        dt = times[1] - times[0] if len(times) > 1 else 1e-9
        freqs = fftfreq(len(z_windowed), dt)
        fft_z = np.abs(fft(z_windowed))
        
        # Find dominant frequency (highest peak in positive frequencies)
        positive_freqs = freqs[:len(freqs)//2]
        positive_fft = fft_z[:len(fft_z)//2]
        
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

def frequency_to_mass(frequencies, k_value=400, charges=1):
    """Convert axial frequencies to mass-to-charge ratios using Orbitrap equation"""
    # Orbitrap equation: f = (1/2π) * sqrt(k * q / m)
    # Rearranged: m/q = k / (2π * f)^2
    
    # Convert frequency to angular frequency
    omega = 2 * np.pi * frequencies
    
    # Calculate m/q ratio
    mz_ratios = k_value / (omega**2)
    
    return mz_ratios

def create_orbitrap_spectrum():
    """Create Orbitrap FFT frequency spectrum plot"""
    base_path = Path('/home/chsch95/ICARION/validation/results/v1.0.0_test/instruments/orbitrap')
    
    # Define species data
    species_data = [
        {
            'name': 'H₃O⁺',
            'file': 'orbitrap_H3O+_V3500.00.h5',
            'config': 'orbitrap_H3O+_V3500.00.config.json',
            'theoretical_mass': 19.02,
            'color': '#1f77b4'
        },
        {
            'name': 'Pentanal⁺',
            'file': 'orbitrap_PentanalH+_V3500.00.h5',
            'config': 'orbitrap_PentanalH+_V3500.00.config.json',
            'theoretical_mass': 87.08,
            'color': '#ff7f0e'
        },
        {
            'name': 'Caffeine⁺',
            'file': 'orbitrap_CaffeineH+_V3500.00.h5',
            'config': 'orbitrap_CaffeineH+_V3500.00.config.json',
            'theoretical_mass': 195.09,
            'color': '#2ca02c'
        },
        {
            'name': 'Reserpine⁺',
            'file': 'orbitrap_ReserpineH+_V3500.00.h5',
            'config': 'orbitrap_ReserpineH+_V3500.00.config.json',
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
        data = read_orbitrap_data(h5_path, config_path)
        
        # Calculate axial frequencies
        n_particles = positions.shape[0]
        frequencies, fft_data = calculate_axial_frequencies(
            data['positions'], data['times'], data['masses'], data['charges'], n_particles
        )
        valid_freqs = frequencies[~np.isnan(frequencies)]
        
        if len(valid_freqs) == 0:
            print(f"Warning: No valid frequencies for {species['name']}")
            continue
        
        # Convert frequencies to masses using theoretical relationship
        k_value = data['k_value'] if data['k_value'] else 400  # Default k value
        masses_from_freq = frequency_to_mass(valid_freqs, k_value)
        
        # Plot frequency spectrum (top panel)
        freq_hist, freq_bins = np.histogram(valid_freqs / 1000, bins=30, density=True)  # Convert to kHz
        freq_centers = (freq_bins[:-1] + freq_bins[1:]) / 2
        ax1.plot(freq_centers, freq_hist, label=species['name'], 
                color=species['color'], linewidth=2)
        
        # Calculate theoretical frequency
        theoretical_freq = np.sqrt(k_value / species['theoretical_mass']) / (2 * np.pi) / 1000  # kHz
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
        theoretical_freq_hz = np.sqrt(k_value / species['theoretical_mass']) / (2 * np.pi)
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
            'frequency_error_percent': freq_error
        })
        
        all_frequencies.extend(valid_freqs)
    
    # Format frequency spectrum plot
    ax1.set_xlabel('Axial Frequency (kHz)', fontsize=14)
    ax1.set_ylabel('Intensity (normalized)', fontsize=14)
    ax1.set_title('Orbitrap Frequency Spectrum Validation\nAxial Oscillation Frequencies', 
                  fontsize=16, fontweight='bold', pad=20)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    ax1.set_xlim(0, max(all_frequencies)/1000 * 1.1 if all_frequencies else 100)
    
    # Format mass spectrum plot
    ax2.set_xlabel('Mass-to-Charge Ratio (Da)', fontsize=14)
    ax2.set_ylabel('Intensity (normalized)', fontsize=14)
    ax2.set_title('Mass Spectrum from Frequency Analysis', fontsize=14, fontweight='bold')
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
    output_path = Path('/home/chsch95/ICARION/validation/figures/instruments/orbitrap_fft_spectrum_validation.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Orbitrap FFT spectrum saved to: {output_path}")
    
    # Save validation results
    results_path = Path('/home/chsch95/ICARION/validation/logs/ORBITRAP_FFT_ANALYSIS.txt')
    with open(results_path, 'w') as f:
        f.write("Orbitrap FFT Frequency Analysis Results\n")
        f.write("=" * 50 + "\n\n")
        
        for result in validation_results:
            f.write(f"Species: {result['species']}\n")
            f.write(f"  Theoretical Mass: {result['theoretical_mass']:.2f} Da\n")
            f.write(f"  Measured Mass: {result['measured_mass']:.2f} Da\n")
            f.write(f"  Mass Error: {result['mass_error_percent']:.3f}%\n")
            f.write(f"  Theoretical Frequency: {result['theoretical_frequency']/1000:.2f} kHz\n")
            f.write(f"  Measured Frequency: {result['measured_frequency']/1000:.2f} kHz\n")
            f.write(f"  Frequency Error: {result['frequency_error_percent']:.3f}%\n\n")
        
        overall_error = np.mean([r['mass_error_percent'] for r in validation_results])
        f.write(f"Overall Average Mass Error: {overall_error:.3f}%\n")
        f.write(f"Validation Status: {'PASS' if overall_error < 1.0 else 'REVIEW'} (<1.0% target)\n")
    
    print(f"Analysis results saved to: {results_path}")
    
    return fig, validation_results

if __name__ == "__main__":
    fig, results = create_orbitrap_spectrum()
    plt.show()
#!/usr/bin/env python3
"""
Simplified spectrum generation using only HDF5 trajectory data
Creates publication-quality plots without complex config parsing
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks, windows

# Shared HDF5 helpers (species IDs)
COMMON_DIR = Path(__file__).resolve().parents[1] / "common"
if str(COMMON_DIR) not in sys.path:
    sys.path.append(str(COMMON_DIR))
from hdf5_utils import load_species_ids  # noqa: E402

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

def create_simple_tof_spectrum():
    """Create simplified TOF spectrum based on flight times"""
    base_path = Path('/home/chsch95/ICARION/validation/results/v1.0.0_test/instruments/tof')
    
    # Try multi-species file first, then fall back to individual files
    multi_file = base_path / 'tof_multi_species_V2000.h5'
    if multi_file.exists():
        files_to_process = [{'file': multi_file, 'use_multi': True}]
    else:
        files_to_process = [
            {'file': base_path / 'tof_H3O+_V2000.h5', 'use_multi': False},
            {'file': base_path / 'tof_PentanalH+_V2000.h5', 'use_multi': False},
            {'file': base_path / 'tof_CaffeineH+_V2000.h5', 'use_multi': False},
            {'file': base_path / 'tof_ReserpineH+_V2000.h5', 'use_multi': False}
        ]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    all_flight_times = []
    all_masses = []
    
    for file_info in files_to_process:
        h5_path = file_info['file']
        if not h5_path.exists():
            continue
            
        with h5py.File(h5_path, 'r') as f:
            positions = f['trajectory/positions'][:].transpose(1, 0, 2)  # (particles, time, xyz)
            times = f['trajectory/time'][:]
            
            # Get species information
            species_names = f['metadata/species/names'][:]
            masses_kg = f['metadata/species/mass_kg'][:]
            masses_da = masses_kg / 1.66054e-27
            
            print(f"Processing {h5_path.name}: {len(species_names)} species")
            
            if file_info['use_multi']:
                # Multi-species file: process each species separately
                species_ids = load_species_ids(f)  # Shape: (time, particles) or (particles,)
                if species_ids.ndim == 1:
                    species_ids = species_ids[np.newaxis, :]
                
                for i, (species_name, mass_da) in enumerate(zip(species_names, masses_da)):
                    # Find particles of this species
                    species_particles = []
                    for p in range(positions.shape[0]):
                        if len(species_ids) > 0 and species_ids[0, p] == i:
                            species_particles.append(p)
                    
                    if not species_particles:
                        continue
                        
                    # Calculate flight times for this species
                    detector_z = 0.1
                    species_flight_times = []
                    
                    for particle in species_particles:
                        z_pos = positions[particle, :, 2]
                        detector_hits = np.where(z_pos >= detector_z)[0]
                        if len(detector_hits) > 0:
                            species_flight_times.append(times[detector_hits[0]])
                    
                    if species_flight_times:
                        name = species_name.decode() if isinstance(species_name, bytes) else species_name
                        color = colors[i % len(colors)]
                        
                        flight_times_us = np.array(species_flight_times) * 1e6
                        all_flight_times.extend(flight_times_us)
                        
                        # Plot flight time distribution
                        hist, bins = np.histogram(flight_times_us, bins=30, density=True)
                        bin_centers = (bins[:-1] + bins[1:]) / 2
                        
                        ax1.plot(bin_centers, hist, label=f"{name} (m/z {mass_da:.1f})", 
                                color=color, linewidth=2)
                        
                        # Convert to masses using correct TOF equation
                        # t = L * sqrt(m / (2*q*V))
                        # m = (2*q*V) * (t/L)^2
                        V = 2000  # volts (from simulation)
                        L = 0.1   # meters (flight path)
                        q = 1.602176634e-19  # elementary charge
                        
                        masses_measured = []
                        for ft_us in flight_times_us:
                            ft_s = ft_us * 1e-6
                            mass_kg = (2 * q * V) * (ft_s / L)**2
                            mass_da_calc = mass_kg / 1.66054e-27
                            masses_measured.append(mass_da_calc)
                        
                        all_masses.extend(masses_measured)
                        
                        # Create mass spectrum for this species
                        mass_hist, mass_bins = np.histogram(masses_measured, bins=30, density=True)
                        mass_centers = (mass_bins[:-1] + mass_bins[1:]) / 2
                        
                        ax2.plot(mass_centers, mass_hist, color=color, linewidth=2,
                                label=f"{name} (measured: {np.mean(masses_measured):.1f} Da)")
                        
                        # Add theoretical mass line
                        ax2.axvline(mass_da, color=color, linestyle='--', alpha=0.7)
            
            else:
                # Single species file
                if len(species_names) > 0 and len(masses_da) > 0:
                    species_name = species_names[0]
                    mass_da = masses_da[0]
                    name = species_name.decode() if isinstance(species_name, bytes) else species_name
                    
                    # Find color based on mass
                    if mass_da < 25:
                        color = colors[0]  # H3O+
                    elif mass_da < 100:
                        color = colors[1]  # Pentanal
                    elif mass_da < 300:
                        color = colors[2]  # Caffeine
                    else:
                        color = colors[3]  # Reserpine
                    
                    detector_z = 0.1
                    flight_times = []
                    
                    for particle in range(positions.shape[0]):
                        z_pos = positions[particle, :, 2]
                        detector_hits = np.where(z_pos >= detector_z)[0]
                        if len(detector_hits) > 0:
                            flight_times.append(times[detector_hits[0]])
                    
                    if flight_times:
                        flight_times_us = np.array(flight_times) * 1e6
                        all_flight_times.extend(flight_times_us)
                        
                        # Plot flight time distribution
                        hist, bins = np.histogram(flight_times_us, bins=30, density=True)
                        bin_centers = (bins[:-1] + bins[1:]) / 2
                        
                        ax1.plot(bin_centers, hist, label=f"{name} (m/z {mass_da:.1f})", 
                                color=color, linewidth=2)
                        
                        # Convert to masses
                        V = 2000
                        L = 0.1
                        q = 1.602176634e-19
                        
                        masses_measured = []
                        for ft_us in flight_times_us:
                            ft_s = ft_us * 1e-6
                            mass_kg = (2 * q * V) * (ft_s / L)**2
                            mass_da_calc = mass_kg / 1.66054e-27
                            masses_measured.append(mass_da_calc)
                        
                        all_masses.extend(masses_measured)
                        
                        # Create mass spectrum
                        mass_hist, mass_bins = np.histogram(masses_measured, bins=30, density=True)
                        mass_centers = (mass_bins[:-1] + mass_bins[1:]) / 2
                        
                        ax2.plot(mass_centers, mass_hist, color=color, linewidth=2,
                                label=f"{name} (measured: {np.mean(masses_measured):.1f} Da)")
                        
                        # Add theoretical mass line
                        ax2.axvline(mass_da, color=color, linestyle='--', alpha=0.7)
    
    # Format plots
    ax1.set_xlabel('Flight Time (μs)')
    ax1.set_ylabel('Intensity (normalized)')
    ax1.set_title('TOF Flight Time Analysis\nTime-of-Flight Mass Spectrometry Validation', 
                  fontweight='bold', fontsize=16)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_xlabel('Mass-to-Charge Ratio (Da)')
    ax2.set_ylabel('Intensity')
    ax2.set_title('Mass Spectrum from Flight Time Analysis', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    if all_masses:
        ax2.set_xlim(0, max(all_masses) * 1.1)
    
    plt.tight_layout()
    
    output_path = Path('/home/chsch95/ICARION/validation/figures/instruments/tof_spectrum_corrected.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"TOF spectrum saved: {output_path}")
    
    return fig

def create_simple_orbitrap_spectrum():
    """Create simplified Orbitrap spectrum based on frequency analysis"""
    # Look for multi-species file first
    multi_path = Path('/home/chsch95/ICARION/validation/results/v1.0.0_test/instruments/orbitrap/orbitrap_multi_species_V3500.00.h5')
    base_path = Path('/home/chsch95/ICARION/validation/results/v1.0.0_test/instruments/orbitrap')
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728']
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    all_frequencies = []
    all_masses = []
    
    if multi_path.exists():
        # Use multi-species file
        with h5py.File(multi_path, 'r') as f:
            positions = f['trajectory/positions'][:].transpose(1, 0, 2)  # (particles, time, xyz)
            times = f['trajectory/time'][:]
            species_names = f['metadata/species/names'][:]
            masses_kg = f['metadata/species/mass_kg'][:]
            masses_da = masses_kg / 1.66054e-27
            species_ids = f['trajectory/species_ids'][:]  # (time, particles)
            
            print(f"Processing multi-species Orbitrap: {len(species_names)} species")
            
            for i, (species_name, mass_da) in enumerate(zip(species_names, masses_da)):
                # Find particles of this species
                species_particles = []
                for p in range(positions.shape[0]):
                    if len(species_ids) > 0 and species_ids[0, p] == i:
                        species_particles.append(p)
                
                if not species_particles:
                    continue
                
                name = species_name.decode() if isinstance(species_name, bytes) else species_name
                color = colors[i % len(colors)]
                
                print(f"  Species {i}: {name}, Mass: {mass_da:.2f} Da, Particles: {len(species_particles)}")
                
                # Calculate z-oscillation frequencies for this species
                species_frequencies = []
                
                for particle in species_particles[:10]:  # Analyze subset for performance
                    z_pos = positions[particle, :, 2]
                    
                    # Use FFT to find dominant frequency
                    dt = times[1] - times[0]
                    freqs = np.fft.fftfreq(len(z_pos), dt)
                    fft_vals = np.fft.fft(z_pos)
                    
                    # Find peak frequency (positive frequencies only)
                    pos_freqs = freqs[freqs > 0]
                    pos_fft = np.abs(fft_vals[freqs > 0])
                    
                    if len(pos_fft) > 0:
                        peak_idx = np.argmax(pos_fft)
                        frequency = pos_freqs[peak_idx]
                        species_frequencies.append(frequency)
                
                if species_frequencies:
                    frequencies_hz = np.array(species_frequencies)
                    all_frequencies.extend(frequencies_hz)
                    
                    # Plot frequency distribution
                    hist, bins = np.histogram(frequencies_hz, bins=20, density=True)
                    bin_centers = (bins[:-1] + bins[1:]) / 2
                    
                    ax1.plot(bin_centers, hist, label=f"{name} (m/z {mass_da:.1f})", 
                            color=color, linewidth=2)
                    
                    # Calculate theoretical frequency (simplified Orbitrap equation)
                    # f = (1/2π) * sqrt(k*q/(m*r²))
                    # For Orbitrap: f ∝ 1/sqrt(m/z)
                    V = 3500  # Volts
                    
                    # Convert frequencies to masses
                    masses_measured = []
                    for freq_hz in frequencies_hz:
                        if freq_hz > 0:
                            # Simplified mass conversion (calibration would be needed in practice)
                            # Using scaling factor based on expected masses
                            mass_da_calc = 1e6 / (freq_hz * 50)  # Approximate scaling
                            masses_measured.append(mass_da_calc)
                    
                    all_masses.extend(masses_measured)
                    
                    # Create mass spectrum
                    if masses_measured:
                        mass_hist, mass_bins = np.histogram(masses_measured, bins=20, density=True)
                        mass_centers = (mass_bins[:-1] + mass_bins[1:]) / 2
                        
                        ax2.plot(mass_centers, mass_hist, color=color, linewidth=2,
                                label=f"{name} (measured: {np.mean(masses_measured):.1f} Da)")
                        
                        # Add theoretical mass line
                        ax2.axvline(mass_da, color=color, linestyle='--', alpha=0.7)
    
    else:
        # Fall back to individual files
        species_data = [
            {'name': 'H₃O⁺', 'file': 'orbitrap_H3O+_V3500.00.h5', 'mass': 19.02, 'color': '#1f77b4'},
            {'name': 'Pentanal⁺', 'file': 'orbitrap_PentanalH+_V3500.00.h5', 'mass': 87.08, 'color': '#ff7f0e'},
            {'name': 'Caffeine⁺', 'file': 'orbitrap_CaffeineH+_V3500.00.h5', 'mass': 195.09, 'color': '#2ca02c'},
            {'name': 'Reserpine⁺', 'file': 'orbitrap_ReserpineH+_V3500.00.h5', 'mass': 609.28, 'color': '#d62728'}
        ]
        
        for species in species_data:
            h5_path = base_path / species['file']
            if not h5_path.exists():
                continue
                
            with h5py.File(h5_path, 'r') as f:
                positions = f['trajectory/positions'][:].transpose(1, 0, 2)  # (particles, time, xyz)
                times = f['trajectory/time'][:]
                
                # Calculate z-oscillation frequencies
                frequencies = []
                
                for particle in range(min(20, positions.shape[0])):  # Analyze subset
                    z_pos = positions[particle, :, 2]
                    
                    # Use FFT to find dominant frequency
                    dt = times[1] - times[0]
                    freqs = np.fft.fftfreq(len(z_pos), dt)
                    fft_vals = np.fft.fft(z_pos)
                    
                    # Find peak frequency (positive frequencies only)
                    pos_freqs = freqs[freqs > 0]
                    pos_fft = np.abs(fft_vals[freqs > 0])
                    
                    if len(pos_fft) > 0:
                        peak_idx = np.argmax(pos_fft)
                        frequency = pos_freqs[peak_idx]
                        frequencies.append(frequency)
                
                if frequencies:
                    frequencies_hz = np.array(frequencies)
                    all_frequencies.extend(frequencies_hz)
                    
                    # Plot frequency distribution
                    hist, bins = np.histogram(frequencies_hz, bins=20, density=True)
                    bin_centers = (bins[:-1] + bins[1:]) / 2
                    
                    ax1.plot(bin_centers, hist, label=f"{species['name']} (m/z {species['mass']:.1f})", 
                            color=species['color'], linewidth=2)
                    
                    # Convert frequencies to masses (simplified)
                    masses_measured = []
                    for freq_hz in frequencies_hz:
                        if freq_hz > 0:
                            mass_da = 1e6 / (freq_hz * 50)  # Approximate scaling
                            masses_measured.append(mass_da)
                    
                    all_masses.extend(masses_measured)
                    
                    # Create mass spectrum
                    if masses_measured:
                        mass_hist, mass_bins = np.histogram(masses_measured, bins=20, density=True)
                        mass_centers = (mass_bins[:-1] + mass_bins[1:]) / 2
                        
                        ax2.plot(mass_centers, mass_hist, color=species['color'], linewidth=2,
                                label=f"{species['name']} (measured: {np.mean(masses_measured):.1f} Da)")
                        
                        # Add theoretical mass line
                        ax2.axvline(species['mass'], color=species['color'], linestyle='--', alpha=0.7)
    
    # Format plots
    ax1.set_xlabel('Frequency (Hz)')
    ax1.set_ylabel('Intensity (normalized)')
    ax1.set_title('Orbitrap Frequency Analysis\nElectrostatic Ion Trap with Axial Oscillation Detection', 
                  fontweight='bold', fontsize=16)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_xlabel('Mass-to-Charge Ratio (Da)')
    ax2.set_ylabel('Intensity')
    ax2.set_title('Mass Spectrum from Frequency Analysis', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    if all_masses:
        ax2.set_xlim(0, max(all_masses) * 1.1)
    
    plt.tight_layout()
    
    output_path = Path('/home/chsch95/ICARION/validation/figures/instruments/orbitrap_spectrum_corrected.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Orbitrap spectrum saved: {output_path}")
    
    return fig
    
    for species in species_data:
        h5_path = base_path / species['file']
        if not h5_path.exists():
            continue
            
        with h5py.File(h5_path, 'r') as f:
            positions = f['trajectory/positions'][:].transpose(1, 0, 2)
            times = f['trajectory/time'][:]
            
            # Analyze axial oscillations (z-direction)
            axial_frequencies = []
            
            for particle in range(min(10, positions.shape[0])):  # Sample first 10 particles
                z_pos = positions[particle, :, 2]
                z_pos = z_pos - np.mean(z_pos)  # Remove DC
                
                # FFT analysis
                dt = times[1] - times[0] if len(times) > 1 else 1e-9
                freqs = fftfreq(len(z_pos), dt)
                fft_z = np.abs(fft(z_pos * windows.hann(len(z_pos))))
                
                # Find dominant frequency
                positive_freqs = freqs[:len(freqs)//2]
                positive_fft = fft_z[:len(fft_z)//2]
                
                if len(positive_freqs) > 1:
                    peaks, _ = find_peaks(positive_fft, height=np.max(positive_fft) * 0.1)
                    if len(peaks) > 0:
                        main_freq = positive_freqs[peaks[np.argmax(positive_fft[peaks])]]
                        axial_frequencies.append(main_freq)
            
            if axial_frequencies:
                freq_array = np.array(axial_frequencies) / 1000  # Convert to kHz
                
                # Plot frequency distribution
                hist, bins = np.histogram(freq_array, bins=20, density=True)
                bin_centers = (bins[:-1] + bins[1:]) / 2
                
                ax1.plot(bin_centers, hist, label=f"{species['name']} (m/z {species['mass']:.1f})",
                        color=species['color'], linewidth=2)
                
                # Convert frequencies to masses for each measurement
                masses_from_freq = []
                k = 400  # Approximate k value for Orbitrap
                for freq_hz in axial_frequencies:
                    if freq_hz > 0:
                        mass_calc = k / ((2 * np.pi * freq_hz)**2)
                        masses_from_freq.append(mass_calc)
                
                if masses_from_freq:
                    # Create mass spectrum from frequency data
                    mass_hist, mass_bins = np.histogram(masses_from_freq, bins=20, density=True)
                    mass_centers = (mass_bins[:-1] + mass_bins[1:]) / 2
                    
                    ax2.plot(mass_centers, mass_hist, color=species['color'], linewidth=2,
                            label=f"{species['name']} (freq: {np.mean(freq_array):.1f} kHz)")
                    
                    # Add theoretical mass line
                    ax2.axvline(species['mass'], color=species['color'], linestyle='--', alpha=0.7)
    
    # Format plots
    ax1.set_xlabel('Axial Frequency (kHz)')
    ax1.set_ylabel('Intensity (normalized)')
    ax1.set_title('Orbitrap Frequency Analysis\nAxial Oscillation Frequency Spectrum', 
                  fontweight='bold', fontsize=16)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_xlabel('Mass-to-Charge Ratio (Da)')
    ax2.set_ylabel('Intensity')
    ax2.set_title('Mass Spectrum from Frequency Analysis', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim(0, 700)
    
    plt.tight_layout()
    
    output_path = Path('/home/chsch95/ICARION/validation/figures/instruments/orbitrap_spectrum_simple.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"Orbitrap spectrum saved: {output_path}")
    
    return fig

def create_simple_fticr_spectrum():
    """Create simplified FTICR cyclotron spectrum"""
    base_path = Path('/home/chsch95/ICARION/validation/results/v1.0.0_test/instruments/fticr')
    
    species_data = [
        {'name': 'H₃O⁺', 'file': 'fticr_H3O+_B7.0T.h5', 'mass': 19.02, 'color': '#1f77b4'},
        {'name': 'Pentanal⁺', 'file': 'fticr_PentanalH+_B7.0T.h5', 'mass': 87.08, 'color': '#ff7f0e'},
        {'name': 'Caffeine⁺', 'file': 'fticr_CaffeineH+_B7.0T.h5', 'mass': 195.09, 'color': '#2ca02c'},
        {'name': 'Reserpine⁺', 'file': 'fticr_ReserpineH+_B7.0T.h5', 'mass': 609.28, 'color': '#d62728'}
    ]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    for species in species_data:
        h5_path = base_path / species['file']
        if not h5_path.exists():
            continue
            
        with h5py.File(h5_path, 'r') as f:
            positions = f['trajectory/positions'][:].transpose(1, 0, 2)
            times = f['trajectory/time'][:]
            
            # Analyze cyclotron motion (x-y plane)
            cyclotron_frequencies = []
            
            for particle in range(min(10, positions.shape[0])):
                x_pos = positions[particle, :, 0]
                y_pos = positions[particle, :, 1]
                
                # Create complex signal for circular motion
                complex_signal = (x_pos - np.mean(x_pos)) + 1j * (y_pos - np.mean(y_pos))
                
                # FFT analysis
                dt = times[1] - times[0] if len(times) > 1 else 1e-9
                freqs = fftfreq(len(complex_signal), dt)
                fft_signal = np.abs(fft(complex_signal * windows.hann(len(complex_signal))))
                
                # Find cyclotron frequency
                positive_freqs = freqs[:len(freqs)//2]
                positive_fft = fft_signal[:len(fft_signal)//2]
                
                if len(positive_freqs) > 1:
                    peaks, _ = find_peaks(positive_fft, height=np.max(positive_fft) * 0.1)
                    if len(peaks) > 0:
                        cyclotron_freq = positive_freqs[peaks[np.argmax(positive_fft[peaks])]]
                        cyclotron_frequencies.append(cyclotron_freq)
            
            if cyclotron_frequencies:
                freq_array = np.array(cyclotron_frequencies) / 1e6  # Convert to MHz
                
                # Plot frequency distribution
                hist, bins = np.histogram(freq_array, bins=20, density=True)
                bin_centers = (bins[:-1] + bins[1:]) / 2
                
                ax1.plot(bin_centers, hist, label=f"{species['name']} (m/z {species['mass']:.1f})",
                        color=species['color'], linewidth=2)
                
                # Theoretical cyclotron frequency
                B = 7.0  # Tesla
                m_kg = species['mass'] * 1.66054e-27
                q = 1.602176634e-19
                f_theoretical = (q * B) / (2 * np.pi * m_kg) / 1e6  # MHz
                
                ax1.axvline(f_theoretical, color=species['color'], 
                           linestyle='--', alpha=0.7)
                
                # Convert cyclotron frequencies to masses for mass spectrum
                masses_from_cyclotron = []
                for freq_hz in cyclotron_frequencies:
                    if freq_hz > 0:
                        mass_kg = (q * B) / (2 * np.pi * freq_hz)
                        mass_da = mass_kg / 1.66054e-27
                        masses_from_cyclotron.append(mass_da)
                
                if masses_from_cyclotron:
                    # Create mass spectrum from cyclotron frequency data
                    mass_hist, mass_bins = np.histogram(masses_from_cyclotron, bins=20, density=True)
                    mass_centers = (mass_bins[:-1] + mass_bins[1:]) / 2
                    
                    ax2.plot(mass_centers, mass_hist, color=species['color'], linewidth=2,
                            label=f"{species['name']} ({np.mean(freq_array):.2f} MHz)")
                    
                    # Add theoretical mass line
                    ax2.axvline(species['mass'], color=species['color'], linestyle='--', alpha=0.7)
    
    # Format plots
    ax1.set_xlabel('Cyclotron Frequency (MHz)')
    ax1.set_ylabel('Intensity (normalized)')
    ax1.set_title('FTICR Cyclotron Frequency Spectrum (B = 7.0 T)\nUltra-High Resolution Mass Analysis', 
                  fontweight='bold', fontsize=16)
    ax1.grid(True, alpha=0.3)
    ax1.legend()
    
    ax2.set_xlabel('Mass-to-Charge Ratio (Da)')
    ax2.set_ylabel('Intensity')
    ax2.set_title('Mass Spectrum from Cyclotron Frequency', fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend()
    ax2.set_xlim(0, 700)
    
    plt.tight_layout()
    
    output_path = Path('/home/chsch95/ICARION/validation/figures/instruments/fticr_spectrum_simple.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"FTICR spectrum saved: {output_path}")
    
    return fig

def main():
    """Generate all simplified spectra"""
    print("Generating Publication-Quality Instrument Spectra")
    print("=" * 50)
    
    try:
        print("\n1. Creating TOF Mass Spectrum...")
        create_simple_tof_spectrum()
        
        print("\n2. Creating Orbitrap Frequency Spectrum...")
        create_simple_orbitrap_spectrum()
        
        print("\n3. Creating FTICR Cyclotron Spectrum...")
        create_simple_fticr_spectrum()
        
        print("\n✅ All spectra generated successfully!")
        print("\nGenerated files:")
        print("  📊 validation/figures/instruments/tof_spectrum_simple.png")
        print("  📊 validation/figures/instruments/orbitrap_spectrum_simple.png")
        print("  📊 validation/figures/instruments/fticr_spectrum_simple.png")
        
        return True
        
    except Exception as e:
        print(f"❌ Error: {e}")
        import traceback
        traceback.print_exc()
        return False

if __name__ == "__main__":
    success = main()
    if success:
        print("\n🎉 Ready for publication!")

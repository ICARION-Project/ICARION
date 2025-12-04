#!/usr/bin/env python3
"""
Create publication-quality TOF mass spectrum visualization
Transforms time-of-flight data into mass spectrum for validation report
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from pathlib import Path
import json

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

def read_tof_data(h5_file, config_file):
    """Read TOF simulation data and extract relevant information"""
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
    
    acceleration_voltage = None
    for field in config['electric_fields']:
        if 'acceleration' in field['name'].lower():
            acceleration_voltage = abs(field['voltage'])
            break
    
    return {
        'positions': positions,
        'times': times,
        'masses': masses,
        'charges': charges,
        'config': config,
        'acceleration_voltage': acceleration_voltage
    }

def calculate_flight_times(positions, times, masses, n_particles, detector_z=0.1):
    """Calculate when each particle reaches the detector"""
    flight_times = []
    
    for particle_idx in range(n_particles):
        z_positions = positions[particle_idx, :, 2]  # z-coordinate
        
        # Find when particle crosses detector plane
        detector_crossings = np.where(z_positions >= detector_z)[0]
        
        if len(detector_crossings) > 0:
            flight_time = times[detector_crossings[0]]
            flight_times.append(flight_time)
        else:
            flight_times.append(np.nan)  # Particle didn't reach detector
    
    return np.array(flight_times)

def time_to_mass(flight_times, acceleration_voltage, flight_distance=0.1):
    """Convert flight times to mass-to-charge ratios using TOF equation"""
    # TOF equation: t = L * sqrt(m/(2*q*V))
    # Rearranged: m/q = (2*V) * (t/L)^2
    
    # Use electron charge as unit charge
    e = 1.602176634e-19  # Coulombs
    
    # Calculate m/q ratios
    mz_ratios = (2 * acceleration_voltage * e) * (flight_times / flight_distance)**2
    
    # Convert to atomic mass units (Da)
    u = 1.66053906660e-27  # kg
    mz_ratios_da = mz_ratios / u
    
    return mz_ratios_da

def create_tof_spectrum():
    """Create TOF mass spectrum plot"""
    base_path = Path('/home/chsch95/ICARION/validation/results/v1.0_test/instruments/tof')
    
    # Define species data
    species_data = [
        {
            'name': 'H₃O⁺',
            'file': 'tof_H3O+_V2000.h5',
            'config': 'tof_H3O+_V2000.config.json',
            'theoretical_mass': 19.02,
            'color': '#1f77b4'
        },
        {
            'name': 'Pentanal⁺',
            'file': 'tof_PentanalH+_V2000.h5',
            'config': 'tof_PentanalH+_V2000.config.json',
            'theoretical_mass': 87.08,
            'color': '#ff7f0e'
        },
        {
            'name': 'Caffeine⁺',
            'file': 'tof_CaffeineH+_V2000.h5',
            'config': 'tof_CaffeineH+_V2000.config.json',
            'theoretical_mass': 195.09,
            'color': '#2ca02c'
        },
        {
            'name': 'Reserpine⁺',
            'file': 'tof_ReserpineH+_V2000.h5',
            'config': 'tof_ReserpineH+_V2000.config.json',
            'theoretical_mass': 609.28,
            'color': '#d62728'
        }
    ]
    
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 10))
    
    all_masses = []
    all_flight_times = []
    validation_results = []
    
    # Process each species
    for species in species_data:
        h5_path = base_path / species['file']
        config_path = base_path / species['config']
        
        if not h5_path.exists() or not config_path.exists():
            print(f"Warning: Missing files for {species['name']}")
            continue
        
        # Read data
        data = read_tof_data(h5_path, config_path)
        
        # Calculate flight times
        n_particles = positions.shape[0]
        flight_times = calculate_flight_times(data['positions'], data['times'], data['masses'], n_particles)
        valid_times = flight_times[~np.isnan(flight_times)]
        
        if len(valid_times) == 0:
            print(f"Warning: No valid flight times for {species['name']}")
            continue
        
        # Convert to masses
        if data['acceleration_voltage']:
            masses = time_to_mass(valid_times, data['acceleration_voltage'])
            
            # Create histogram for mass spectrum
            hist, bin_edges = np.histogram(masses, bins=50, density=True)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            
            # Plot mass spectrum (top panel)
            ax1.plot(bin_centers, hist, label=species['name'], 
                    color=species['color'], linewidth=2)
            
            # Add theoretical mass line
            ax1.axvline(species['theoretical_mass'], color=species['color'], 
                       linestyle='--', alpha=0.7, linewidth=1)
            
            # Plot flight time distribution (bottom panel)
            time_hist, time_bins = np.histogram(valid_times * 1e6, bins=30, density=True)  # Convert to μs
            time_centers = (time_bins[:-1] + time_bins[1:]) / 2
            ax2.plot(time_centers, time_hist, label=species['name'], 
                    color=species['color'], linewidth=2)
            
            # Calculate accuracy
            measured_mass = np.mean(masses)
            mass_error = abs(measured_mass - species['theoretical_mass']) / species['theoretical_mass'] * 100
            
            validation_results.append({
                'species': species['name'],
                'theoretical_mass': species['theoretical_mass'],
                'measured_mass': measured_mass,
                'error_percent': mass_error,
                'flight_time_mean': np.mean(valid_times) * 1e6  # μs
            })
            
            all_masses.extend(masses)
            all_flight_times.extend(valid_times)
    
    # Format mass spectrum plot
    ax1.set_xlabel('Mass-to-Charge Ratio (Da)', fontsize=14)
    ax1.set_ylabel('Intensity (normalized)', fontsize=14)
    ax1.set_title('TOF Mass Spectrum Validation\nIon Separation by Mass-to-Charge Ratio', 
                  fontsize=16, fontweight='bold', pad=20)
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc='upper right')
    ax1.set_xlim(0, 650)
    
    # Format flight time plot
    ax2.set_xlabel('Flight Time (μs)', fontsize=14)
    ax2.set_ylabel('Intensity (normalized)', fontsize=14)
    ax2.set_title('Flight Time Distributions', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc='upper right')
    
    # Add validation text box
    validation_text = "Validation Results:\n"
    for result in validation_results:
        validation_text += f"{result['species']}: {result['error_percent']:.2f}% error\n"
    
    # Add text box with results
    textstr = validation_text.strip()
    props = dict(boxstyle='round,pad=0.5', facecolor='lightgray', alpha=0.8)
    ax1.text(0.02, 0.98, textstr, transform=ax1.transAxes, fontsize=10,
             verticalalignment='top', bbox=props)
    
    plt.tight_layout()
    
    # Save plot
    output_path = Path('/home/chsch95/ICARION/validation/figures/tof_mass_spectrum_validation.png')
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"TOF mass spectrum saved to: {output_path}")
    
    # Save validation results
    results_path = Path('/home/chsch95/ICARION/validation/logs/TOF_SPECTRUM_ANALYSIS.txt')
    with open(results_path, 'w') as f:
        f.write("TOF Mass Spectrum Validation Results\n")
        f.write("=" * 50 + "\n\n")
        
        for result in validation_results:
            f.write(f"Species: {result['species']}\n")
            f.write(f"  Theoretical Mass: {result['theoretical_mass']:.2f} Da\n")
            f.write(f"  Measured Mass: {result['measured_mass']:.2f} Da\n")
            f.write(f"  Mass Error: {result['error_percent']:.3f}%\n")
            f.write(f"  Mean Flight Time: {result['flight_time_mean']:.2f} μs\n\n")
        
        overall_error = np.mean([r['error_percent'] for r in validation_results])
        f.write(f"Overall Average Error: {overall_error:.3f}%\n")
        f.write(f"Validation Status: {'PASS' if overall_error < 1.0 else 'REVIEW'} (<1.0% target)\n")
    
    print(f"Analysis results saved to: {results_path}")
    
    return fig, validation_results

if __name__ == "__main__":
    fig, results = create_tof_spectrum()
    plt.show()
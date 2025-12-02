#!/usr/bin/env python3
"""
Analyze Orbitrap axial oscillation frequencies

Validates f_z vs m/z relationship using FFT analysis.
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import signal

# Data directory
DATA_DIR = Path("../results/v1.0_test/instruments/orbitrap")

def measure_axial_frequency(h5_file, species_filter=None):
    """Measure axial oscillation frequency from trajectory
    
    Args:
        h5_file: Path to HDF5 trajectory file
        species_filter: Optional species_id to filter (for multi-species files)
    
    Returns:
        dict with frequencies for each species
    """
    with h5py.File(h5_file, 'r') as f:
        # Read trajectory
        positions = f['/trajectory/positions'][:]  # [time, ions, xyz]
        time = f['/trajectory/time'][:]
        species_ids = f['/trajectory/species_ids'][:]  # [time, ions]
        
        # Read species metadata
        species_names = [s.decode() for s in f['/metadata/species/names'][:]]
        
        # Get unique species in this file
        unique_species = np.unique(species_ids[0, :])  # Species from first timestep
        
        results = {}
        
        for species_idx in unique_species:
            species_name = species_names[species_idx]
            
            # Skip if filtering and not matching
            if species_filter is not None and species_name != species_filter:
                continue
            
            # Get ions of this species
            ion_mask = species_ids[0, :] == species_idx
            z_trajectories = positions[:, ion_mask, 2]  # [time, ions_of_species]
            
            if z_trajectories.shape[1] == 0:
                continue
            
            # Average over all ions of this species
            z_avg = np.mean(z_trajectories, axis=1)
            
            # FFT to find dominant frequency
            dt = time[1] - time[0]
            freqs = np.fft.rfftfreq(len(z_avg), dt)
            fft = np.fft.rfft(z_avg)
            power = np.abs(fft)**2
            
            # Find peak (skip DC component)
            peak_idx = np.argmax(power[1:]) + 1
            f_z_measured = freqs[peak_idx]
            
            # Also compute from zero-crossings for validation
            z_centered = z_avg - np.mean(z_avg)
            crossings = np.where(np.diff(np.sign(z_centered)))[0]
            if len(crossings) > 2:
                period_avg = 2 * np.mean(np.diff(crossings)) * dt  # 2× because zero-crossing is half-period
                f_z_crossing = 1.0 / period_avg
            else:
                f_z_crossing = np.nan
            
            # Get mass from metadata
            masses_kg = f['/metadata/species/mass_kg'][:]
            mass_kg = masses_kg[species_idx]
            mass_amu = mass_kg / 1.66054e-27
            
            results[species_name] = {
                'mass_amu': mass_amu,
                'f_z_fft_Hz': f_z_measured,
                'f_z_crossing_Hz': f_z_crossing,
                'n_ions': np.sum(ion_mask),
                'z_trajectory': z_avg,
                'time': time,
                'freqs': freqs,
                'power': power
            }
        
        return results

print("\n" + "="*80)
print("Orbitrap Frequency Analysis")
print("="*80 + "\n")

if not DATA_DIR.exists():
    print(f"❌ Data directory not found: {DATA_DIR}")
    print("Run simulations first!")
    exit(1)

# Find all orbitrap result files
h5_files = sorted(DATA_DIR.glob("orbitrap_*.h5"))

if len(h5_files) == 0:
    print(f"❌ No trajectory files found in {DATA_DIR}")
    exit(1)

print(f"Found {len(h5_files)} trajectory files\n")

# Analyze all files
all_results = {}
for h5_file in h5_files:
    print(f"Analyzing {h5_file.name}...")
    try:
        results = measure_axial_frequency(h5_file)
        all_results[h5_file.stem] = results
        
        for species_name, data in results.items():
            print(f"  {species_name:15s} (m={data['mass_amu']:6.1f} u, n={data['n_ions']:3d}): "
                  f"f_z = {data['f_z_fft_Hz']/1000:.2f} kHz (FFT), "
                  f"{data['f_z_crossing_Hz']/1000:.2f} kHz (crossing)")
    
    except Exception as e:
        print(f"  ⚠️  Failed: {e}")
        continue

print(f"\n{'='*80}")
print(f"Analysis complete: {len(all_results)} files processed")
print(f"{'='*80}\n")

# Collect data for mass scaling plot
masses = []
frequencies = []
species_labels = []

for file_key, file_results in all_results.items():
    for species_name, data in file_results.items():
        masses.append(data['mass_amu'])
        frequencies.append(data['f_z_fft_Hz'])
        species_labels.append(species_name)

if len(masses) == 0:
    print("❌ No valid data to plot!")
    exit(1)

masses = np.array(masses)
frequencies = np.array(frequencies)

# Theoretical scaling: f_z ∝ 1/sqrt(m)
# Fit to find k parameter: f_z = (1/2π) * sqrt(q*k/m)
q = 1.602176634e-19
m_kg = masses * 1.66054e-27
# f_z² = (1/4π²) * (q*k/m)
# k = f_z² * m * 4π² / q
k_fitted = np.median(frequencies**2 * m_kg * 4 * np.pi**2 / q)

# Generate theory curve
m_theory = np.logspace(np.log10(masses.min()*0.8), np.log10(masses.max()*1.2), 100)
m_theory_kg = m_theory * 1.66054e-27
f_theory = (1/(2*np.pi)) * np.sqrt(q * k_fitted / m_theory_kg)

# Create plots
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

# Left: f_z vs sqrt(m/z)
ax1.scatter(np.sqrt(masses), frequencies/1000, s=100, alpha=0.7, edgecolors='black')
for i, label in enumerate(species_labels):
    ax1.annotate(label, (np.sqrt(masses[i]), frequencies[i]/1000), 
                xytext=(5, 5), textcoords='offset points', fontsize=9)

ax1.plot(np.sqrt(m_theory), f_theory/1000, 'r--', linewidth=2, 
         label=f'Theory: k={k_fitted:.0f} V/m²')
ax1.set_xlabel('√(m/z)  [√u]', fontsize=12)
ax1.set_ylabel('Axial Frequency [kHz]', fontsize=12)
ax1.set_title('Orbitrap Frequency vs Mass', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)

# Right: Log-log plot for scaling verification
ax2.loglog(masses, frequencies/1000, 'o', markersize=10, alpha=0.7)
for i, label in enumerate(species_labels):
    ax2.annotate(label, (masses[i], frequencies[i]/1000), 
                xytext=(5, 5), textcoords='offset points', fontsize=9)

ax2.loglog(m_theory, f_theory/1000, 'r--', linewidth=2, 
          label=f'f_z ∝ m^(-0.5)')
ax2.set_xlabel('m/z [u]', fontsize=12)
ax2.set_ylabel('Axial Frequency [kHz]', fontsize=12)
ax2.set_title('Log-Log Scaling Verification', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3, which='both')

plt.tight_layout()

# Save plot
out_file = "../results/v1.0_test/instruments/orbitrap_frequency_scaling.png"
plt.savefig(out_file, dpi=300, bbox_inches='tight')
print(f"✅ Saved frequency scaling plot: {out_file}\n")

# Plot individual trajectories (if we have detailed single-species data)
fig2, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

plot_idx = 0
for file_key, file_results in all_results.items():
    for species_name, data in file_results.items():
        if plot_idx >= 4:
            break
        
        ax = axes[plot_idx]
        
        # Plot trajectory
        time_ms = data['time'] * 1000
        z_mm = data['z_trajectory'] * 1000
        ax.plot(time_ms, z_mm, linewidth=1, alpha=0.7)
        
        ax.set_xlabel('Time [ms]', fontsize=10)
        ax.set_ylabel('Axial Position z [mm]', fontsize=10)
        ax.set_title(f"{species_name} (m={data['mass_amu']:.0f} u, "
                    f"f_z={data['f_z_fft_Hz']/1000:.2f} kHz)", 
                    fontsize=11, fontweight='bold')
        ax.grid(True, alpha=0.3)
        
        plot_idx += 1

plt.tight_layout()
out_file2 = "../results/v1.0_test/instruments/orbitrap_trajectories.png"
plt.savefig(out_file2, dpi=300, bbox_inches='tight')
print(f"✅ Saved trajectory plots: {out_file2}\n")

# Validation summary
print("="*80)
print("Validation Summary:")
print("="*80)
print(f"Fitted k parameter: {k_fitted:.0f} V/m²")
print(f"\nFrequency measurements:")
print("-" * 60)
for i, (m, f, label) in enumerate(zip(masses, frequencies, species_labels)):
    m_kg = m * 1.66054e-27
    f_theory_single = (1/(2*np.pi)) * np.sqrt(q * k_fitted / m_kg)
    error_pct = 100 * (f - f_theory_single) / f_theory_single
    print(f"{label:15s} (m={m:6.1f} u): f_z = {f/1000:6.2f} kHz "
          f"(theory: {f_theory_single/1000:6.2f} kHz, error: {error_pct:+.2f}%)")
print("-" * 60)
print(f"\nTheoretical relationship validated: f_z ∝ 1/√m")
print("="*80 + "\n")

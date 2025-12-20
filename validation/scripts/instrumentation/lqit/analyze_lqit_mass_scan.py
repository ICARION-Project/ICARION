#!/usr/bin/env python3
"""
Analyze LQIT mass scan results

Tracks when ions are ejected during AC frequency sweep and validates
mass-dependent resonant ejection frequencies.
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
import sys

# Shared HDF5 helpers (species IDs)
COMMON_DIR = Path(__file__).resolve().parents[2] / "common"
if str(COMMON_DIR) not in sys.path:
    sys.path.append(str(COMMON_DIR))
from hdf5_utils import load_species_ids  # noqa: E402

def analyze_mass_scan(h5_file):
    """Analyze mass scan trajectory to find ejection frequencies"""
    
    with h5py.File(h5_file, 'r') as f:
        # Read trajectory data
        positions = f['/trajectory/positions'][:]  # (time, ion, xyz)
        times = f['/trajectory/time'][:]
        
        # Read domain parameters for each species
        domains = []
        for i in range(3):  # 3 domains/species
            domain_group = f[f'domains/domain_{i}']
            
            r0 = domain_group['geometry/radius_m'][()]
            L = domain_group['geometry/length_m'][()]
            origin = domain_group['geometry/origin_m'][:]
            
            rf_V = domain_group['fields/rf/voltage_V'][()]
            rf_f = domain_group['fields/rf/frequency_Hz'][()]
            
            # AC sweep parameters (read from inline waveform)
            ac_V = domain_group['fields/ac/voltage_V'][()]
            
            domains.append({
                'index': i,
                'r0': r0,
                'L': L,
                'origin': origin,
                'rf_V': rf_V,
                'rf_f': rf_f,
                'ac_V': ac_V,
                'z_min': origin[2],
                'z_max': origin[2] + L
            })
        
        # Read ion metadata
        species_names = load_species_ids(f)
        if species_names.ndim == 2:
            species_names = species_names[0, :]
        species_masses = f['metadata/species/mass_kg'][:]
        species_ids = [s.decode('utf-8') if isinstance(s, (bytes, bytearray)) else str(s) for s in f['metadata/species/names'][:]]
    
    # Group ions by species
    unique_species = {}
    for i, species in enumerate(species_names):
        if species not in unique_species:
            # Find which domain this species belongs to (match by index)
            species_idx = species_ids.index(species)
            domain_idx = species_idx  # Assuming domain_i corresponds to species_i
            
            unique_species[species] = {
                'indices': [],
                'domain': domains[domain_idx],
                'mass_kg': species_masses[species_idx]
            }
        unique_species[species]['indices'].append(i)
    
    print("\n" + "="*80)
    print("LQIT Mass Scan Analysis")
    print("="*80)
    
    # AC frequency sweep parameters (hardcoded, should match config)
    f_start = 30e3  # Hz
    f_end = 500e3   # Hz
    t_scan = times[-1]
    
    print(f"AC Sweep: {f_start/1e3:.0f} - {f_end/1e3:.0f} kHz over {t_scan*1e6:.0f} µs")
    print(f"Sweep rate: {(f_end - f_start)/t_scan/1e9:.2f} GHz/s")
    print()
    
    results = []
    
    for species_name, data in unique_species.items():
        ion_indices = data['indices']
        domain = data['domain']
        mass_kg = data['mass_kg']
        mass_amu = mass_kg / 1.66054e-27
        
        # Calculate theoretical secular frequency
        q_C = 1.602176634e-19
        omega_rf = 2 * np.pi * domain['rf_f']
        q_param = (4 * q_C * domain['rf_V']) / (mass_kg * omega_rf**2 * domain['r0']**2)
        omega_sec_theory = (q_param * omega_rf) / (2 * np.sqrt(2))
        f_sec_theory = omega_sec_theory / (2 * np.pi)
        
        # Track confinement for this species
        n_ions = len(ion_indices)
        confined = np.zeros((len(times), n_ions), dtype=bool)
        
        for i, ion_idx in enumerate(ion_indices):
            z = positions[:, ion_idx, 2]
            r = np.sqrt(positions[:, ion_idx, 0]**2 + positions[:, ion_idx, 1]**2)
            
            confined[:, i] = (z >= domain['z_min']) & (z <= domain['z_max']) & (r < domain['r0'])
        
        # Find ejection times (when confined drops below 50%)
        frac_confined = confined.sum(axis=1) / n_ions
        
        # Find first time when <50% confined
        ejection_idx = np.where(frac_confined < 0.5)[0]
        if len(ejection_idx) > 0:
            t_eject = times[ejection_idx[0]]
            # Calculate AC frequency at ejection time
            f_eject = f_start + (f_end - f_start) * (t_eject / t_scan)
            ejected = True
        else:
            t_eject = None
            f_eject = None
            ejected = False
        
        # Final confinement
        final_confined = confined[-1, :].sum()
        
        results.append({
            'species': species_name,
            'mass_amu': mass_amu,
            'q_param': q_param,
            'f_sec_theory': f_sec_theory,
            'n_ions': n_ions,
            'final_confined': final_confined,
            't_eject': t_eject,
            'f_eject': f_eject,
            'ejected': ejected,
            'frac_confined': frac_confined,
            'times': times
        })
    
    # Sort by mass
    results.sort(key=lambda x: x['mass_amu'])
    
    # Print results table
    print(f"{'Species':<20} {'m/z':>8} {'q':>8} {'f_sec(theo)':>12} {'f_eject':>12} {'Error':>10} {'Status'}")
    print("-"*80)
    
    for r in results:
        if r['ejected']:
            error = 100 * (r['f_eject'] - r['f_sec_theory']) / r['f_sec_theory']
            status = f"✅ EJECTED ({r['final_confined']}/{r['n_ions']})"
            f_eject_str = f"{r['f_eject']/1e3:6.1f} kHz"
            error_str = f"{error:+.1f}%"
        else:
            status = f"⚠️  CONFINED ({r['final_confined']}/{r['n_ions']})"
            f_eject_str = "N/A"
            error_str = "N/A"
        
        print(f"{r['species']:<20} {r['mass_amu']:8.1f} {r['q_param']:8.3f} "
              f"{r['f_sec_theory']/1e3:9.1f} kHz {f_eject_str:>12} {error_str:>10} {status}")
    
    print("="*80 + "\n")
    
    # Plot confinement vs time/frequency
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 8))
    
    # Plot 1: Confinement vs time
    for r in results:
        ax1.plot(r['times'] * 1e6, r['frac_confined'] * 100, 
                label=f"{r['species']} (m/z={r['mass_amu']:.0f})", linewidth=2)
    
    ax1.set_xlabel('Time [µs]', fontsize=12)
    ax1.set_ylabel('Ions Confined [%]', fontsize=12)
    ax1.set_ylim(-5, 105)
    ax1.grid(alpha=0.3)
    ax1.legend(fontsize=10)
    ax1.set_title('LQIT Mass Scan: Ion Confinement', fontsize=14, fontweight='bold')
    
    # Plot 2: Confinement vs AC frequency
    for r in results:
        # Calculate instantaneous AC frequency for each time point
        f_ac_instant = f_start + (f_end - f_start) * (r['times'] / t_scan)
        ax2.plot(f_ac_instant / 1e3, r['frac_confined'] * 100,
                label=f"{r['species']}", linewidth=2)
        
        # Mark theoretical resonance
        ax2.axvline(r['f_sec_theory'] / 1e3, color='gray', linestyle='--', alpha=0.5)
        ax2.text(r['f_sec_theory'] / 1e3, 95, f"{r['f_sec_theory']/1e3:.0f}",
                ha='center', fontsize=9, color='gray')
    
    ax2.set_xlabel('AC Excitation Frequency [kHz]', fontsize=12)
    ax2.set_ylabel('Ions Confined [%]', fontsize=12)
    ax2.set_ylim(-5, 105)
    ax2.set_xlim(f_start/1e3, f_end/1e3)
    ax2.grid(alpha=0.3)
    ax2.legend(fontsize=10)
    ax2.set_title('Mass-Selective Resonant Ejection', fontsize=14, fontweight='bold')
    
    plt.tight_layout()
    
    # Save plot
    out_dir = Path(h5_file).parent
    plot_file = out_dir / "lqit_mass_scan_analysis.png"
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"📊 Plot saved: {plot_file}")
    
    return results

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python3 analyze_lqit_mass_scan.py <h5_file>")
        sys.exit(1)
    
    h5_file = sys.argv[1]
    if not Path(h5_file).exists():
        print(f"Error: File not found: {h5_file}")
        sys.exit(1)
    
    analyze_mass_scan(h5_file)

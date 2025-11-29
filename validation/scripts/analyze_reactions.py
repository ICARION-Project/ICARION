#!/usr/bin/env python3
"""
Analyze reaction kinetics validation results

Validates:
1. First-order decay: N(t) = N₀ * exp(-k*t)
2. Bimolecular: N(t) = N₀ * exp(-k*[N2]*t)

Extracts rate constants from exponential fits.
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.optimize import curve_fit

# Physical constants
K_B = 1.380649e-23  # J/K

def read_species_counts(trajectory_file):
    """
    Read species counts vs time from HDF5 trajectory.
    
    Returns: (times, species_counts_dict)
    """
    with h5py.File(trajectory_file, 'r') as f:
        times = f['/trajectory/times'][:]
        species_ids = f['/trajectory/species_ids'][:]
        
        # Get unique species and count at each timestep
        unique_species = np.unique(species_ids)
        
        species_counts = {}
        for species in unique_species:
            species_name = species.decode() if isinstance(species, bytes) else species
            # Count at each timestep (assuming species_ids doesn't change during simulation)
            # Note: This is simplified - actual implementation depends on HDF5 structure
            count = np.sum(species_ids == species)
            species_counts[species_name] = count
        
        # For time evolution, we need to track active ions
        # This requires reading per-timestep data if available
        # Simplified: assume species_ids is constant and use trajectory length
        
    return times, species_counts

def exponential_decay(t, N0, k):
    """Exponential decay model: N(t) = N0 * exp(-k*t)"""
    return N0 * np.exp(-k * t)

def analyze_first_order(results_dir, k_expected_s):
    """Analyze first-order decay reaction"""
    results_dir = Path(results_dir)
    
    print(f"\nAnalyzing first-order decay: k = {k_expected_s:.1e} /s")
    print("-" * 60)
    
    traj_file = results_dir / f"first_order_decay_k{k_expected_s:.0e}.h5"
    
    if not traj_file.exists():
        print(f"❌ File not found: {traj_file}")
        return None
    
    try:
        with h5py.File(traj_file, 'r') as f:
            times = f['/trajectory/times'][:]
            species_ids = f['/trajectory/species_ids'][:]
            
            # Count ReactantA+ and ProductB+ at each timestep
            # NOTE: This assumes species tracking per timestep
            # Actual implementation depends on HDF5 structure
            
            # Simplified analysis: use final species distribution
            unique, counts = np.unique(species_ids, return_counts=True)
            species_dict = dict(zip(unique, counts))
            
            print(f"  Species distribution:")
            for species, count in species_dict.items():
                species_name = species.decode() if isinstance(species, bytes) else species
                print(f"    {species_name}: {count}")
            
            # For proper analysis, we need time-resolved species counts
            # This would require modifications to HDF5 output format
            
            print(f"  ⚠️  Time-resolved species tracking not yet implemented")
            print(f"  Expected: Exponential decay with τ = {1/k_expected_s*1e6:.1f} µs")
            
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return None
    
    return None

def analyze_bimolecular(results_dir, k_expected_cm3_s, pressure_Pa, temperature_K):
    """Analyze bimolecular reaction"""
    results_dir = Path(results_dir)
    
    # Calculate effective rate
    n_N2 = pressure_Pa / (K_B * temperature_K)
    k_eff_s = k_expected_cm3_s * 1e-6 * n_N2
    
    print(f"\nAnalyzing bimolecular: k = {k_expected_cm3_s:.1e} cm³/s")
    print(f"  [N2] = {n_N2:.2e} m⁻³")
    print(f"  k_eff = {k_eff_s:.2e} /s")
    print(f"  τ_eff = {1/k_eff_s*1e6:.1f} µs")
    print("-" * 60)
    
    traj_file = results_dir / f"bimolecular_k{k_expected_cm3_s:.0e}.h5"
    
    if not traj_file.exists():
        print(f"❌ File not found: {traj_file}")
        return None
    
    try:
        with h5py.File(traj_file, 'r') as f:
            times = f['/trajectory/times'][:]
            species_ids = f['/trajectory/species_ids'][:]
            
            unique, counts = np.unique(species_ids, return_counts=True)
            species_dict = dict(zip(unique, counts))
            
            print(f"  Species distribution:")
            for species, count in species_dict.items():
                species_name = species.decode() if isinstance(species, bytes) else species
                print(f"    {species_name}: {count}")
            
            print(f"  ⚠️  Time-resolved species tracking not yet implemented")
            print(f"  Expected: Pseudo-first-order decay")
            
    except Exception as e:
        print(f"  ❌ Error: {e}")
        return None
    
    return None

def analyze_all_reactions(results_dir):
    """Analyze all reaction kinetics results"""
    results_dir = Path(results_dir)
    
    print("\n" + "="*80)
    print("Reaction Kinetics Validation Analysis")
    print("="*80)
    print(f"Results directory: {results_dir}")
    print("="*80)
    
    # First-order reactions
    print("\n" + "="*80)
    print("FIRST-ORDER DECAY: A+ → B+")
    print("="*80)
    
    first_order_rates = [1e3, 5e3, 1e4]
    for k_s in first_order_rates:
        analyze_first_order(results_dir, k_s)
    
    # Bimolecular reactions
    print("\n" + "="*80)
    print("BIMOLECULAR: A+ + N2 → B+ + N2")
    print("="*80)
    
    bimolecular_rates = [1e-10, 5e-10, 1e-9]
    pressure_Pa = 1000.0
    temperature_K = 300.0
    
    for k_cm3_s in bimolecular_rates:
        analyze_bimolecular(results_dir, k_cm3_s, pressure_Pa, temperature_K)
    
    # Summary
    print("\n" + "="*80)
    print("VALIDATION NOTE")
    print("="*80)
    print("⚠️  Proper validation requires time-resolved species tracking!")
    print("\nCurrent HDF5 structure stores:")
    print("  - /trajectory/species_ids: [ions] (static)")
    print("\nRequired for reaction validation:")
    print("  - /trajectory/species_ids: [time, ions] (dynamic)")
    print("  OR")
    print("  - /species_counts: [time, species] (aggregated)")
    print("\nAlternative: Use reaction logger output if available")
    print("="*80 + "\n")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        results_dir = Path(sys.argv[1])
    else:
        results_dir = Path("../results/v1.0_test/physics/reactions")
    
    analyze_all_reactions(results_dir)

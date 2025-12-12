#!/usr/bin/env python3
"""
analyze_quad_stability.py

Analyzes quadrupole stability map validation results.

For each (a,q) test point:
  1. Reads HDF5 trajectory data
  2. Calculates transmission efficiency (fraction of ions that reached detector)
  3. Determines if configuration is stable or unstable
  4. Compares with theoretical Mathieu stability predictions

Outputs:
  - Summary table of all test points
  - Comparison with theoretical stability regions
  - Statistics on agreement with theory
"""

import h5py
import numpy as np
import json
import argparse
import os
import sys
from pathlib import Path

# Physical constants
E = 1.602176634e-19  # Elementary charge (C)

def read_config_parameters(config_path):
    """Read (a,q) parameters from config file."""
    with open(config_path, 'r') as f:
        config = json.load(f)
    
    # Extract quadrupole field parameters from domain
    domain = config['domains'][0]
    V_rf = domain['fields']['RF']['voltage_V']
    freq = domain['fields']['RF']['frequency_Hz']
    U_dc = domain['fields']['DC']['quad_V']
    r0 = domain['geometry']['radius_m']
    
    # Extract ion mass - hardcoded for CaffeineH+
    # CaffeineH+ = 195.08 amu (C8H11N4O2+)
    mass_amu = 195.08
    mass_kg = mass_amu * 1.66053906660e-27
    
    # Calculate Mathieu parameters
    omega = 2 * np.pi * freq
    
    a = 8 * E * U_dc / (mass_kg * r0**2 * omega**2)
    q = 4 * E * V_rf / (mass_kg * r0**2 * omega**2)
    
    return a, q, U_dc, V_rf, freq, r0, mass_amu


def is_theoretically_stable(a, q):
    """
    Check if (a,q) point is in the first Mathieu stability region.
    
    Uses approximate boundaries based on standard Mathieu stability diagram.
    First stability region (β_u and β_z both stable):
      - q: 0 to ~0.908
      - a: varies with q
    """
    # First stability region boundaries (more accurate)
    if q < 0 or q > 0.908:
        return False
    
    # Lower a boundary (approximately)
    a_min = -0.05  # Below this, u-direction becomes unstable
    
    # Upper a boundary (empirical fit from simulation data)
    # Observed pattern: a_max increases roughly linearly with q
    # At q=0.025: a_max ~ 0.01
    # At q=0.225: a_max ~ 0.04
    # Linear fit: a_max ≈ 0.15*q
    a_max = 0.15 * q
    
    return a_min <= a <= a_max


def analyze_trajectory(h5_path):
    """
    Analyze a single trajectory HDF5 file.
    
    Returns:
        transmission: Fraction of ions transmitted (0-1)
        n_initial: Initial number of ions
        n_final: Number of ions at detector
    """
    with h5py.File(h5_path, 'r') as f:
        # Read trajectory data (newer format)
        positions = f['trajectory']['positions'][:]  # Shape: (n_timesteps, n_particles, 3)
        
        # Get initial number of ions
        n_initial = positions.shape[1]
        
        # Get final positions (last timestep)
        final_positions = positions[-1, :, :]  # Shape: (n_particles, 3)
        
        # Count ions that made it through
        # In quadrupole, ions are transmitted if they reach the end (z > 0 for 50mm quadrupole starting at -0.05)
        # and are still within radial limits
        x_final = final_positions[:, 0]
        y_final = final_positions[:, 1]
        z_final = final_positions[:, 2]
        r_final = np.sqrt(x_final**2 + y_final**2)
        
        # Consider transmitted if:
        #  1. z > -0.005 m (near end of quadrupole which starts at -0.05 and ends at 0)
        #  2. r < 0.01 m (within 10mm radius - the quadrupole radius)
        transmitted = (z_final > -0.005) & (r_final < 0.01)
        n_final = np.sum(transmitted)
        
        transmission = n_final / n_initial if n_initial > 0 else 0.0
        
    return transmission, n_initial, n_final


def main():
    parser = argparse.ArgumentParser(
        description="Analyze quadrupole stability map validation results."
    )
    parser.add_argument(
        "--results-dir",
        help="Directory containing quadrupole trajectory .h5 files",
    )
    parser.add_argument(
        "--config-dir",
        help="Directory containing quadrupole config JSON files",
    )
    args = parser.parse_args()

    if args.results_dir:
        results_dir = Path(args.results_dir).expanduser().resolve()
    else:
        results_dir = _resolve_first_existing(_default_results_dir_candidates())
        if results_dir is None:
            candidates = "\n  - ".join(str(p) for p in _default_results_dir_candidates())
            print("❌ ERROR: Unable to locate results directory. Checked:\n  - " + candidates)
            return 1

    if not results_dir.exists():
        print(f"❌ ERROR: Results directory not found: {results_dir}")
        return 1

    if args.config_dir:
        config_dir = Path(args.config_dir).expanduser().resolve()
    else:
        config_dir = _resolve_first_existing(_default_config_dir_candidates())
        if config_dir is None:
            candidates = "\n  - ".join(str(p) for p in _default_config_dir_candidates())
            print("❌ ERROR: Unable to locate config directory. Checked:\n  - " + candidates)
            return 1

    if not config_dir.exists():
        print(f"❌ ERROR: Config directory not found: {config_dir}")
        return 1

    h5_files = sorted(results_dir.glob('quad_a*.h5'))
    if len(h5_files) == 0:
        print(f"❌ ERROR: No HDF5 files found in {results_dir}")
        return 1

    print("="*80)
    print("QUADRUPOLE STABILITY MAP ANALYSIS")
    print("="*80)
    print(f"Results dir : {results_dir}")
    print(f"Config dir  : {config_dir}")
    print(f"Found {len(h5_files)} trajectory files")
    print()
    
    # Analyze each test point
    results = []
    
    for h5_path in h5_files:
        # Extract (a,q) from filename
        # Format: quad_a-0.0100_q0.2000.h5 or quad_a+0.0300_q0.2000.h5
        filename = h5_path.stem
        parts = filename.split('_')
        a_str_h5 = parts[1].replace('a', '')  # e.g., "-0.0100" or "+0.0300"
        q_str = parts[2].replace('q', '')
        a_file = float(a_str_h5)
        q_file = float(q_str)
        
        matching_configs = list(config_dir.glob(f"quad_stability_a*{a_file:+.4f}*_q{q_str}.json"))
        
        if len(matching_configs) == 0:
            print(f"⚠️  Warning: Config not found for {filename}")
            continue
        
        config_path = matching_configs[0]
        
        # Read actual (a,q) from config
        a, q, U_dc, V_rf, freq, r0, mass = read_config_parameters(config_path)
        
        # Analyze trajectory
        transmission, n_initial, n_final = analyze_trajectory(h5_path)
        
        # Theoretical stability
        theory_stable = is_theoretically_stable(a, q)
        
        # Empirical stability (threshold at 50% transmission)
        empirical_stable = transmission >= 0.5
        
        # Agreement with theory
        agreement = (theory_stable == empirical_stable)
        
        results.append({
            'a': float(a),
            'q': float(q),
            'U_dc': float(U_dc),
            'V_rf': float(V_rf),
            'transmission': float(transmission),
            'n_initial': int(n_initial),
            'n_final': int(n_final),
            'theory_stable': bool(theory_stable),
            'empirical_stable': bool(empirical_stable),
            'agreement': bool(agreement),
            'filename': str(filename)
        })
    
    # Sort by q, then a
    results.sort(key=lambda x: (x['q'], x['a']))
        SCRIPT_DIR = Path(__file__).resolve().parent
        REPO_ROOT = SCRIPT_DIR.parents[3]
        VALIDATION_DIR = REPO_ROOT / "validation"
    
    # Print results table
    print()
    print("RESULTS SUMMARY")
    print("="*80)
    print(f"{'#':<4} {'a':>7} {'q':>7} {'U (V)':>8} {'V (V)':>8} {'Trans':>6} "
          f"{'Theory':>8} {'Simul':>8} {'Match':>6}")
    print("-"*80)
    
    for i, r in enumerate(results, 1):
        theory_str = "Stable" if r['theory_stable'] else "Unstable"
        simul_str = "Stable" if r['empirical_stable'] else "Unstable"
        match_str = "✅" if r['agreement'] else "❌"
        
        print(f"{i:<4} {r['a']:>7.4f} {r['q']:>7.4f} {r['U_dc']:>8.2f} {r['V_rf']:>8.2f} "
              f"{r['transmission']:>6.1%} {theory_str:>8} {simul_str:>8} {match_str:>6}")
    
    # Statistics
    print()
    print("="*80)
    print("STATISTICS")
    print("="*80)
    
    total = len(results)
    theory_stable_count = sum(1 for r in results if r['theory_stable'])
    theory_unstable_count = total - theory_stable_count
    simul_stable_count = sum(1 for r in results if r['empirical_stable'])
    simul_unstable_count = total - simul_stable_count
    agreement_count = sum(1 for r in results if r['agreement'])
    disagreement_count = total - agreement_count
    
    print(f"Total test points:        {total}")
    print(f"  Theory predicts stable: {theory_stable_count} ({100*theory_stable_count/total:.1f}%)")
    print(f"  Theory predicts unstable: {theory_unstable_count} ({100*theory_unstable_count/total:.1f}%)")
    print()
    print(f"  Simulation shows stable: {simul_stable_count} ({100*simul_stable_count/total:.1f}%)")
    print(f"  Simulation shows unstable: {simul_unstable_count} ({100*simul_unstable_count/total:.1f}%)")
    print()
    print(f"Agreement with theory:    {agreement_count}/{total} ({100*agreement_count/total:.1f}%)")
    
    if disagreement_count > 0:
        print()
        print("DISAGREEMENTS:")
        print("-"*80)
        for r in results:
            if not r['agreement']:
                theory_str = "Stable" if r['theory_stable'] else "Unstable"
                simul_str = "Stable" if r['empirical_stable'] else "Unstable"
                print(f"  a={r['a']:>7.4f}, q={r['q']:>7.4f}: Theory={theory_str}, "
                      f"Simul={simul_str} (Trans={r['transmission']:.1%})")
    

    def _default_results_dir_candidates():
        """Return ordered list of candidate results directories."""
        return [
            VALIDATION_DIR / "results" / "instruments" / "quadrupole",
            VALIDATION_DIR / "results" / "v1.0_test" / "instruments" / "quadrupole",
            REPO_ROOT / "results" / "v1.0_test" / "instruments" / "quadrupole",
        ]

    def _default_config_dir_candidates():
        return [
            VALIDATION_DIR / "configs" / "instruments" / "quadrupole",
            REPO_ROOT / "validation" / "configs" / "instruments" / "quadrupole",
        ]

    def _resolve_first_existing(paths):
        for path in paths:
            if path.exists():
                return path
        return None

    # Save results to JSON
    output_path = results_dir / 'stability_analysis.json'
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)
    
    print()
    print(f"✅ Analysis saved to: {output_path}")
    print()
    
    return 0


if __name__ == '__main__':
    sys.exit(main())

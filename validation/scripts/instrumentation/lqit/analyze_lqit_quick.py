#!/usr/bin/env python3
import h5py
import numpy as np
from pathlib import Path

results_dir = Path('results/v1.0_test/instruments/lqit')
h5_files = sorted([f for f in results_dir.glob('lqit_*.h5') if 'mass_scan' not in f.name])

print("="*80)
print("LQIT STABILITY ANALYSIS")
print("="*80)
print(f"Found {len(h5_files)} test files\n")

for h5_path in h5_files:
    with h5py.File(h5_path, 'r') as f:
        # Read basic info
        pos = f['positions']
        n_initial = pos.shape[1]
        
        # Count surviving ions
        final_pos = pos[-1, :, :]
        n_final = np.sum(~np.isnan(final_pos[:, 0]))
        retention = n_final / n_initial * 100
        
        # Get RF parameters from attrs
        q = f.attrs.get('q_parameter', 0.0)
        
        # Status
        status = "✅ STABLE" if retention >= 50 else "❌ UNSTABLE"
        
        print(f"{h5_path.name:<45} q={q:.3f}  {n_final:4d}/{n_initial:4d} ions ({retention:5.1f}%) {status}")

print("\n✅ Analysis complete!")

#!/usr/bin/env python3
"""
Analyze Quadrupole Stability Map

Reads trajectory data from (a,q) parameter sweep and generates stability map.
Compares with theoretical Mathieu stability boundaries.
"""

import numpy as np
import h5py
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.patches import Polygon
import json

# Data directory
DATA_DIR = Path("../results/v1.0.0_test/instruments/quadrupole_stability")

def load_simulation_data(h5_file):
    """Load final ion state and measure confinement
    
    Uses last timestep from trajectory data to determine stability.
    """
    with h5py.File(h5_file, 'r') as f:
        # Read geometry for bounds
        origin = f['domains/domain_0/geometry/origin_m'][:]
        r0 = f['domains/domain_0/geometry/radius_m'][()]
        L = f['domains/domain_0/geometry/length_m'][()]
        
        # Read trajectory (shape: [time, ions, xyz])
        positions = f['/trajectory/positions'][:]
        time = f['/trajectory/time'][:]
        
        # Extract final positions (last timestep)
        x = positions[-1, :, 0]
        y = positions[-1, :, 1]
        z = positions[-1, :, 2]
        
        # Check confinement
        r = np.sqrt(x**2 + y**2)
        z_min = origin[2]
        z_max = origin[2] + L
        
        confined = (r < r0) & (z >= z_min) & (z <= z_max)
        confined_pct = 100 * np.sum(confined) / len(confined)
        
        # Also check if any ions reached large radius (clear instability)
        max_r = np.max(r)
        escaped = max_r > 2*r0
        
        return {
            'confined_pct': confined_pct,
            'max_r': max_r,
            'escaped': escaped,
            'final_time_us': time[-1] * 1e6
        }

def theoretical_stability_boundary():
    """Generate theoretical stability boundaries for first region
    
    These are approximate - exact boundaries require solving Mathieu equations
    """
    # Lower q boundary (always a=0 line)
    q_lower = np.array([0.0, 0.908])
    a_lower = np.array([0.0, 0.0])
    
    # Upper boundary (empirical approximation)
    q_upper = np.linspace(0, 0.908, 50)
    # From Mathieu stability chart literature:
    # βz = 0 boundary approximately: a ≈ -q²/2 + q/sqrt(2)
    # βz = 1 boundary approximately: a ≈ q²/2 - q/sqrt(2)
    a_upper_pos = 0.706*q_upper - 0.708*q_upper**2  # Upper boundary (approximate)
    a_upper_neg = -0.706*q_upper + 0.708*q_upper**2  # Lower boundary (approximate)
    
    return {
        'q_boundary': q_upper,
        'a_upper': a_upper_pos,
        'a_lower': a_upper_neg
    }

print("\n" + "="*80)
print("Quadrupole Stability Map Analysis")
print("="*80 + "\n")

if not DATA_DIR.exists():
    print(f"❌ Data directory not found: {DATA_DIR}")
    print("Run simulations first!")
    exit(1)

# Find all trajectory files
h5_files = sorted(DATA_DIR.glob("quad_q*.h5"))

if len(h5_files) == 0:
    print(f"❌ No trajectory files found in {DATA_DIR}")
    exit(1)

print(f"Found {len(h5_files)} trajectory files\n")

# Parse (a,q) from filenames and analyze
results = []
for h5_file in h5_files:
    # Parse: quad_q0.400_a+0.050_stable.h5 or quad_q0.400_a-0.050_unstable.h5
    stem = h5_file.stem  # Remove .h5
    parts = stem.split('_')
    
    q_str = parts[1]  # q0.400
    a_str = parts[2]  # a+0.050 or a-0.050
    
    q = float(q_str[1:])  # Remove 'q'
    a = float(a_str[1:])  # Remove 'a'
    
    # Load data
    try:
        data = load_simulation_data(h5_file)
        
        # Classify: stable if >50% confined at end
        is_stable = data['confined_pct'] > 50.0
        
        results.append({
            'q': q,
            'a': a,
            'confined_pct': data['confined_pct'],
            'max_r': data['max_r'],
            'stable': is_stable,
            'escaped': data['escaped']
        })
        
        status = "✅ STABLE" if is_stable else "❌ UNSTABLE"
        print(f"q={q:.3f} a={a:+.3f}: {data['confined_pct']:5.1f}% confined → {status}")
        
    except Exception as e:
        print(f"⚠️  Failed to process {h5_file.name}: {e}")
        continue

print(f"\n{'='*80}")
print(f"Analysis complete: {len(results)} points")
print(f"Stable:   {sum(r['stable'] for r in results)}")
print(f"Unstable: {sum(not r['stable'] for r in results)}")
print(f"{'='*80}\n")

# Convert to arrays for plotting
results_arr = np.array([(r['q'], r['a'], r['confined_pct']) for r in results])
q_sim = results_arr[:, 0]
a_sim = results_arr[:, 1]
confined = results_arr[:, 2]

# Generate stability map plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

# Left: Binary stability map with theory overlay
stable_mask = confined > 50
unstable_mask = confined <= 50

ax1.scatter(q_sim[stable_mask], a_sim[stable_mask], 
           c='green', marker='o', s=100, label='Stable (sim)', alpha=0.7, edgecolors='black')
ax1.scatter(q_sim[unstable_mask], a_sim[unstable_mask], 
           c='red', marker='x', s=100, label='Unstable (sim)', alpha=0.7)

# Overlay theoretical boundaries
theory = theoretical_stability_boundary()
ax1.plot(theory['q_boundary'], theory['a_upper'], 'b--', linewidth=2, label='Theory (approx)')
ax1.plot(theory['q_boundary'], theory['a_lower'], 'b--', linewidth=2)
ax1.axvline(0.908, color='blue', linestyle='--', linewidth=2)

# Shade theoretical stable region
stable_region = np.column_stack([
    theory['q_boundary'],
    theory['a_lower']
])
stable_region = np.vstack([
    stable_region,
    np.column_stack([theory['q_boundary'][::-1], theory['a_upper'][::-1]])
])
poly = Polygon(stable_region, alpha=0.1, facecolor='green', edgecolor='none')
ax1.add_patch(poly)

ax1.set_xlabel('q parameter', fontsize=12)
ax1.set_ylabel('a parameter', fontsize=12)
ax1.set_title('Quadrupole Stability Map (Binary Classification)', fontsize=14, fontweight='bold')
ax1.legend(fontsize=10)
ax1.grid(True, alpha=0.3)
ax1.set_xlim(-0.05, 1.1)
ax1.set_ylim(-0.25, 0.25)

# Right: Continuous confinement percentage
scatter = ax2.scatter(q_sim, a_sim, c=confined, cmap='RdYlGn', 
                     s=100, vmin=0, vmax=100, edgecolors='black')
ax2.plot(theory['q_boundary'], theory['a_upper'], 'b--', linewidth=2, label='Theory')
ax2.plot(theory['q_boundary'], theory['a_lower'], 'b--', linewidth=2)
ax2.axvline(0.908, color='blue', linestyle='--', linewidth=2)

cbar = plt.colorbar(scatter, ax=ax2)
cbar.set_label('Confined %', fontsize=10)

ax2.set_xlabel('q parameter', fontsize=12)
ax2.set_ylabel('a parameter', fontsize=12)
ax2.set_title('Confinement Percentage Map', fontsize=14, fontweight='bold')
ax2.legend(fontsize=10)
ax2.grid(True, alpha=0.3)
ax2.set_xlim(-0.05, 1.1)
ax2.set_ylim(-0.25, 0.25)

plt.tight_layout()

# Save plot
out_file = "../results/v1.0.0_test/instruments/quadrupole_stability_map.png"
plt.savefig(out_file, dpi=300, bbox_inches='tight')
print(f"✅ Saved stability map: {out_file}\n")

# Generate scan line plots (a=0 and q=0.4)
fig2, (ax3, ax4) = plt.subplots(1, 2, figsize=(14, 5))

# a=0 scan (RF-only line)
a0_mask = np.abs(a_sim) < 0.01
if np.any(a0_mask):
    q_a0 = q_sim[a0_mask]
    conf_a0 = confined[a0_mask]
    sort_idx = np.argsort(q_a0)
    
    ax3.plot(q_a0[sort_idx], conf_a0[sort_idx], 'o-', linewidth=2, markersize=8)
    ax3.axvline(0.908, color='red', linestyle='--', linewidth=2, label='Theory boundary')
    ax3.axhline(50, color='gray', linestyle=':', linewidth=1, label='50% threshold')
    ax3.set_xlabel('q parameter', fontsize=12)
    ax3.set_ylabel('Confined %', fontsize=12)
    ax3.set_title('Scan along a=0 (RF-only)', fontsize=13, fontweight='bold')
    ax3.grid(True, alpha=0.3)
    ax3.legend()
    ax3.set_ylim(-5, 105)

# q≈0.4 scan (vary DC)
q04_mask = np.abs(q_sim - 0.4) < 0.05
if np.any(q04_mask):
    a_q04 = a_sim[q04_mask]
    conf_q04 = confined[q04_mask]
    sort_idx = np.argsort(a_q04)
    
    ax4.plot(a_q04[sort_idx], conf_q04[sort_idx], 'o-', linewidth=2, markersize=8)
    ax4.axhline(50, color='gray', linestyle=':', linewidth=1, label='50% threshold')
    ax4.set_xlabel('a parameter', fontsize=12)
    ax4.set_ylabel('Confined %', fontsize=12)
    ax4.set_title('Scan along q≈0.4 (vary DC)', fontsize=13, fontweight='bold')
    ax4.grid(True, alpha=0.3)
    ax4.legend()
    ax4.set_ylim(-5, 105)

plt.tight_layout()
out_file2 = "../results/v1.0.0_test/instruments/quadrupole_scan_lines.png"
plt.savefig(out_file2, dpi=300, bbox_inches='tight')
print(f"✅ Saved scan line plots: {out_file2}\n")

print("="*80)
print("Validation Summary:")
print("="*80)
print(f"Grid coverage: q ∈ [{q_sim.min():.2f}, {q_sim.max():.2f}], "
      f"a ∈ [{a_sim.min():.2f}, {a_sim.max():.2f}]")
print(f"Stable points:   {np.sum(stable_mask)} / {len(results)}")
print(f"Unstable points: {np.sum(unstable_mask)} / {len(results)}")
print(f"\nTheoretical first stability region: 0 < q < 0.908, |a| < ~0.2")
print("Compare simulated boundaries with theoretical curves in plots!")
print("="*80 + "\n")

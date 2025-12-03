#!/usr/bin/env python3
"""
Fit the empirical stability boundary from actual simulation data.
"""

import json
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy.interpolate import interp1d

# Load results
results_dir = Path(__file__).parent.parent.parent.parent / 'results' / 'v1.0_test' / 'instruments' / 'quadrupole'
with open(results_dir / 'stability_analysis.json', 'r') as f:
    results = json.load(f)

# Extract data
data = [(r['q'], r['a'], r['transmission']) for r in results]
data.sort()

# Group by q value and find stability boundary (50% transmission)
q_values = sorted(set(r['q'] for r in results))
boundary_points = []

for q in q_values:
    # Get all a values for this q
    q_data = [(r['a'], r['transmission']) for r in results if r['q'] == q]
    q_data.sort()
    
    # Find where transmission crosses 50%
    a_vals = np.array([x[0] for x in q_data])
    trans = np.array([x[1] for x in q_data])
    
    # Find upper boundary (where transmission drops below 50%)
    stable_mask = trans >= 0.5
    if np.any(stable_mask):
        # Get highest a value with >50% transmission
        stable_a = a_vals[stable_mask]
        unstable_a = a_vals[~stable_mask]
        
        if len(stable_a) > 0:
            max_stable_a = np.max(stable_a)
            
            # Try to interpolate between last stable and first unstable
            if len(unstable_a) > 0 and np.min(unstable_a) > max_stable_a:
                # There's a boundary
                next_unstable_a = np.min(unstable_a[unstable_a > max_stable_a])
                
                # Get transmissions
                trans_stable = trans[a_vals == max_stable_a][0]
                trans_unstable = trans[a_vals == next_unstable_a][0]
                
                # Linear interpolation to 50%
                if trans_stable != trans_unstable:
                    a_boundary = max_stable_a + (next_unstable_a - max_stable_a) * (0.5 - trans_stable) / (trans_unstable - trans_stable)
                else:
                    a_boundary = (max_stable_a + next_unstable_a) / 2
            else:
                # Boundary is at or above max tested a
                a_boundary = max_stable_a
            
            boundary_points.append((q, a_boundary))
            print(f"q={q:.3f}: boundary at a={a_boundary:.4f}")

# Convert to arrays
q_boundary = np.array([p[0] for p in boundary_points])
a_boundary = np.array([p[1] for p in boundary_points])

print(f"\n✅ Found {len(boundary_points)} boundary points")
print("\nBoundary data:")
for q, a in boundary_points:
    print(f"  q={q:.4f} -> a_max={a:.4f}")

# Try different polynomial fits
print("\n" + "="*60)
print("FITTING DIFFERENT MODELS:")
print("="*60)

# Model 1: Linear
p_linear = np.polyfit(q_boundary, a_boundary, 1)
a_linear = np.poly1d(p_linear)
residuals_linear = a_boundary - a_linear(q_boundary)
rms_linear = np.sqrt(np.mean(residuals_linear**2))
print(f"\n1. Linear: a = {p_linear[0]:.4f}*q + {p_linear[1]:.4f}")
print(f"   RMS error: {rms_linear:.4f}")

# Model 2: Quadratic
p_quad = np.polyfit(q_boundary, a_boundary, 2)
a_quad = np.poly1d(p_quad)
residuals_quad = a_boundary - a_quad(q_boundary)
rms_quad = np.sqrt(np.mean(residuals_quad**2))
print(f"\n2. Quadratic: a = {p_quad[0]:.4f}*q² + {p_quad[1]:.4f}*q + {p_quad[2]:.4f}")
print(f"   RMS error: {rms_quad:.4f}")

# Model 3: Cubic
p_cubic = np.polyfit(q_boundary, a_boundary, 3)
a_cubic = np.poly1d(p_cubic)
residuals_cubic = a_boundary - a_cubic(q_boundary)
rms_cubic = np.sqrt(np.mean(residuals_cubic**2))
print(f"\n3. Cubic: a = {p_cubic[0]:.4f}*q³ + {p_cubic[1]:.4f}*q² + {p_cubic[2]:.4f}*q + {p_cubic[3]:.4f}")
print(f"   RMS error: {rms_cubic:.4f}")

# Model 4: Square root form (like theoretical)
from scipy.optimize import curve_fit

def sqrt_model(q, a0, q_max):
    """a_max = a0 * sqrt(1 - q/q_max)"""
    return a0 * np.sqrt(np.maximum(0, 1 - q/q_max))

try:
    popt_sqrt, _ = curve_fit(sqrt_model, q_boundary, a_boundary, p0=[0.2, 0.9], bounds=([0, 0.8], [0.5, 1.2]))
    a_sqrt = sqrt_model(q_boundary, *popt_sqrt)
    residuals_sqrt = a_boundary - a_sqrt
    rms_sqrt = np.sqrt(np.mean(residuals_sqrt**2))
    print(f"\n4. Square root: a = {popt_sqrt[0]:.4f} * sqrt(1 - q/{popt_sqrt[1]:.4f})")
    print(f"   RMS error: {rms_sqrt:.4f}")
except:
    print("\n4. Square root fit failed")
    rms_sqrt = float('inf')

# Choose best model
best_model = np.argmin([rms_linear, rms_quad, rms_cubic, rms_sqrt])
model_names = ["Linear", "Quadratic", "Cubic", "Square root"]
print(f"\n{'='*60}")
print(f"✅ BEST MODEL: {model_names[best_model]}")
print(f"{'='*60}")

# Create plot
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 7))

# Left plot: Full stability map with all models
q_vals = np.array([r['q'] for r in results])
a_vals = np.array([r['a'] for r in results])
trans = np.array([r['transmission'] for r in results])

scatter = ax1.scatter(q_vals, a_vals, c=trans*100, s=150, 
                     cmap='RdYlGn', vmin=0, vmax=100,
                     edgecolors='black', linewidth=1.5, zorder=10)

# Plot all fitted models
q_plot = np.linspace(0.05, 0.95, 200)
ax1.plot(q_plot, a_linear(q_plot), 'b--', linewidth=2, alpha=0.6, label=f'Linear (RMS={rms_linear:.3f})')
ax1.plot(q_plot, a_quad(q_plot), 'g--', linewidth=2, alpha=0.6, label=f'Quadratic (RMS={rms_quad:.3f})')
ax1.plot(q_plot, a_cubic(q_plot), 'm--', linewidth=2, alpha=0.6, label=f'Cubic (RMS={rms_cubic:.3f})')
if rms_sqrt != float('inf'):
    ax1.plot(q_plot, sqrt_model(q_plot, *popt_sqrt), 'r--', linewidth=2, alpha=0.6, label=f'Sqrt (RMS={rms_sqrt:.3f})')

# Plot boundary points
ax1.plot(q_boundary, a_boundary, 'ko', markersize=10, label='Boundary (50% trans)', zorder=15)

ax1.set_xlabel('q (Trapping Parameter)', fontsize=14, fontweight='bold')
ax1.set_ylabel('a (Stability Parameter)', fontsize=14, fontweight='bold')
ax1.set_title('Stability Boundary Fitting', fontsize=16, fontweight='bold')
ax1.set_xlim(-0.02, 1.05)
ax1.set_ylim(-0.06, 0.30)
ax1.grid(True, alpha=0.3)
ax1.legend(loc='upper left', fontsize=10)

cbar = plt.colorbar(scatter, ax=ax1, label='Transmission (%)')

# Right plot: Residuals
ax2.axhline(0, color='k', linestyle='-', linewidth=1)
ax2.plot(q_boundary, residuals_linear, 'bo-', label='Linear', markersize=6)
ax2.plot(q_boundary, residuals_quad, 'go-', label='Quadratic', markersize=6)
ax2.plot(q_boundary, residuals_cubic, 'mo-', label='Cubic', markersize=6)
if rms_sqrt != float('inf'):
    ax2.plot(q_boundary, residuals_sqrt, 'ro-', label='Sqrt', markersize=6)

ax2.set_xlabel('q', fontsize=12, fontweight='bold')
ax2.set_ylabel('Residual (a_measured - a_fit)', fontsize=12, fontweight='bold')
ax2.set_title('Fit Quality Comparison', fontsize=14, fontweight='bold')
ax2.grid(True, alpha=0.3)
ax2.legend()

plt.tight_layout()
plt.savefig(results_dir / 'stability_boundary_fit.png', dpi=300, bbox_inches='tight')
print(f"\n✅ Plot saved to: {results_dir / 'stability_boundary_fit.png'}")

# Save fit parameters
fit_params = {
    "boundary_points": [(float(q), float(a)) for q, a in boundary_points],
    "linear": {"coeffs": [float(c) for c in p_linear], "rms": float(rms_linear)},
    "quadratic": {"coeffs": [float(c) for c in p_quad], "rms": float(rms_quad)},
    "cubic": {"coeffs": [float(c) for c in p_cubic], "rms": float(rms_cubic)},
    "best_model": model_names[best_model]
}
if rms_sqrt != float('inf'):
    fit_params["sqrt"] = {"a0": float(popt_sqrt[0]), "q_max": float(popt_sqrt[1]), "rms": float(rms_sqrt)}

with open(results_dir / 'stability_fit_params.json', 'w') as f:
    json.dump(fit_params, f, indent=2)

print(f"✅ Fit parameters saved to: {results_dir / 'stability_fit_params.json'}")

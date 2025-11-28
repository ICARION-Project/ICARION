#!/usr/bin/env python3
"""
Quick IMS drift velocity analysis

Calculates drift velocity from position data and compares to Mason-Schamp prediction.
"""

import h5py
import numpy as np
import sys

# Mason-Schamp parameters
K0_H3Op_He = 10.5e-4  # m²/(V·s) - H3O+ in He reduced mobility
P0 = 101325.0  # Pa
T0 = 273.15  # K
T = 300.0  # K

def analyze_drift(h5_file, E_Vm, P_Pa):
    """Analyze drift velocity from HDF5 trajectory data"""
    
    with h5py.File(h5_file, 'r') as f:
        positions = f['/trajectory/positions'][:]  # (time, ion, xyz)
        times = f['/trajectory/time'][:]
        
    # Get z-positions (drift direction)
    z_pos = positions[:, :, 2]
    
    # Find ions that drifted (moved significantly in +z)
    z_start = z_pos[0, :]
    z_end = z_pos[-1, :]
    drift_distance = z_end - z_start
    
    # Filter ions that moved > 1mm in +z (actual drift, not noise)
    drifted_mask = drift_distance > 1e-3
    n_drifted = np.sum(drifted_mask)
    
    if n_drifted == 0:
        print(f"  ⚠️  No ions drifted > 1mm")
        return None
    
    # Calculate drift velocity for drifted ions
    t_total = times[-1] - times[0]
    drift_velocities = drift_distance[drifted_mask] / t_total
    
    v_drift_mean = np.mean(drift_velocities)
    v_drift_std = np.std(drift_velocities)
    
    # Mason-Schamp prediction
    N_ratio = (P0 / P_Pa) * (T / T0)
    v_expected = K0_H3Op_He * E_Vm * N_ratio
    
    error_pct = 100 * (v_drift_mean - v_expected) / v_expected
    
    print(f"  Drifted ions: {n_drifted}/{len(z_start)}")
    print(f"  Drift distance: {np.mean(drift_distance[drifted_mask])*1e3:.1f} ± {np.std(drift_distance[drifted_mask])*1e3:.1f} mm")
    print(f"  Drift velocity: {v_drift_mean:.1f} ± {v_drift_std:.1f} m/s")
    print(f"  Mason-Schamp:   {v_expected:.1f} m/s")
    print(f"  Error:          {error_pct:+.1f}%")
    
    return v_drift_mean, v_expected, error_pct

if __name__ == "__main__":
    print("=" * 70)
    print("IMS Drift Velocity Analysis")
    print("=" * 70)
    
    tests = [
        ("Friction", "ims_friction_5000Vm_100Pa.h5", 5000, 100),
        ("HSS", "ims_hss_5000Vm_100Pa.h5", 5000, 100),
        ("EHSS", "ims_ehss_5000Vm_100Pa.h5", 5000, 100),
    ]
    
    results = []
    for model, filename, E, P in tests:
        print(f"\n{model} @ E={E} V/m, P={P} Pa:")
        h5_path = f"../results/v1.0_test/instruments/ims/{filename}"
        try:
            result = analyze_drift(h5_path, E, P)
            if result:
                results.append((model, result))
        except Exception as e:
            print(f"  ❌ Error: {e}")
    
    print("\n" + "=" * 70)
    print("SUMMARY")
    print("=" * 70)
    for model, (v_meas, v_exp, err) in results:
        status = "✅" if abs(err) < 10 else "⚠️"
        print(f"{status} {model:10s}: {v_meas:6.1f} m/s (expected {v_exp:.1f} m/s, error {err:+.1f}%)")

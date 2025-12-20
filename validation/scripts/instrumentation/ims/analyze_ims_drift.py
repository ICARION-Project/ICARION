#!/usr/bin/env python3
"""
Quick IMS drift velocity analysis

Calculates drift velocity from position data and compares to Mason-Schamp prediction
with effective temperature correction.
"""

import h5py
import numpy as np
import sys

# Physical constants
kB = 1.380649e-23  # J/K - Boltzmann constant
amu = 1.66053906660e-27  # kg - atomic mass unit
e = 1.602176634e-19  # C - elementary charge

# Mason-Schamp parameters
K0_H3Op_He = 24.1e-4  # m²/(V·s) - H3O+ in He reduced mobility (from species_database_v1.json)
M_H3Op = 19.02 * amu  # kg - H3O+ mass
M_Hep = 4.003 * amu  # kg - He mass
LOSCHMIDT = 2.6867774e25  # m^-3 - number density at STP
T0 = 273.15  # K

def analyze_drift(h5_file):
    """Analyze drift velocity from HDF5 trajectory data"""
    
    with h5py.File(h5_file, 'r') as f:
        # Read simulation parameters from HDF5
        axial_V = f['domains/domain_0/fields/dc/axial_V'][()]
        length_m = f['domains/domain_0/geometry/length_m'][()]
        E_Vm = axial_V / length_m
        P_Pa = f['domains/domain_0/environment/pressure_Pa'][()]
        T_K = f['domains/domain_0/environment/temperature_K'][()]
        N_m3 = f['domains/domain_0/environment/particle_density_m3'][()]
        
        # Read trajectory data
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
    
    # Calculate E/N in Townsend units
    E_N_Td = (E_Vm / N_m3) * 1e21  # 1 Td = 1e-21 V·m²
    
    # Effective temperature correction (matches future damping force scaling)
    N_ratio = LOSCHMIDT / N_m3
    v_drift_field = K0_H3Op_He * E_Vm * N_ratio
    T_heating = (M_Hep / (3 * kB)) * v_drift_field**2
    T_eff = T_K + T_heating

    # Mason-Schamp with effective-temperature-adjusted mobility
    K_eff = K0_H3Op_He * np.sqrt(T0 / T_eff)
    v_expected = K_eff * E_Vm * N_ratio
    
    error_pct = 100 * (v_drift_mean - v_expected) / v_expected
    
    print(f"  E-field: {E_Vm:.1f} V/m @ P={P_Pa:.0f} Pa, T={T_K:.0f} K")
    print(f"  E/N: {E_N_Td:.1f} Td")
    print(f"  T_eff: {T_eff:.1f} K (T_heating: {T_heating:.1f} K)")
    print(f"  Drifted ions: {n_drifted}/{len(z_start)}")
    print(f"  Drift distance: {np.mean(drift_distance[drifted_mask])*1e3:.1f} ± {np.std(drift_distance[drifted_mask])*1e3:.1f} mm")
    print(f"  Drift velocity: {v_drift_mean:.1f} ± {v_drift_std:.1f} m/s")
    print(f"  Mason-Schamp (T_eff):   {v_expected:.1f} m/s")
    print(f"  Error:          {error_pct:+.1f}%")
    
    return v_drift_mean, v_expected, error_pct

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: analyze_ims_drift.py <h5_file1> [h5_file2] ...")
        sys.exit(1)
    
    print("=" * 70)
    print("IMS Drift Velocity Analysis")
    print("=" * 70)
    
    results = []
    for h5_file in sys.argv[1:]:
        # Extract model name from filename (e.g., ims_hss_966Vm_100Pa.h5 -> HSS)
        import os
        basename = os.path.basename(h5_file)
        model = basename.split('_')[1].upper()
        
        print(f"\n{model}:")
        try:
            result = analyze_drift(h5_file)
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

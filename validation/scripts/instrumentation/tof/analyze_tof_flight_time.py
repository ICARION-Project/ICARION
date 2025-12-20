#!/usr/bin/env python3
"""
Analyze Time-of-Flight (TOF) mass spectrometer flight times.

Theory:
  Ion acceleration: E_kin = q*V_acc = 0.5*m*v²
  Final velocity: v = sqrt(2*q*V_acc / m)
  
  Flight time:
    Acceleration phase (uniformly accelerated motion):
      t_acc = 2*L_acc / v  (average velocity = v/2)
    
    Drift phase (constant velocity):
      t_drift = L_drift / v
    
  Total: t = (2*L_acc + L_drift) / v
           = (2*L_acc + L_drift) * sqrt(m / (2*q*V_acc))
         
  NOT simply t = L_total * sqrt(m / (2*q*V_acc)) ← WRONG!
"""

import h5py
import numpy as np
import sys
from pathlib import Path
import json

# Shared HDF5 helpers (species IDs, embedded inputs)
COMMON_DIR = Path(__file__).resolve().parents[2] / "common"
if str(COMMON_DIR) not in sys.path:
    sys.path.append(str(COMMON_DIR))
from hdf5_utils import load_species_ids, load_config_json, load_species_db_json  # noqa: E402


def _decode_species_id(value):
    """Ensure species identifiers coming from HDF5 are strings."""
    if isinstance(value, bytes):
        return value.decode()
    return str(value)

def calculate_theoretical_flight_time(mass_u, charge_e, V_acc, L_acc, L_drift):
    """Calculate theoretical TOF flight time.
    
    Corrected formula accounting for acceleration phase:
      t_acc = 2*L_acc / v  (uniformly accelerated motion)
      t_drift = L_drift / v (constant velocity drift)
      t_total = (2*L_acc + L_drift) / v
    """
    q_C = charge_e * 1.602176634e-19  # C
    m_kg = mass_u * 1.66053906660e-27  # kg
    
    # Final velocity after acceleration: v = sqrt(2*q*V/m)
    v_final = np.sqrt(2 * q_C * V_acc / m_kg)
    
    # Acceleration phase: uniformly accelerated motion, average velocity = v/2
    # Distance L_acc traveled with average velocity v/2: t_acc = L_acc / (v/2) = 2*L_acc/v
    t_acc = 2 * L_acc / v_final
    
    # Drift phase: constant velocity v
    t_drift = L_drift / v_final
    
    # Total time: t = (2*L_acc + L_drift) / v = (2*L_acc + L_drift) * sqrt(m/(2*q*V))
    return t_acc + t_drift

def analyze_tof_trajectory(h5_path, config_path=None):
    """Extract flight times per species from a TOF HDF5 trajectory."""
    with h5py.File(h5_path, 'r') as f:
        positions = f['trajectory/positions'][:]  # (n_steps, n_ions, 3)
        time = f['trajectory/time'][:]  # (n_steps,)

        # Number of ions
        n_ions = positions.shape[1] if positions.ndim == 3 else 1
        if positions.ndim == 2:
            positions = positions[:, np.newaxis, :]

        raw_species = load_species_ids(f)
        if raw_species.ndim == 2:
            raw_species = raw_species[0, :]

    ion_species = [_decode_species_id(s) for s in raw_species]
    if len(ion_species) != n_ions:
        raise ValueError(f"Species metadata length ({len(ion_species)}) does not match ion count ({n_ions}).")

    # Read config
    if config_path is None:
        # Prefer embedded config if present
        with h5py.File(h5_path, 'r') as f:
            embedded_cfg = load_config_json(f)
        if embedded_cfg:
            config = embedded_cfg
            config_path = None
        else:
            config_path = str(h5_path).replace('.h5', '.json').replace('results/v1.0_test/instruments/', 'configs/instruments/')
    if config_path:
        with open(config_path, 'r') as f:
            config = json.load(f)

    domain0 = config['domains'][0]
    fields = domain0.get('fields', {})
    dc = fields.get('DC') or fields.get('dc') or {}
    V_acc = dc.get('axial_V')
    L_acc = domain0.get('geometry', {}).get('acc_length_m')
    L_total = domain0.get('geometry', {}).get('length_m')

    species_db = {}
    with h5py.File(h5_path, 'r') as f:
        embedded_species_db = load_species_db_json(f)
    if embedded_species_db:
        species_db = embedded_species_db.get('species', {})
    elif config.get('species_database_path') or config.get('species_database'):
        db_path = config.get('species_database_path') or config.get('species_database')
        with open(db_path, 'r') as f:
            raw_db = json.load(f)
        species_db = raw_db.get('species', {})

    species_configs = []
    for entry in config['ions']['species']:
        sid = entry.get('species_id', 'Unknown')
        db_entry = species_db.get(sid, {})
        mass_u = db_entry.get('mass_amu', entry.get('mass_Da', 0))
        charge = db_entry.get('charge', entry.get('charge', entry.get('charge_state', 1)))
        center = entry.get('position', {}).get('center', [0, 0, 0])
        z_start = center[2] if len(center) > 2 else 0.0
        species_configs.append({
            'species_id': sid,
            'mass_u': mass_u,
            'charge': charge,
            'z_start': z_start
        })

    z_detector = 0.999  # m
    results = []

    for spec in species_configs:
        sid = spec['species_id']
        ion_indices = [idx for idx, label in enumerate(ion_species) if label == sid]
        if not ion_indices:
            continue

        flight_times = []
        for ion_idx in ion_indices:
            z_pos = positions[:, ion_idx, 2]
            crossed = np.where(z_pos >= z_detector)[0]
            if len(crossed) > 0:
                flight_times.append(time[crossed[0]])

        if not flight_times:
            # No ions of this species reached the detector
            continue

        t_mean = float(np.mean(flight_times) * 1e6)
        t_std = float(np.std(flight_times) * 1e6)

        z_start = spec['z_start']
        L_acc_eff = max(L_acc - z_start, 1e-12)
        V_eff = V_acc * L_acc_eff / L_acc
        L_drift = L_total - L_acc
        t_theory = calculate_theoretical_flight_time(
            spec['mass_u'],
            spec['charge'],
            V_eff,
            L_acc_eff,
            L_drift
        ) * 1e6

        error_pct = 100 * (t_mean - t_theory) / t_theory if t_theory else float('inf')

        results.append({
            'species': sid,
            'mass_u': spec['mass_u'],
            'charge': spec['charge'],
            'V_acc': V_acc,
            'L_total_m': L_total,
            't_theory_us': t_theory,
            't_mean_us': t_mean,
            't_std_us': t_std,
            'error_pct': error_pct,
            'n_ions': len(flight_times),
            'resolution': t_mean / t_std if t_std > 0 else np.inf
        })

    return results

def main():
    if len(sys.argv) < 2:
        print("Usage: python analyze_tof_flight_time.py <trajectory.h5> [<trajectory2.h5> ...]")
        sys.exit(1)
    
    results = []
    for h5_path in sys.argv[1:]:
        if not Path(h5_path).exists():
            print(f"⚠ File not found: {h5_path}")
            continue
        
        try:
            species_results = analyze_tof_trajectory(h5_path)
            if not species_results:
                print(f"⚠ No ions reached detector: {h5_path}")
                continue
            results.extend(species_results)
        except Exception as e:
            print(f"✗ Error analyzing {h5_path}: {e}")
            import traceback
            traceback.print_exc()
            continue
    
    if not results:
        print("No valid results to display")
        return
    
    # Print results table
    print("\n" + "="*110)
    print("TIME-OF-FLIGHT (TOF) VALIDATION")
    print("="*110)
    print(f"{'Species':<15} {'Mass [u]':<10} {'V_acc [V]':<12} {'t_theory [µs]':<15} {'t_meas [µs]':<15} {'σ [ns]':<10} {'Error [%]':<12} {'Status':<10}")
    print("-"*110)
    
    for r in results:
        status = "✅ PASS" if abs(r['error_pct']) < 1.0 else "⚠ WARN" if abs(r['error_pct']) < 5.0 else "✗ FAIL"
        sigma_ns = r['t_std_us'] * 1000  # convert µs to ns
        print(f"{r['species']:<15} {r['mass_u']:<10.2f} {r['V_acc']:<12.0f} {r['t_theory_us']:<15.3f} {r['t_mean_us']:<15.3f} {sigma_ns:<10.2f} {r['error_pct']:>+11.2f} {status:<10}")
    
    print("-"*110)
    print(f"\nFlight tube: L = {results[0]['L_total_m']*100:.0f} cm")
    print(f"Theory: t ∝ √m  →  Mass scaling test:")
    
    # Verify mass scaling: t1/t2 = sqrt(m1/m2)
    if len(results) >= 2:
        for i in range(len(results)-1):
            r1, r2 = results[i], results[i+1]
            ratio_measured = r1['t_mean_us'] / r2['t_mean_us']
            ratio_theory = np.sqrt(r1['mass_u'] / r2['mass_u'])
            ratio_error = 100 * (ratio_measured - ratio_theory) / ratio_theory
            print(f"  {r1['species']}/{r2['species']}: measured={ratio_measured:.4f}, theory={ratio_theory:.4f}, Δ={ratio_error:+.2f}%")
    
    print("\nMass resolution (m/Δm = t/Δt):")
    for r in results:
        print(f"  {r['species']}: {r['resolution']:.0f}")
    
    print("="*110)

if __name__ == '__main__':
    main()

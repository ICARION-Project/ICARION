#!/usr/bin/env python3
"""
Analyze transport drift velocity validation results

Extracts drift velocity from trajectories and validates:
  v_drift = K₀ × E
  
Compares measured K₀ with literature values.
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from scipy import stats

K_BOLTZMANN = 1.380649e-23  # J/K

# Literature reduced mobility values (cm²/(V·s) in He at 300K)
K0_LITERATURE = {
    "H3O+": 2.8,
    "PentanalH+": 1.9,
    "ReserpineH+": 0.85,
}


def _decode_if_bytes(value):
    if isinstance(value, (bytes, np.bytes_)):
        return value.decode('utf-8')
    return value


def _read_dataset(h5file, path, *, first=False, default=None):
    try:
        data = h5file[path][()]
    except KeyError:
        return default

    if isinstance(data, np.ndarray):
        if data.shape == ():
            data = data.item()
        elif first:
            data = data.reshape(-1)[0]
        else:
            data = data.tolist()

    if isinstance(data, list):
        return [_decode_if_bytes(v) for v in data]

    return _decode_if_bytes(data)


def _read_time_axis(h5file):
    if '/trajectory/time' in h5file:
        return h5file['/trajectory/time'][:]
    if '/trajectory/times' in h5file:
        return h5file['/trajectory/times'][:]
    raise KeyError("Trajectory file is missing '/trajectory/time' dataset")


def extract_metadata(h5file, filename):
    """Gather contextual metadata from the trajectory file."""
    species = _read_dataset(h5file, '/metadata/species/names', first=True, default='Unknown')
    model = _read_dataset(h5file, '/metadata/config/collision_model', default='Unknown')
    temp_K = _read_dataset(h5file, '/domains/domain_0/environment/temperature_K')
    pressure_Pa = _read_dataset(h5file, '/domains/domain_0/environment/pressure_Pa')
    gas_species = _read_dataset(h5file, '/domains/domain_0/environment/gas_species', first=True, default='Unknown')
    en_td = _read_dataset(h5file, '/domains/domain_0/fields/dc/EN_Td')
    axial_V = _read_dataset(h5file, '/domains/domain_0/fields/dc/axial_V')
    length_m = _read_dataset(h5file, '/domains/domain_0/geometry/length_m')

    computed_en_td = None
    try:
        if axial_V is not None and length_m and temp_K and pressure_Pa:
            e_field = float(axial_V) / float(length_m)
            number_density = float(pressure_Pa) / (K_BOLTZMANN * float(temp_K))
            if number_density > 0.0:
                computed_en_td = (e_field / number_density) / 1e-21
    except ZeroDivisionError:
        computed_en_td = None

    # Fall back to filename hints if metadata missing
    if temp_K is None:
        temp_K = _infer_numeric_from_name(filename, 'K')
    if en_td is None:
        en_td = computed_en_td if computed_en_td is not None else _infer_numeric_from_name(filename, 'Td')
    elif isinstance(en_td, (int, float)) and en_td == 0.0 and computed_en_td is not None:
        en_td = computed_en_td

    return {
        'species': species,
        'model': model,
        'temp_K': float(temp_K) if temp_K is not None else None,
        'en_td': float(en_td) if en_td is not None else None,
        'pressure_Pa': float(pressure_Pa) if pressure_Pa is not None else None,
        'gas_species': gas_species,
    }


def _infer_numeric_from_name(filename, suffix):
    try:
        stem = Path(filename).stem
        token = stem.split('_')[-1]
        if token.endswith(suffix):
            return float(token[:-len(suffix)])
    except Exception:
        return None
    return None


def extract_drift_velocity(h5file):
    """Extract drift velocity from a trajectory HDF5 file handle."""
    pos = h5file['/trajectory/positions'][:]  # [time, ions, xyz]
    times = _read_time_axis(h5file)

    # Extract z-coordinate (drift direction)
    z = pos[:, :, 2]  # [time, ions]

    # Calculate mean position vs time
    z_mean = np.mean(z, axis=1)

    # Use second half of trajectory (after initial transients)
    n_half = len(times) // 2
    t_fit = times[n_half:]
    z_fit = z_mean[n_half:]

    # Linear fit: z(t) = z₀ + v_drift · t
    slope, _, r_value, _, std_err = stats.linregress(t_fit, z_fit)

    v_drift = slope  # m/s
    v_drift_std = std_err  # Standard error

    return v_drift, v_drift_std, times[-1], r_value**2


def analyze_drift_results(results_dir):
    results_dir = Path(results_dir)

    print("\n" + "="*80)
    print("Transport Drift Velocity Analysis")
    print("="*80)
    print(f"Results directory: {results_dir}")
    print("="*80 + "\n")
    
    # Find trajectory files
    traj_files = list(results_dir.glob("drift_*.h5"))
    
    if not traj_files:
        print(f"❌ No trajectory files found in {results_dir}")
        return
    
    print(f"Found {len(traj_files)} trajectory files\n")
    
    results = []
    
    for traj_file in sorted(traj_files):
        try:
            with h5py.File(traj_file, 'r') as f:
                metadata = extract_metadata(f, traj_file.name)
                v_drift, v_drift_std, t_final, r2 = extract_drift_velocity(f)

            results.append({
                **metadata,
                'v_drift': v_drift,
                'v_drift_std': v_drift_std,
                't_sim': t_final,
                'r2': r2,
                'file': traj_file.name,
            })

            print(f"✓ {traj_file.name}")
            print(f"    v_drift = {v_drift:.2f} ± {v_drift_std:.2f} m/s")
            print(f"    R² = {r2:.6f}")

        except Exception as e:
            print(f"❌ Error processing {traj_file.name}: {e}")
    
    if not results:
        print("❌ No valid results")
        return
    
    print(f"\n{'='*80}")
    print(f"Successfully analyzed {len(results)} trajectories")
    print("="*80 + "\n")
    
    # Group results by species and model for plotting
    species_list = sorted({r['species'] for r in results}) or ['Unknown']
    models = sorted({r['model'] for r in results}) or ['Unknown']

    n_rows = max(1, len(species_list))
    n_cols = max(1, len(models))
    fig, axes = plt.subplots(n_rows, n_cols, figsize=(5 * n_cols, 4 * n_rows))
    fig.suptitle('Transport Physics: Drift Velocity vs E/N', fontsize=14, fontweight='bold')
    axes_array = np.array(axes).reshape(n_rows, n_cols)

    for i, species in enumerate(species_list):
        for j, model in enumerate(models):
            ax = axes_array[i, j]

            subset = [r for r in results if r['species'] == species and r['model'] == model]

            if not subset:
                ax.text(0.5, 0.5, 'No data', ha='center', va='center', transform=ax.transAxes)
                ax.set_title(f'{species} - {model}')
                ax.set_xlabel('E/N (Td)')
                ax.set_ylabel('Drift Velocity (m/s)')
                continue

            temps = sorted({r['temp_K'] for r in subset if r['temp_K'] is not None}) or [None]

            for temp in temps:
                temp_subset = [r for r in subset if r['temp_K'] == temp]
                if not temp_subset:
                    continue

                valid_subset = [r for r in temp_subset if r['en_td'] is not None]
                if not valid_subset:
                    continue

                en_vals = np.array([r['en_td'] for r in valid_subset])
                v_vals = np.array([r['v_drift'] for r in valid_subset])
                v_std = np.array([r['v_drift_std'] for r in valid_subset])

                sort_idx = np.argsort(en_vals)
                en_vals = en_vals[sort_idx]
                v_vals = v_vals[sort_idx]
                v_std = v_std[sort_idx]

                label = f'{temp:.0f} K' if temp is not None else 'Unknown T'
                ax.errorbar(en_vals, v_vals, yerr=v_std,
                            marker='o', label=label, capsize=3)

            ax.set_xlabel('E/N (Td)')
            ax.set_ylabel('Drift Velocity (m/s)')
            ax.set_title(f'{species} - {model}')
            if ax.get_legend_handles_labels()[0]:
                ax.legend()
            ax.grid(True, alpha=0.3)
    
    plt.tight_layout()
    
    # Save plot
    plot_file = results_dir / 'drift_velocity_analysis.png'
    plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    print(f"📊 Saved plot: {plot_file}")
    
    # Summary statistics
    print("\n" + "="*80)
    print("VALIDATION SUMMARY")
    print("="*80)
    
    for species in species_list:
        print(f"\n{species} (K₀_lit = {K0_LITERATURE.get(species, 'N/A')} cm²/(V·s)):")
        for model in models:
            subset = [r for r in results if r['species'] == species and r['model'] == model]
            if subset:
                r2_mean = np.mean([r['r2'] for r in subset])
                en_range = [r['en_td'] for r in subset if r['en_td'] is not None]
                en_summary = f"E/N span {min(en_range):.1f}-{max(en_range):.1f} Td" if en_range else "E/N unknown"
                print(f"  {model:10s}: {len(subset):3d} configs, R² = {r2_mean:.6f} ({en_summary})")
    
    print("\n" + "="*80 + "\n")

if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        results_dir = Path(sys.argv[1])
    else:
        results_dir = Path("../results/v1.0_test/physics/transport/drift")
    
    analyze_drift_results(results_dir)

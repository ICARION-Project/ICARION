#!/usr/bin/env python3
"""
Comprehensive Thermalization Analysis for ICARION v1.0
======================================================
Analyzes thermalization test results including:
- Temperature evolution and final accuracy
- Maxwell-Boltzmann velocity distribution fitting
- Statistical tests (Chi-squared, Kolmogorov-Smirnov)
- Visualization of velocity distributions
- Comprehensive summary statistics

Author: ICARION Validation Suite
"""

import h5py
import numpy as np
import json
from pathlib import Path
from scipy import stats
from scipy.optimize import curve_fit
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

# Physical constants
K_BOLTZMANN = 1.380649e-23  # J/K
AMU_TO_KG = 1.66053906660e-27

# Ion species data
SPECIES_DATA = {
    "H3O+": {"mass_amu": 19.0, "color": "#1f77b4"},
    "PentanalH+": {"mass_amu": 87.0, "color": "#ff7f0e"},
    "26DTBPH+": {"mass_amu": 192.0, "color": "#2ca02c"}
}

def load_hdf5_data(filepath):
    """Load trajectory data from HDF5 file"""
    with h5py.File(filepath, 'r') as f:
        velocities = f['/trajectory/velocities'][:]  # (time, ions, 3)
        times = f['/trajectory/time'][:]
        
        # Get metadata from HDF5 structure (per HDF5_OUTPUT_STRUCTURE.md)
        target_temp = float(f['/domains/domain_0/environment/temperature_K'][()])
        pressure = float(f['/domains/domain_0/environment/pressure_Pa'][()])
        
        # Get species info from /ions/initial_species_id
        species_ids_raw = f['/ions/initial_species_id'][:]
        # Decode bytes to string and get first (they should all be the same for thermalization tests)
        if isinstance(species_ids_raw[0], bytes):
            species_id = species_ids_raw[0].decode('utf-8')
        else:
            species_id = str(species_ids_raw[0])
        
        # Get collision model from /metadata/config/collision_model
        if '/metadata/config/collision_model' in f:
            collision_model_raw = f['/metadata/config/collision_model'][()]
            if isinstance(collision_model_raw, bytes):
                collision_model = collision_model_raw.decode('utf-8')
            else:
                collision_model = str(collision_model_raw)
        else:
            # Fallback: parse from filename
            filename = Path(filepath).stem
            if filename.startswith('therm_hss_'):
                collision_model = 'HSS'
            elif filename.startswith('therm_ehss_'):
                collision_model = 'EHSS'
            else:
                collision_model = 'Unknown'
        
    return {
        'velocities': velocities,
        'times': times,
        'species_id': species_id,
        'target_temp': target_temp,
        'pressure': pressure,
        'collision_model': collision_model,
        'mass_amu': SPECIES_DATA.get(species_id, {'mass_amu': 19.0})['mass_amu']
    }

def calculate_temperature(velocities, mass_kg):
    """
    Calculate temperature from velocity distribution.
    T = m * <v²> / (3 * k_B)
    """
    v_squared = np.sum(velocities**2, axis=-1)  # |v|² for each ion
    mean_v_squared = np.mean(v_squared)
    temperature = (mass_kg * mean_v_squared) / (3 * K_BOLTZMANN)
    return temperature, mean_v_squared

def maxwell_boltzmann_speed(v, T, mass_kg):
    """
    Maxwell-Boltzmann speed distribution (3D)
    P(v) = 4π * (m/(2πkT))^(3/2) * v² * exp(-mv²/(2kT))
    """
    factor = np.sqrt(2 / np.pi) * (mass_kg / (K_BOLTZMANN * T))**(3/2)
    return factor * v**2 * np.exp(-mass_kg * v**2 / (2 * K_BOLTZMANN * T))

def analyze_velocity_distribution(velocities, mass_kg, target_temp, n_bins=50):
    """
    Analyze velocity distribution and compare to Maxwell-Boltzmann.
    Returns chi-squared and KS test statistics.
    """
    # Calculate speed magnitudes
    speeds = np.sqrt(np.sum(velocities**2, axis=-1)).flatten()
    
    # Create histogram
    counts, bin_edges = np.histogram(speeds, bins=n_bins, density=True)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_width = bin_edges[1] - bin_edges[0]
    
    # Theoretical MB distribution at target temperature
    theoretical = maxwell_boltzmann_speed(bin_centers, target_temp, mass_kg)
    
    # Chi-squared test (where counts are sufficient)
    valid_bins = counts > 5  # Standard chi-squared requirement
    if np.sum(valid_bins) > 10:
        observed = counts[valid_bins] * len(speeds) * bin_width
        expected = theoretical[valid_bins] * len(speeds) * bin_width
        chi_squared = np.sum((observed - expected)**2 / expected)
        dof = np.sum(valid_bins) - 1
        p_value_chi = 1 - stats.chi2.cdf(chi_squared, dof)
    else:
        chi_squared, p_value_chi, dof = None, None, None
    
    # Kolmogorov-Smirnov test (non-parametric)
    # Generate theoretical samples from MB distribution
    def mb_cdf(v, T, m):
        """Cumulative distribution function for MB speed distribution"""
        # Using erf function for analytical CDF
        a = np.sqrt(m / (2 * K_BOLTZMANN * T))
        return stats.gamma.cdf(a * v, 3/2)
    
    # KS test
    theoretical_cdf = mb_cdf(speeds, target_temp, mass_kg)
    empirical_cdf = np.arange(1, len(speeds) + 1) / len(speeds)
    speeds_sorted = np.sort(speeds)
    
    # Simple KS statistic
    ks_statistic = np.max(np.abs(np.sort(theoretical_cdf) - empirical_cdf))
    
    # Alternative: use scipy's KS test with custom distribution
    try:
        from scipy.stats import kstest
        # Generate samples from theoretical MB
        v_most_probable = np.sqrt(2 * K_BOLTZMANN * target_temp / mass_kg)
        theoretical_samples = np.random.rayleigh(v_most_probable / np.sqrt(2), size=len(speeds))
        ks_result = stats.ks_2samp(speeds, theoretical_samples)
        ks_statistic = ks_result.statistic
        p_value_ks = ks_result.pvalue
    except:
        p_value_ks = None
    
    # Calculate RMS difference
    rms_diff = np.sqrt(np.mean((counts[valid_bins] - theoretical[valid_bins])**2))
    
    # Most probable speed comparison
    v_mp_theoretical = np.sqrt(2 * K_BOLTZMANN * target_temp / mass_kg)
    v_mp_observed = bin_centers[np.argmax(counts)]
    v_mp_error = abs(v_mp_observed - v_mp_theoretical) / v_mp_theoretical * 100
    
    # Mean speed comparison
    v_mean_theoretical = np.sqrt(8 * K_BOLTZMANN * target_temp / (np.pi * mass_kg))
    v_mean_observed = np.mean(speeds)
    v_mean_error = abs(v_mean_observed - v_mean_theoretical) / v_mean_theoretical * 100
    
    # RMS speed comparison
    v_rms_theoretical = np.sqrt(3 * K_BOLTZMANN * target_temp / mass_kg)
    v_rms_observed = np.sqrt(np.mean(speeds**2))
    v_rms_error = abs(v_rms_observed - v_rms_theoretical) / v_rms_theoretical * 100
    
    return {
        'speeds': speeds,
        'bin_centers': bin_centers,
        'counts': counts,
        'theoretical': theoretical,
        'chi_squared': chi_squared,
        'p_value_chi': p_value_chi,
        'dof': dof,
        'ks_statistic': ks_statistic,
        'p_value_ks': p_value_ks,
        'rms_diff': rms_diff,
        'v_mp_theoretical': v_mp_theoretical,
        'v_mp_observed': v_mp_observed,
        'v_mp_error': v_mp_error,
        'v_mean_theoretical': v_mean_theoretical,
        'v_mean_observed': v_mean_observed,
        'v_mean_error': v_mean_error,
        'v_rms_theoretical': v_rms_theoretical,
        'v_rms_observed': v_rms_observed,
        'v_rms_error': v_rms_error
    }

def analyze_single_test(filepath, verbose=True):
    """Analyze a single thermalization test"""
    
    # Load data
    data = load_hdf5_data(filepath)
    mass_kg = data['mass_amu'] * AMU_TO_KG
    
    # Analyze initial state
    v_initial = data['velocities'][0]
    temp_initial, v2_initial = calculate_temperature(v_initial, mass_kg)
    
    # Analyze final state
    v_final = data['velocities'][-1]
    temp_final, v2_final = calculate_temperature(v_final, mass_kg)
    
    # Temperature metrics
    temp_error_initial = abs(temp_initial - data['target_temp']) / data['target_temp'] * 100
    temp_error_final = abs(temp_final - data['target_temp']) / data['target_temp'] * 100
    temp_change = temp_final - temp_initial
    
    # Analyze velocity distribution
    dist_analysis = analyze_velocity_distribution(v_final, mass_kg, data['target_temp'])
    
    # Determine quality (scientific thresholds)
    if temp_error_final < 2.5 and dist_analysis['v_rms_error'] < 2.5:
        quality = "EXCELLENT"
        emoji = "✅"
    elif temp_error_final < 5.0 and dist_analysis['v_rms_error'] < 5.0:
        quality = "GOOD"
        emoji = "✅"
    elif temp_error_final < 10.0:
        quality = "ACCEPTABLE"
        emoji = "⚠️ "
    else:
        quality = "POOR"
        emoji = "❌"
    
    if verbose:
        print(f"\n{'='*80}")
        print(f"File: {Path(filepath).name}")
        print(f"{'='*80}")
        print(f"Species: {data['species_id']} | Model: {data['collision_model']} | "
              f"T={data['target_temp']:.0f}K | P={data['pressure']:.1f}Pa")
        print(f"\n--- Temperature Analysis ---")
        print(f"  Initial:  {temp_initial:7.1f} K  (error: {temp_error_initial:5.1f}%)")
        print(f"  Final:    {temp_final:7.1f} K  (error: {temp_error_final:5.1f}%)")
        print(f"  Target:   {data['target_temp']:7.1f} K")
        print(f"  Change:   {temp_change:+7.1f} K")
        
        print(f"\n--- Maxwell-Boltzmann Distribution Analysis ---")
        print(f"  Most Probable Speed:")
        print(f"    Theory:   {dist_analysis['v_mp_theoretical']:7.1f} m/s")
        print(f"    Observed: {dist_analysis['v_mp_observed']:7.1f} m/s")
        print(f"    Error:    {dist_analysis['v_mp_error']:7.2f}%")
        
        print(f"\n  Mean Speed:")
        print(f"    Theory:   {dist_analysis['v_mean_theoretical']:7.1f} m/s")
        print(f"    Observed: {dist_analysis['v_mean_observed']:7.1f} m/s")
        print(f"    Error:    {dist_analysis['v_mean_error']:7.2f}%")
        
        print(f"\n  RMS Speed:")
        print(f"    Theory:   {dist_analysis['v_rms_theoretical']:7.1f} m/s")
        print(f"    Observed: {dist_analysis['v_rms_observed']:7.1f} m/s")
        print(f"    Error:    {dist_analysis['v_rms_error']:7.2f}%")
        
        print(f"\n--- Statistical Tests ---")
        if dist_analysis['chi_squared'] is not None:
            print(f"  Chi-squared: χ² = {dist_analysis['chi_squared']:.2f} (dof={dist_analysis['dof']})")
            print(f"               p-value = {dist_analysis['p_value_chi']:.4f}")
        if dist_analysis['p_value_ks'] is not None:
            print(f"  KS Test:     D = {dist_analysis['ks_statistic']:.4f}")
            print(f"               p-value = {dist_analysis['p_value_ks']:.4f}")
        
        print(f"\n  Status: {emoji} {quality}")
        print(f"{'='*80}")
    
    return {
        'filepath': str(filepath),
        'species': data['species_id'],
        'collision_model': data['collision_model'],
        'target_temp': data['target_temp'],
        'pressure': data['pressure'],
        'temp_initial': temp_initial,
        'temp_final': temp_final,
        'temp_error_initial': temp_error_initial,
        'temp_error_final': temp_error_final,
        'temp_change': temp_change,
        'v_mp_error': dist_analysis['v_mp_error'],
        'v_mean_error': dist_analysis['v_mean_error'],
        'v_rms_error': dist_analysis['v_rms_error'],
        'chi_squared': dist_analysis['chi_squared'],
        'p_value_chi': dist_analysis['p_value_chi'],
        'ks_statistic': dist_analysis['ks_statistic'],
        'p_value_ks': dist_analysis['p_value_ks'],
        'quality': quality
    }

def analyze_batch(results_dir="results/v1.0_test/physics/thermalization", pattern="*.h5"):
    """Analyze all thermalization tests in batch"""
    
    results_path = Path(results_dir)
    h5_files = sorted(results_path.glob(pattern))
    
    if not h5_files:
        print(f"No HDF5 files found in {results_dir}")
        return
    
    # Prepare log file
    validation_dir = Path(__file__).parent.parent.parent
    log_dir = validation_dir / "logs"
    log_dir.mkdir(exist_ok=True)
    log_file = log_dir / "THERMALIZATION_ANALYSIS_LOG.txt"
    log_lines = []
    
    def log_print(msg):
        """Print and log simultaneously"""
        print(msg)
        log_lines.append(msg)
    
    log_print(f"\n{'='*80}")
    log_print(f"ICARION v1.0 - Comprehensive Thermalization Validation")
    log_print(f"{'='*80}")
    log_print(f"Analyzing {len(h5_files)} thermalization tests...")
    log_print(f"Location: {results_dir}")
    log_print(f"Analysis Date: {Path(__file__).stat().st_mtime}")
    log_print(f"{'='*80}\n")
    
    all_results = []
    
    for h5_file in h5_files:
        try:
            result = analyze_single_test(h5_file, verbose=False)
            all_results.append(result)
            
            # Compact output
            status_emoji = {"EXCELLENT": "✅", "GOOD": "✅", "ACCEPTABLE": "⚠️ ", "POOR": "❌"}[result['quality']]
            msg = (f"{status_emoji} {Path(h5_file).stem:45s} | "
                   f"T_err: {result['temp_error_final']:5.1f}% | "
                   f"MB_err: {result['v_rms_error']:5.1f}% | "
                   f"{result['quality']}")
            log_print(msg)
            
        except Exception as e:
            msg = f"❌ {Path(h5_file).stem:45s} | ERROR: {e}"
            log_print(msg)
    
    # Summary statistics
    if all_results:
        log_print(f"\n{'='*80}")
        log_print(f"SUMMARY STATISTICS")
        log_print(f"{'='*80}")
        
        temp_errors = [r['temp_error_final'] for r in all_results]
        v_rms_errors = [r['v_rms_error'] for r in all_results]
        v_mean_errors = [r['v_mean_error'] for r in all_results]
        v_mp_errors = [r['v_mp_error'] for r in all_results]
        
        log_print(f"\nTotal tests analyzed: {len(all_results)}")
        
        log_print(f"\n--- Temperature Accuracy ---")
        log_print(f"  Mean error:    {np.mean(temp_errors):5.2f}%")
        log_print(f"  Std deviation: {np.std(temp_errors):5.2f}%")
        log_print(f"  Min error:     {np.min(temp_errors):5.2f}%")
        log_print(f"  Max error:     {np.max(temp_errors):5.2f}%")
        
        log_print(f"\n--- Maxwell-Boltzmann Distribution Accuracy ---")
        log_print(f"  RMS Speed Error:")
        log_print(f"    Mean:        {np.mean(v_rms_errors):5.2f}%")
        log_print(f"    Std:         {np.std(v_rms_errors):5.2f}%")
        log_print(f"    Range:       {np.min(v_rms_errors):5.2f}% to {np.max(v_rms_errors):5.2f}%")
        
        log_print(f"\n  Mean Speed Error:")
        log_print(f"    Mean:        {np.mean(v_mean_errors):5.2f}%")
        log_print(f"    Std:         {np.std(v_mean_errors):5.2f}%")
        log_print(f"    Range:       {np.min(v_mean_errors):5.2f}% to {np.max(v_mean_errors):5.2f}%")
        
        log_print(f"\n  Most Probable Speed Error:")
        log_print(f"    Mean:        {np.mean(v_mp_errors):5.2f}%")
        log_print(f"    Std:         {np.std(v_mp_errors):5.2f}%")
        log_print(f"    Range:       {np.min(v_mp_errors):5.2f}% to {np.max(v_mp_errors):5.2f}%")
        
        # Quality distribution
        quality_counts = {}
        for r in all_results:
            quality_counts[r['quality']] = quality_counts.get(r['quality'], 0) + 1
        
        log_print(f"\n--- Quality Distribution ---")
        for quality in ["EXCELLENT", "GOOD", "ACCEPTABLE", "POOR"]:
            count = quality_counts.get(quality, 0)
            percent = count / len(all_results) * 100
            emoji = {"EXCELLENT": "✅", "GOOD": "✅", "ACCEPTABLE": "⚠️ ", "POOR": "❌"}[quality]
            log_print(f"  {emoji} {quality:12s}: {count:3d} ({percent:5.1f}%)")
        
        # Overall assessment
        excellent_rate = quality_counts.get("EXCELLENT", 0) / len(all_results) * 100
        good_or_better = (quality_counts.get("EXCELLENT", 0) + quality_counts.get("GOOD", 0)) / len(all_results) * 100
        acceptable_or_better = good_or_better + quality_counts.get("ACCEPTABLE", 0) / len(all_results) * 100
        
        log_print(f"\n--- Overall Assessment ---")
        if excellent_rate >= 95:
            log_print(f"OUTSTANDING: {excellent_rate:.1f}% excellent results (< 2.5% error)")
        elif excellent_rate >= 75:
            log_print(f"VERY GOOD: {excellent_rate:.1f}% excellent, {good_or_better:.1f}% good or better")
        elif good_or_better >= 90:
            log_print(f"GOOD: {good_or_better:.1f}% within 5% tolerance")
        elif acceptable_or_better >= 90:
            log_print(f"ACCEPTABLE: {acceptable_or_better:.1f}% within 10% tolerance")
        else:
            log_print(f"NEEDS IMPROVEMENT: Only {acceptable_or_better:.1f}% within acceptable range")
        
        log_print(f"\n{'='*80}")
        
        # Write log file
        with open(log_file, 'w') as f:
            f.write('\n'.join(log_lines))
        print(f"\n📄 Full analysis log saved to: {log_file.relative_to(Path.cwd())}")

def plot_thermalization_curves(results_dir="results/v1.0_test/physics/thermalization", 
                               output_dir="validation/figures"):
    """
    Generate thermalization curve plots for validation report.
    Creates two figures: one for HSS, one for EHSS showing all three ion species.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # Select representative tests: 300K, 20Pa for all three ions
    test_configs = [
        ("HSS", "H3O+", "therm_hss_H3Op_300K_20.0Pa.h5"),
        ("HSS", "PentanalH+", "therm_hss_PentanalHp_300K_20.0Pa.h5"),
        ("HSS", "26DTBPH+", "therm_hss_26DTBPHp_300K_20.0Pa.h5"),
        ("EHSS", "H3O+", "therm_ehss_H3Op_300K_20.0Pa.h5"),
        ("EHSS", "PentanalH+", "therm_ehss_PentanalHp_300K_20.0Pa.h5"),
        ("EHSS", "26DTBPH+", "therm_ehss_26DTBPHp_300K_20.0Pa.h5"),
    ]
    
    results_path = Path(results_dir)
    
    # Create figure for each collision model
    for model in ["HSS", "EHSS"]:
        fig, ax = plt.subplots(figsize=(10, 6))
        
        for cm, species, filename in test_configs:
            if cm != model:
                continue
                
            filepath = results_path / filename
            if not filepath.exists():
                print(f"Warning: File not found: {filepath}")
                continue
            
            try:
                data = load_hdf5_data(filepath)
                mass_kg = data['mass_amu'] * AMU_TO_KG
                
                # Calculate temperature at each timestep
                temps = []
                for i in range(len(data['times'])):
                    v = data['velocities'][i]
                    temp, _ = calculate_temperature(v, mass_kg)
                    temps.append(temp)
                
                # Plot
                times_ms = data['times'] * 1000  # Convert to ms
                color = SPECIES_DATA[species]['color']
                ax.plot(times_ms, temps, label=species, color=color, linewidth=2)
                
                # Target line
                if species == "H3O+":  # Only draw once
                    ax.axhline(data['target_temp'], color='black', linestyle='--', 
                              linewidth=1.5, label='Target (300 K)', alpha=0.7)
                
            except Exception as e:
                print(f"Error processing {filename}: {e}")
        
        ax.set_xlabel('Time (ms)', fontsize=12)
        ax.set_ylabel('Temperature (K)', fontsize=12)
        ax.set_title(f'Thermalization Curves - {model} Model (300 K, 20 Pa)', fontsize=14, fontweight='bold')
        ax.legend(fontsize=11, loc='best')
        ax.grid(True, alpha=0.3)
        ax.set_ylim([0, 350])
        
        plt.tight_layout()
        output_file = output_path / f'thermalization_{model.lower()}_300K_20Pa.png'
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        print(f"Saved: {output_file}")
        plt.close()

def plot_velocity_distributions(results_dir="results/v1.0_test/physics/thermalization",
                                output_dir="validation/figures"):
    """
    Generate velocity distribution plots comparing simulation to Maxwell-Boltzmann.
    Creates figure with 3 subplots (one per ion species) for EHSS at 300K, 20Pa.
    """
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)
    
    # EHSS tests at 300K, 20Pa
    test_configs = [
        ("H3O+", "therm_ehss_H3Op_300K_20.0Pa.h5"),
        ("PentanalH+", "therm_ehss_PentanalHp_300K_20.0Pa.h5"),
        ("26DTBPH+", "therm_ehss_26DTBPHp_300K_20.0Pa.h5"),
    ]
    
    results_path = Path(results_dir)
    
    # Create figure with 3 subplots
    fig = plt.figure(figsize=(15, 5))
    gs = gridspec.GridSpec(1, 3, figure=fig)
    
    for idx, (species, filename) in enumerate(test_configs):
        ax = fig.add_subplot(gs[0, idx])
        
        filepath = results_path / filename
        if not filepath.exists():
            print(f"Warning: File not found: {filepath}")
            continue
        
        try:
            data = load_hdf5_data(filepath)
            mass_kg = data['mass_amu'] * AMU_TO_KG
            
            # Get final velocities
            v_final = data['velocities'][-1]
            speeds = np.sqrt(np.sum(v_final**2, axis=-1)).flatten()
            
            # Create histogram
            n_bins = 50
            counts, bin_edges = np.histogram(speeds, bins=n_bins, density=True)
            bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
            bin_width = bin_edges[1] - bin_edges[0]
            
            # Theoretical Maxwell-Boltzmann
            v_range = np.linspace(0, np.max(speeds), 200)
            mb_theory = maxwell_boltzmann_speed(v_range, data['target_temp'], mass_kg)
            
            # Plot
            color = SPECIES_DATA[species]['color']
            ax.bar(bin_centers, counts, width=bin_width, alpha=0.6, 
                   color=color, label=f'{species} Simulation', edgecolor='black', linewidth=0.5)
            ax.plot(v_range, mb_theory, 'k-', linewidth=2.5, 
                   label='Maxwell-Boltzmann', zorder=10)
            
            # Calculate and display metrics
            temp_final, _ = calculate_temperature(v_final, mass_kg)
            temp_error = abs(temp_final - data['target_temp']) / data['target_temp'] * 100
            
            v_rms_theory = np.sqrt(3 * K_BOLTZMANN * data['target_temp'] / mass_kg)
            v_rms_obs = np.sqrt(np.mean(speeds**2))
            v_rms_error = abs(v_rms_obs - v_rms_theory) / v_rms_theory * 100
            
            # Add text box with metrics
            textstr = f'T = {temp_final:.1f} K ({temp_error:.1f}% error)\n'
            textstr += f'$v_{{rms}}$ error: {v_rms_error:.2f}%'
            props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
            ax.text(0.97, 0.97, textstr, transform=ax.transAxes, fontsize=10,
                   verticalalignment='top', horizontalalignment='right', bbox=props)
            
            ax.set_xlabel('Speed (m/s)', fontsize=11)
            ax.set_ylabel('Probability Density', fontsize=11)
            ax.set_title(f'{species}', fontsize=12, fontweight='bold')
            ax.legend(fontsize=9, loc='upper left')
            ax.grid(True, alpha=0.3)
            
        except Exception as e:
            print(f"Error processing {filename}: {e}")
            import traceback
            traceback.print_exc()
    
    fig.suptitle('Velocity Distribution - EHSS Model (300 K, 20 Pa)', 
                fontsize=14, fontweight='bold', y=1.02)
    plt.tight_layout()
    
    output_file = output_path / 'velocity_distributions_ehss_300K_20Pa.png'
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    print(f"Saved: {output_file}")
    plt.close()

def plot_error_heatmap():
    """Generate heatmap of final temperature errors across pressure/temperature space"""
    validation_dir = Path(__file__).parent.parent.parent
    # Use the actual results location
    results_dir = Path("results/v1.0_test/physics/thermalization")
    figures_dir = validation_dir / "figures"
    
    # Collect all test results
    results_by_model = {"hss": {}, "ehss": {}}
    
    for model in ["hss", "ehss"]:
        for species in SPECIES_DATA.keys():
            for temp_K in [150, 300, 1000]:
                for pressure_Pa in [0.2, 2, 20, 200, 2000]:
                    # Find matching HDF5 file - use actual naming pattern (always .0Pa format)
                    filename = f"therm_{model}_{species.replace('+', 'p')}_{temp_K}K_{pressure_Pa:.1f}Pa.h5"
                    filepath = results_dir / filename
                    
                    if not filepath.exists():
                        continue
                    
                    try:
                        data = load_hdf5_data(filepath)
                        # Calculate final temperature from velocities
                        mass_kg = data['mass_amu'] * AMU_TO_KG
                        v_final = data['velocities'][-1]  # Last timestep
                        T_final, _ = calculate_temperature(v_final, mass_kg)
                        T_target = data['target_temp']
                        error_pct = abs(T_final - T_target) / T_target * 100
                        
                        key = (temp_K, pressure_Pa)
                        if key not in results_by_model[model]:
                            results_by_model[model][key] = []
                        results_by_model[model][key].append(error_pct)
                    except Exception as e:
                        print(f"   Warning: Could not process {filepath.name}: {e}")
    
    # Create figure with two subplots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    
    temps = [150, 300, 1000]
    pressures = [0.2, 2, 20, 200, 2000]
    pressure_labels = ["0.2", "2", "20", "200", "2000"]
    
    for idx, (model, ax) in enumerate(zip(["hss", "ehss"], axes)):
        # Build error matrix (temps × pressures)
        error_matrix = np.zeros((len(temps), len(pressures)))
        
        for i, temp in enumerate(temps):
            for j, pressure in enumerate(pressures):
                key = (temp, pressure)
                if key in results_by_model[model] and results_by_model[model][key]:
                    # Average over all species at this condition
                    error_matrix[i, j] = np.mean(results_by_model[model][key])
                else:
                    error_matrix[i, j] = np.nan
        
        # Plot heatmap (vmax=2.5 matches excellent threshold)
        im = ax.imshow(error_matrix, aspect='auto', cmap='RdYlGn_r', 
                      vmin=0, vmax=2.5, interpolation='nearest')
        
        # Set ticks and labels
        ax.set_xticks(np.arange(len(pressures)))
        ax.set_yticks(np.arange(len(temps)))
        ax.set_xticklabels(pressure_labels)
        ax.set_yticklabels(temps)
        
        ax.set_xlabel('Pressure [Pa]', fontsize=11, fontweight='bold')
        ax.set_ylabel('Temperature [K]', fontsize=11, fontweight='bold')
        ax.set_title(f'{model.upper()} Model', fontsize=12, fontweight='bold')
        
        # Add text annotations
        for i in range(len(temps)):
            for j in range(len(pressures)):
                if not np.isnan(error_matrix[i, j]):
                    text = ax.text(j, i, f'{error_matrix[i, j]:.2f}%',
                                 ha="center", va="center", color="black", 
                                 fontsize=9, fontweight='bold')
        
        # Colorbar
        cbar = plt.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label('Temperature Error [%]', rotation=270, labelpad=20, 
                      fontsize=10, fontweight='bold')
    
    plt.suptitle('Final Temperature Accuracy Across Parameter Space', 
                fontsize=14, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    
    output_file = figures_dir / "temperature_error_heatmap.png"
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"   Saved: {output_file.name}")


def generate_validation_figures():
    """Generate all figures for validation report"""
    print("\n" + "="*80)
    print("GENERATING VALIDATION REPORT FIGURES")
    print("="*80)
    
    print("\n1. Thermalization Curves...")
    plot_thermalization_curves()
    
    print("\n2. Velocity Distributions...")
    plot_velocity_distributions()
    
    print("\n3. Temperature Error Heatmap...")
    plot_error_heatmap()
    
    print("\n" + "="*80)
    print("✅ All validation figures generated!")
    print("="*80 + "\n")

def main():
    if len(sys.argv) > 1:
        if sys.argv[1] == "--figures":
            # Generate validation figures
            generate_validation_figures()
        else:
            # Single file analysis
            filepath = sys.argv[1]
            if not Path(filepath).exists():
                print(f"Error: File not found: {filepath}")
                sys.exit(1)
            analyze_single_test(filepath, verbose=True)
    else:
        # Batch analysis
        analyze_batch()

if __name__ == "__main__":
    main()

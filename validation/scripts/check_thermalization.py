#!/usr/bin/env python3
"""
Quick thermalization check from HDF5 output using h5dump
Analyzes velocity distribution evolution to verify thermalization behavior
"""

import subprocess
import json
import numpy as np
import re
import sys
from pathlib import Path

def run_h5dump(file_path, dataset):
    """Run h5dump and parse output"""
    try:
        result = subprocess.run(['h5dump', '-d', dataset, file_path], 
                              capture_output=True, text=True, check=True)
        return result.stdout
    except subprocess.CalledProcessError as e:
        print(f"Error running h5dump: {e}")
        return None

def parse_velocities(h5dump_output):
    """Parse velocity data from h5dump output"""
    # Extract data section
    lines = h5dump_output.split('\n')
    data_section = False
    velocities = []
    
    for line in lines:
        if 'DATA {' in line:
            data_section = True
            continue
        elif '}' in line and data_section:
            break
        elif data_section and line.strip():
            # Parse velocity triplets
            # Format: (time,ion,dim): vx, vy, vz,
            matches = re.findall(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', line)
            if matches:
                velocities.extend([float(v) for v in matches])
    
    return np.array(velocities)

def calculate_temperature(velocities, mass_kg):
    """Calculate temperature from velocity distribution using kinetic theory"""
    # velocities shape should be (N_time, N_ions, 3)
    # For H3O+: mass = 19.0233 * 1.66054e-27 kg
    
    # Calculate kinetic energy per particle: (1/2) * m * v²
    # Average over all dimensions: <(1/2) * m * (vx² + vy² + vz²)>
    v_squared = velocities**2
    
    # Group velocities into triplets (vx, vy, vz)
    n_components = len(v_squared) // 3
    if n_components * 3 != len(v_squared):
        print(f"Warning: velocity data not divisible by 3 ({len(v_squared)} components)")
        n_components = len(v_squared) // 3
    
    # Reshape to get (n_particles, 3) for each time step
    v_squared_grouped = v_squared[:n_components*3].reshape(-1, 3)
    
    # Calculate |v|² for each particle
    v_mag_squared = np.sum(v_squared_grouped, axis=1)
    
    # Temperature from equipartition theorem: (1/2) * m * <v²> = (3/2) * k_B * T
    # Therefore: T = m * <v²> / (3 * k_B)
    k_B = 1.380649e-23  # J/K
    mean_v_squared = np.mean(v_mag_squared)
    temperature = (mass_kg * mean_v_squared) / (3 * k_B)
    
    return temperature, mean_v_squared, len(v_squared_grouped)

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 check_thermalization.py <h5_file>")
        sys.exit(1)
    
    h5_file = sys.argv[1]
    if not Path(h5_file).exists():
        print(f"File not found: {h5_file}")
        sys.exit(1)
    
    print(f"Analyzing thermalization in: {h5_file}")
    print("=" * 60)
    
    # Get metadata first
    meta_output = run_h5dump(h5_file, "/metadata/config/dt_s")
    if meta_output:
        dt_match = re.search(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', meta_output)
        if dt_match:
            dt_s = float(dt_match.group(1))
            print(f"Timestep: {dt_s:.2e} s")
    
    temp_output = run_h5dump(h5_file, "/domains/domain_0/environment/temperature_K")
    if temp_output:
        temp_match = re.search(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', temp_output)
        if temp_match:
            target_temp = float(temp_match.group(1))
            print(f"Target temperature: {target_temp:.1f} K")
    
    # Get time points
    time_output = run_h5dump(h5_file, "/trajectory/time")
    if not time_output:
        print("Failed to get time data")
        return
    
    time_matches = re.findall(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', time_output)
    times = np.array([float(t) for t in time_matches])
    print(f"Simulation time: {times[0]:.2e} to {times[-1]:.2e} s ({len(times)} steps)")
    
    # Mass of H3O+ ion
    mass_H3Op = 19.0233 * 1.66054e-27  # kg
    
    # Analyze velocities at first and last time steps
    vel_output = run_h5dump(h5_file, "/trajectory/velocities")
    if not vel_output:
        print("Failed to get velocity data")
        return
    
    # Parse all velocity data
    all_velocities = parse_velocities(vel_output)
    print(f"Total velocity components: {len(all_velocities)}")
    
    # The data structure is (time, ion, component) = (22, 10000, 3)
    n_times = len(times)
    n_ions = 10000  # Fixed from HDF5 structure
    expected_size = n_times * n_ions * 3
    
    if len(all_velocities) != expected_size:
        print(f"Warning: data size mismatch. Expected {expected_size}, got {len(all_velocities)}")
        # Try to determine actual structure
        n_ions = len(all_velocities) // (n_times * 3)
    
    print(f"Analyzing {n_ions} ions over {n_times} time steps")
    
    # Reshape data: (time, ion, component)
    try:
        velocities = all_velocities.reshape(n_times, n_ions, 3)
    except ValueError as e:
        print(f"Reshape error: {e}")
        print(f"Trying alternative: total_size={len(all_velocities)}, n_times={n_times}")
        # Fallback: truncate to exact multiple
        truncate_size = (len(all_velocities) // 3) * 3
        velocities_truncated = all_velocities[:truncate_size].reshape(-1, 3)
        n_samples = len(velocities_truncated)
        print(f"Using {n_samples} velocity samples for analysis")
        # Just analyze initial and final samples
        mid_point = n_samples // 2
        velocities = velocities_truncated.reshape(1, n_samples, 3)  # Fake time dimension
        n_times = 1
    
    # Analyze initial and final distributions
    print("\n--- Initial State (t = 0) ---")
    v_initial = velocities[0].flatten()  # All velocity components at t=0
    temp_initial, mean_v2_initial, n_particles = calculate_temperature(v_initial, mass_H3Op)
    print(f"Temperature: {temp_initial:.1f} K")
    print(f"Mean |v|²: {mean_v2_initial:.2e} m²/s²")
    print(f"RMS velocity: {np.sqrt(mean_v2_initial):.1f} m/s")
    
    print(f"\n--- Final State (t = {times[-1]:.2e} s) ---")
    v_final = velocities[-1].flatten()  # All velocity components at final time
    temp_final, mean_v2_final, _ = calculate_temperature(v_final, mass_H3Op)
    print(f"Temperature: {temp_final:.1f} K")
    print(f"Mean |v|²: {mean_v2_final:.2e} m²/s²")
    print(f"RMS velocity: {np.sqrt(mean_v2_final):.1f} m/s")
    
    # Calculate theoretical thermal velocity at target temperature
    # <v²> = 3kT/m for 3D Maxwell-Boltzmann distribution
    v_thermal_theory = np.sqrt(3 * 1.380649e-23 * target_temp / mass_H3Op)
    print(f"\n--- Theoretical (Maxwell-Boltzmann at {target_temp}K) ---")
    print(f"RMS velocity: {v_thermal_theory:.1f} m/s")
    print(f"Mean |v|²: {v_thermal_theory**2:.2e} m²/s²")
    
    # Temperature evolution analysis
    print(f"\n--- Thermalization Analysis ---")
    temp_change = temp_final - temp_initial
    temp_error_initial = abs(temp_initial - target_temp) / target_temp * 100
    temp_error_final = abs(temp_final - target_temp) / target_temp * 100
    
    print(f"Initial temp error: {temp_error_initial:.1f}%")
    print(f"Final temp error: {temp_error_final:.1f}%")
    print(f"Temperature change: {temp_change:+.1f} K")
    
    if abs(temp_error_final) < 10:
        print("✅ THERMALIZATION SUCCESS: Final temperature within 10% of target")
    else:
        print("❌ THERMALIZATION CONCERN: Final temperature deviates >10% from target")
    
    # Check if thermalization is occurring
    if abs(temp_change) > 50:
        if temp_error_final < temp_error_initial:
            print("📈 THERMALIZATION ACTIVE: Temperature approaching target")
        else:
            print("⚠️  THERMALIZATION ISSUE: Temperature moving away from target")
    else:
        print("📊 TEMPERATURE STABLE: Small change indicates equilibrium")

if __name__ == "__main__":
    main()
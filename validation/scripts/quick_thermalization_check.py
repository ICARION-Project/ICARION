#!/usr/bin/env python3
"""
Quick thermalization check - extracts first/last timestep velocities and calculates temperatures
"""

import subprocess
import re
import numpy as np
import sys
from pathlib import Path

def get_metadata_value(h5_file, dataset_path):
    """Get a single metadata value from HDF5"""
    try:
        result = subprocess.run(['h5dump', '-d', dataset_path, h5_file], 
                              capture_output=True, text=True, check=True)
        match = re.search(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', result.stdout)
        if match:
            return float(match.group(1))
    except:
        pass
    return None

def get_velocity_timestep(h5_file, time_index):
    """Get velocities for a specific time index using h5dump with hyperslabs"""
    try:
        # Use hyperslab selection to get velocities[time_index,:,:]
        result = subprocess.run([
            'h5dump', '-d', '/trajectory/velocities', 
            '-s', f'{time_index},0,0',  # start
            '-c', '1,10000,3',          # count: 1 time step, all ions, all components
            h5_file
        ], capture_output=True, text=True, check=True)
        
        # Parse velocity data
        velocities = []
        lines = result.stdout.split('\n')
        data_section = False
        
        for line in lines:
            if 'DATA {' in line:
                data_section = True
                continue
            elif '}' in line and data_section:
                break
            elif data_section and line.strip():
                # Extract numbers
                matches = re.findall(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', line)
                velocities.extend([float(v) for v in matches])
        
        return np.array(velocities)
    
    except Exception as e:
        print(f"Error extracting velocities for time {time_index}: {e}")
        return None

def calculate_temperature_from_velocities(velocities, mass_kg):
    """Calculate temperature from velocity array"""
    if len(velocities) % 3 != 0:
        print(f"Warning: velocity array size {len(velocities)} not divisible by 3")
        velocities = velocities[:len(velocities)//3*3]
    
    # Reshape to (N_particles, 3)
    vel_3d = velocities.reshape(-1, 3)
    
    # Calculate |v|² for each particle
    v_squared = np.sum(vel_3d**2, axis=1)
    mean_v_squared = np.mean(v_squared)
    
    # Temperature from equipartition: (1/2)m<v²> = (3/2)k_B T
    k_B = 1.380649e-23  # J/K
    temperature = (mass_kg * mean_v_squared) / (3 * k_B)
    
    return temperature, mean_v_squared, len(vel_3d)

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 quick_thermalization_check.py <h5_file>")
        sys.exit(1)
    
    h5_file = sys.argv[1]
    if not Path(h5_file).exists():
        print(f"File not found: {h5_file}")
        sys.exit(1)
    
    print(f"Quick Thermalization Analysis: {Path(h5_file).name}")
    print("=" * 60)
    
    # Get basic metadata
    target_temp = get_metadata_value(h5_file, "/domains/domain_0/environment/temperature_K")
    dt_s = get_metadata_value(h5_file, "/metadata/config/dt_s")
    
    if target_temp:
        print(f"Target temperature: {target_temp:.1f} K")
    if dt_s:
        print(f"Timestep: {dt_s:.2e} s")
    
    # Get time array to know how many timesteps
    try:
        result = subprocess.run(['h5dump', '-d', '/trajectory/time', h5_file], 
                              capture_output=True, text=True, check=True)
        time_matches = re.findall(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', result.stdout)
        times = np.array([float(t) for t in time_matches])
        n_times = len(times)
        print(f"Simulation: {times[0]:.2e} to {times[-1]:.2e} s ({n_times} steps)")
    except:
        print("Could not read time data")
        return
    
    # Mass of H3O+ ion
    mass_H3Op = 19.0233 * 1.66054e-27  # kg
    
    print(f"\n--- Analysis ---")
    
    # Initial velocities (t=0, index=0)
    print("Extracting initial velocities...")
    v_initial = get_velocity_timestep(h5_file, 0)
    if v_initial is not None:
        temp_initial, v2_initial, n_particles = calculate_temperature_from_velocities(v_initial, mass_H3Op)
        print(f"Initial (t=0): {temp_initial:.1f} K, RMS={np.sqrt(v2_initial):.1f} m/s ({n_particles} ions)")
    
    # Final velocities (last timestep)
    print("Extracting final velocities...")
    v_final = get_velocity_timestep(h5_file, n_times-1)
    if v_final is not None:
        temp_final, v2_final, _ = calculate_temperature_from_velocities(v_final, mass_H3Op)
        print(f"Final (t={times[-1]:.2e}): {temp_final:.1f} K, RMS={np.sqrt(v2_final):.1f} m/s")
    
    # Theoretical values
    if target_temp:
        v_thermal_theory = np.sqrt(3 * 1.380649e-23 * target_temp / mass_H3Op)
        print(f"Theory ({target_temp}K): {target_temp:.1f} K, RMS={v_thermal_theory:.1f} m/s")
        
        # Analysis
        if v_initial is not None and v_final is not None:
            temp_change = temp_final - temp_initial
            error_initial = abs(temp_initial - target_temp) / target_temp * 100
            error_final = abs(temp_final - target_temp) / target_temp * 100
            
            print(f"\n--- Thermalization Summary ---")
            print(f"Initial error: {error_initial:.1f}%")
            print(f"Final error:   {error_final:.1f}%")
            print(f"Change:        {temp_change:+.1f} K")
            
            if error_final < 15:
                print("✅ SUCCESS: Temperature within 15% of target")
            elif error_final < error_initial:
                print("📈 IMPROVING: Temperature approaching target")
            else:
                print("⚠️  ISSUE: Temperature not thermalizing correctly")

if __name__ == "__main__":
    main()
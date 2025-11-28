#!/usr/bin/env python3
"""
Simple thermalization check using velocity samples
"""

import subprocess
import re
import numpy as np
import sys

def extract_velocity_sample(h5_file, time_idx, n_sample=1000):
    """Extract sample of velocities for temperature calculation"""
    try:
        result = subprocess.run([
            'h5dump', '-d', '/trajectory/velocities',
            '-s', f'{time_idx},0,0',
            '-c', f'1,{n_sample},3',
            h5_file
        ], capture_output=True, text=True, check=True)
        
        # Parse velocities
        velocities = []
        lines = result.stdout.split('\n')
        in_data = False
        
        for line in lines:
            if 'DATA {' in line:
                in_data = True
                continue
            elif '}' in line and in_data:
                break
            elif in_data and line.strip():
                numbers = re.findall(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', line)
                velocities.extend([float(x) for x in numbers])
        
        return np.array(velocities)
    except Exception as e:
        print(f"Error: {e}")
        return None

def calc_temp(velocities, mass_kg):
    """Calculate temperature from velocity array"""
    if len(velocities) % 3 != 0:
        velocities = velocities[:len(velocities)//3*3]
    
    vel_3d = velocities.reshape(-1, 3)
    v_sq = np.sum(vel_3d**2, axis=1)
    mean_v_sq = np.mean(v_sq)
    
    k_B = 1.380649e-23
    T = (mass_kg * mean_v_sq) / (3 * k_B)
    return T, np.sqrt(mean_v_sq)

def main():
    h5_file = sys.argv[1]
    
    # H3O+ mass
    mass = 19.0233 * 1.66054e-27
    
    print("Thermalization Check")
    print("="*40)
    
    # Get metadata
    result = subprocess.run(['h5dump', '-d', '/domains/domain_0/environment/temperature_K', h5_file], 
                          capture_output=True, text=True)
    target_temp = float(re.search(r'([+-]?\d+\.?\d*)', result.stdout).group(1))
    
    result = subprocess.run(['h5dump', '-d', '/trajectory/time', h5_file], 
                          capture_output=True, text=True)
    times = [float(x) for x in re.findall(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', result.stdout)]
    
    print(f"Target: {target_temp} K")
    print(f"Time: {times[0]:.2e} to {times[-1]:.2e} s ({len(times)} steps)")
    
    # Sample velocities
    print("\nExtracting velocity samples...")
    v_initial = extract_velocity_sample(h5_file, 0)
    v_final = extract_velocity_sample(h5_file, len(times)-1)
    
    if v_initial is not None:
        T_init, rms_init = calc_temp(v_initial, mass)
        print(f"Initial: {T_init:.1f} K, RMS={rms_init:.1f} m/s")
    
    if v_final is not None:
        T_final, rms_final = calc_temp(v_final, mass)
        print(f"Final:   {T_final:.1f} K, RMS={rms_final:.1f} m/s")
    
    # Theoretical
    v_theory = np.sqrt(3 * 1.380649e-23 * target_temp / mass)
    print(f"Theory:  {target_temp:.1f} K, RMS={v_theory:.1f} m/s")
    
    if v_initial is not None and v_final is not None:
        error_init = abs(T_init - target_temp) / target_temp * 100
        error_final = abs(T_final - target_temp) / target_temp * 100
        
        print(f"\nError: {error_init:.1f}% → {error_final:.1f}%")
        
        if error_final < 20:
            print("✅ Good thermalization")
        elif error_final < error_init:
            print("📈 Thermalizing...")
        else:
            print("⚠️ Poor thermalization")

if __name__ == "__main__":
    main()
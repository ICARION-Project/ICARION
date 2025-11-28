#!/usr/bin/env python3

import subprocess
import re
import numpy as np

h5_file = "results/v1.0_test/physics/thermalization/therm_ehss_H3Op_300K_20.0Pa.h5"

print("=== Quick Thermalization Check ===")

# Get target temperature 
result = subprocess.run(['h5dump', '-d', '/domains/domain_0/environment/temperature_K', h5_file], 
                      capture_output=True, text=True)
temp_line = [line for line in result.stdout.split('\n') if '(' in line and '):' in line][0]
target_temp = float(re.search(r':\s*([0-9.]+)', temp_line).group(1))

print(f"Target Temperature: {target_temp} K")

# Sample initial velocities (first 100 ions)
print("\nSampling initial velocities...")
result = subprocess.run([
    'h5dump', '-d', '/trajectory/velocities',
    '-s', '0,0,0', '-c', '1,100,3', h5_file
], capture_output=True, text=True)

# Extract numbers from initial velocities
v_init = []
for line in result.stdout.split('\n'):
    if line.strip() and ('(' in line and '):' in line):
        numbers = re.findall(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', line.split('):')[1])
        v_init.extend([float(x) for x in numbers])

v_init = np.array(v_init)
print(f"Initial velocities extracted: {len(v_init)} components")

# Sample final velocities (last timestep = 21, first 100 ions)  
print("Sampling final velocities...")
result = subprocess.run([
    'h5dump', '-d', '/trajectory/velocities',
    '-s', '21,0,0', '-c', '1,100,3', h5_file
], capture_output=True, text=True)

v_final = []
for line in result.stdout.split('\n'):
    if line.strip() and ('(' in line and '):' in line):
        numbers = re.findall(r'([+-]?\d+\.?\d*(?:[eE][+-]?\d+)?)', line.split('):')[1])
        v_final.extend([float(x) for x in numbers])

v_final = np.array(v_final)
print(f"Final velocities extracted: {len(v_final)} components")

# Calculate temperatures
def calc_temp(velocities, mass_kg):
    if len(velocities) % 3 != 0:
        velocities = velocities[:len(velocities)//3*3]
    
    vel_3d = velocities.reshape(-1, 3)
    v_squared = np.sum(vel_3d**2, axis=1)
    mean_v_sq = np.mean(v_squared)
    
    k_B = 1.380649e-23  # J/K
    T = (mass_kg * mean_v_sq) / (3 * k_B)
    return T, np.sqrt(mean_v_sq)

# H3O+ mass
mass_H3O = 19.0233 * 1.66054e-27  # kg

T_init, rms_init = calc_temp(v_init, mass_H3O)
T_final, rms_final = calc_temp(v_final, mass_H3O)

# Theoretical thermal velocity
v_theory = np.sqrt(3 * 1.380649e-23 * target_temp / mass_H3O)

print(f"\n=== Results ===")
print(f"Initial:  {T_init:.1f} K  (RMS = {rms_init:.1f} m/s)")
print(f"Final:    {T_final:.1f} K  (RMS = {rms_final:.1f} m/s)")
print(f"Target:   {target_temp:.1f} K  (RMS = {v_theory:.1f} m/s)")

error_init = abs(T_init - target_temp) / target_temp * 100
error_final = abs(T_final - target_temp) / target_temp * 100
temp_change = T_final - T_init

print(f"\nInitial error: {error_init:.1f}%")
print(f"Final error:   {error_final:.1f}%")
print(f"Change:        {temp_change:+.1f} K")

if error_final < 15:
    print("✅ GOOD: Final temperature within 15% of target")
elif error_final < error_init:
    print("📈 IMPROVING: Temperature approaching target")
else:
    print("⚠️ ISSUE: Temperature not improving")
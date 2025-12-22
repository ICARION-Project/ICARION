#!/usr/bin/env python3
"""
Plot T_eff vs E/N with effective temperature correction
"""

import numpy as np
import matplotlib.pyplot as plt

# Physical constants
kB = 1.380649e-23  # J/K - Boltzmann constant
amu = 1.66053906660e-27  # kg - atomic mass unit

# Parameters
K0 = 24.1e-4  # m²/(V·s) - H3O+ in He reduced mobility
M_He = 4.003 * amu  # kg - He mass
T = 300.0  # K - gas temperature
LOSCHMIDT = 2.6867774e25  # m^-3 - number density at STP

# E/N range from 1 to 100 Td
E_N_Td = np.logspace(0, 2, 100)  # 1 to 100 Td logarithmically spaced
E_N = E_N_Td * 1e-21  # Convert to V·m²

# Calculate for different pressures
pressures_Pa = [10, 100, 1000, 10000]
colors = ['red', 'orange', 'green', 'blue']

plt.figure(figsize=(12, 8))

for P_Pa, color in zip(pressures_Pa, colors):
    # Calculate number density from ideal gas law
    N_m3 = P_Pa / (kB * T)
    
    # Electric field for each E/N at this pressure
    E_Vm = E_N * N_m3
    
    # Drift velocity at low field
    N_ratio = LOSCHMIDT / N_m3
    v_drift_lowfield = K0 * E_Vm * N_ratio
    
    # Effective temperature heating
    T_heating = (M_He / (3 * kB)) * v_drift_lowfield**2
    T_eff = T + T_heating
    
    plt.loglog(E_N_Td, T_eff, label=f'{P_Pa} Pa', color=color, linewidth=2)

# Add T reference line
plt.axhline(y=T, color='black', linestyle='--', linewidth=1, label='T = 300 K')

plt.xlabel('E/N (Td)', fontsize=14)
plt.ylabel('T_eff (K)', fontsize=14)
plt.title('Effective Temperature vs E/N\n(H3O+ in He, T=300K)', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, which='both', alpha=0.3)
plt.xlim(1, 100)

# Add text box with formula
textstr = r'$T_{eff} = T + \frac{M_{He}}{3k_B}\left(K_0 \frac{E \cdot N_0}{N}\right)^2$'
plt.text(0.02, 0.98, textstr, transform=plt.gca().transAxes,
         fontsize=11, verticalalignment='top',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

# Mark our test points
test_points = [
    (10, 1000, 'green'),
    (40, 1000, 'green'),
    (100, 1000, 'green'),
]
for en, p, col in test_points:
    N = p / (kB * T)
    E = en * 1e-21 * N
    v = K0 * E * (LOSCHMIDT / N)
    T_h = (M_He / (3 * kB)) * v**2
    T_e = T + T_h
    plt.plot(en, T_e, 'o', color=col, markersize=10, markeredgecolor='black', markeredgewidth=2)

plt.tight_layout()
plt.savefig('results/v1.0_test/instruments/ims/Teff_vs_EN.png', dpi=300)
print("Plot saved to: results/v1.0_test/instruments/ims/Teff_vs_EN.png")

# Print values at our test points
print("\n" + "="*60)
print("T_eff values at test conditions:")
print("="*60)
for en, p, _ in test_points:
    N = p / (kB * T)
    E = en * 1e-21 * N
    v = K0 * E * (LOSCHMIDT / N)
    T_h = (M_He / (3 * kB)) * v**2
    T_e = T + T_h
    print(f"E/N = {en:3.0f} Td, P = {p:4.0f} Pa: T_eff = {T_e:7.1f} K (heating: {T_h:6.1f} K)")

plt.show()

#!/usr/bin/env python3
"""
Plot expected mobility vs E/N with effective temperature correction
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
T0 = 273.15  # K - reference temperature
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
    
    # Mobility with T_eff correction: K(T_eff) = K0 * sqrt(T0/T_eff)
    K_eff = K0 * np.sqrt(T0 / T_eff)
    
    # Convert to cm²/(V·s) for plotting
    K_eff_cm2 = K_eff * 1e4
    
    plt.loglog(E_N_Td, K_eff_cm2, label=f'{P_Pa} Pa', color=color, linewidth=2)

# Add K0 reference line
plt.axhline(y=K0*1e4, color='black', linestyle='--', linewidth=1, label='K₀ (low field)')

plt.xlabel('E/N (Td)', fontsize=14)
plt.ylabel('Mobility K (cm²/(V·s))', fontsize=14)
plt.title('Expected Mobility vs E/N with Effective Temperature Correction\n(H3O+ in He, T=300K)', fontsize=16)
plt.legend(fontsize=12)
plt.grid(True, which='both', alpha=0.3)
plt.xlim(1, 100)

# Add text box with formula
textstr = r'$T_{eff} = T + \frac{M_{He}}{3k_B}\left(K_0 \frac{E \cdot N_0}{N}\right)^2$' + '\n' + \
          r'$K(E/N) = K_0 \sqrt{\frac{T_0}{T_{eff}}}$'
plt.text(0.02, 0.02, textstr, transform=plt.gca().transAxes,
         fontsize=11, verticalalignment='bottom',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

plt.tight_layout()
plt.savefig('results/v1.0_test/instruments/ims/mobility_vs_EN.png', dpi=300)
print("Plot saved to: results/v1.0_test/instruments/ims/mobility_vs_EN.png")
plt.show()

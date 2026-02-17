// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#ifndef CONSTANTS_H
#define CONSTANTS_H

// Define M_PI if it's not already defined.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Fundamental physical constants:
inline constexpr double AMU_TO_KG          = 1.66053906660e-27;  // 1 atomic mass unit in kg
inline constexpr double ELEM_CHARGE_C      = 1.602176634e-19;    // Elementary charge (Coulombs)
inline constexpr double BOLTZMANN_CONSTANT = 1.380649e-23;       // Boltzmann constant (J/K)
inline constexpr double LOSCHMIDT_CONSTANT = 2.6867811e25;       // Loschmidt constant (m⁻³ at STP)
inline constexpr double EPSILON_0      = 8.854187817e-12;    // Permittivity of free space (F/m)
inline constexpr double COULOMB_CONST  = 8.987551787e9;      // Coulomb constant k_e = 1/(4πε₀) (N·m²/C²)

// Unit conversion constants:
inline constexpr double CM2_TO_M2          = 1e-4;               // cm² to m² conversion
inline constexpr double ANGSTROM_TO_M      = 1e-10;              // Å (Angstrom) to m
inline constexpr double ANGSTROM2_TO_M2    = 1e-20;              // Ų (Angstrom²) to m²
inline constexpr double ANGSTROM3_TO_M3    = 1e-30;              // ų (Angstrom³) to m³

// STP values
inline constexpr double STP_TEMP     = 273.15;  // STP temperature in K
inline constexpr double STP_PRESSURE = 101325;  // STP pressure in Pa

// Particle-specific constants:
inline constexpr double MOLAR_MASS_HE_KG = 4.002602 * AMU_TO_KG;
inline constexpr double MOLAR_MASS_AR_KG = 39.948 * AMU_TO_KG;
inline constexpr double MOLAR_MASS_CO2_KG = 44.0095 * AMU_TO_KG;
inline constexpr double MOLAR_MASS_NE_KG = 20.1797 * AMU_TO_KG;
inline constexpr double MOLAR_MASS_N2_KG = 28.0134 * AMU_TO_KG;
inline constexpr double MOLAR_MASS_O2_KG = 31.9988 * AMU_TO_KG;
inline constexpr double MOLAR_MASS_H2O_KG = 18.01528 * AMU_TO_KG;

// Polarizabilities (approximate values in cubic meters)
inline constexpr double POLARIZABILITY_HE_SI = 0.195e-30;  //DOI: 10.1088/0953-4075/43/20/202001
inline constexpr double POLARIZABILITY_AR_SI = 1.596e-30;  //DOI: 10.1088/0953-4075/43/20/202001
inline constexpr double POLARIZABILITY_CO2_SI = 2.612e-30; //DOI: 10.1063/1.45874
inline constexpr double POLARIZABILITY_NE_SI = 0.352e-30;  //DOI: 10.1088/0953-4075/43/20/202001
inline constexpr double POLARIZABILITY_N2_SI = 1.74e-30;   //DOI: 10.1063/1.431821
inline constexpr double POLARIZABILITY_O2_SI = 1.567e-30;  //DOI: 10.1098/rspa.1966.0244
inline constexpr double POLARIZABILITY_H2O_SI = 1.45e-30;  // Approximate

// Hard-sphere radii for EHSS collision model (meters)
inline constexpr double RADIUS_HE_M = 1.3e-10;   //ISBN 0470029048
inline constexpr double RADIUS_AR_M = 1.70e-10;  //ISBN 0471099856
inline constexpr double RADIUS_CO2_M = 1.65e-10; //ISBN 3319010956
inline constexpr double RADIUS_NE_M = 1.38e-10;  //ISBN 0471099856
inline constexpr double RADIUS_N2_M = 1.82e-10;  //ISBN 3319010956
inline constexpr double RADIUS_O2_M = 1.73e-10;  //ISBN 3319010956
inline constexpr double RADIUS_H2O_M = 1.58e-10; // Approximate

// Define instrument constants
inline constexpr double EJECTION_SLIT_LENGTH_Z = 20e-3;
inline constexpr double EJECTION_SLIT_LENGTH_Y = 2e-3;

#endif  // CONSTANTS_H

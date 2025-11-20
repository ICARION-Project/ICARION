// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef CONSTANTS_H
#define CONSTANTS_H

// Define M_PI if it's not already defined.
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

// Fundamental physical constants:
const double AMU_TO_KG          = 1.66053906660e-27;  // 1 atomic mass unit in kg
const double ELECTRON_CHARGE    = 1.602176634e-19;    // Elementary charge (Coulombs)
const double BOLTZMANN_CONSTANT = 1.380649e-23;       // Boltzmann constant (J/K)
const double LOSCHMIDT_CONSTANT = 2.6867811e25;       // Loschmidt constant (m⁻³ at STP)
constexpr double EPSILON_0          = 8.854187817e-12;    // Perimittivität (F/m)

// STP values
const double STP_TEMP     = 273.15;  // STP temperature in K
const double STP_PRESSURE = 101325;  // STP pressure in Pa

// Particle-specific constants:
const double MOLAR_MASS_HE_KG = 4.002602 * AMU_TO_KG;
const double MOLAR_MASS_AR_KG = 39.948 * AMU_TO_KG;
const double MOLAR_MASS_CO2_KG = 44.0095 * AMU_TO_KG;
const double MOLAR_MASS_NE_KG = 20.1797 * AMU_TO_KG;
const double MOLAR_MASS_N2_KG = 28.0134 * AMU_TO_KG;
const double MOLAR_MASS_O2_KG = 31.9988 * AMU_TO_KG;

// Polarizabilities (approximate values in cubic meters)
const double POLARIZABILITY_HE_SI = 0.195e-30;  //DOI: 10.1088/0953-4075/43/20/202001
const double POLARIZABILITY_AR_SI = 1.596e-30;  //DOI: 10.1088/0953-4075/43/20/202001
const double POLARIZABILITY_CO2_SI = 2.612e-30; //DOI: 10.1063/1.45874
const double POLARIZABILITY_NE_SI = 0.352e-30;  //DOI: 10.1088/0953-4075/43/20/202001
const double POLARIZABILITY_N2_SI = 1.74e-30;   //DOI: 10.1063/1.431821
const double POLARIZABILITY_O2_SI = 1.567e-30;  //DOI: 10.1098/rspa.1966.0244

// Hard-sphere radii for EHSS collision model (meters)
const double RADIUS_HE_M = 1.3e-10;   //ISBN 0470029048
const double RADIUS_AR_M = 1.70e-10;  //ISBN 0471099856
const double RADIUS_CO2_M = 1.65e-10; //ISBN 3319010956
const double RADIUS_NE_M = 1.38e-10;  //ISBN 0471099856
const double RADIUS_N2_M = 1.82e-10;  //ISBN 3319010956
const double RADIUS_O2_M = 1.73e-10;  //ISBN 3319010956

// Define instrument constants
const double EJECTION_SLIT_LENGTH_Z = 20e-3;
const double EJECTION_SLIT_LENGTH_Y = 2e-3;

#endif  // CONSTANTS_H
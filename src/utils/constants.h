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
// The mass of a Helium atom evaluated from its atomic mass (4.002602 amu).
const double MOLAR_MASS_HE_KG = 4.002602 * AMU_TO_KG;
// The mass of a Nitrogen molecule (N₂), approximately 28.0134 amu.
const double MOLAR_MASS_N2_KG = 28.0134 * AMU_TO_KG;

// Polarizabilities (approximate values in SI units, i.e. cubic meters)
// The polarizability of Helium (~0.205 × 10⁻³⁰ m³).
const double POLARIZABILITY_HE_SI = 0.205e-30;
// The polarizability of Nitrogen (~1.74 × 10⁻³⁰ m³).
const double POLARIZABILITY_N2_SI = 1.74e-30;

// Hard-sphere radii for EHSS collision model (meters)
// Helium radius: LJ sigma = 2.556 Å, radius = sigma/2 ≈ 1.278 Å = 1.278e-10 m
const double RADIUS_HE_M = 1.278e-10;
// Nitrogen (N2) radius: LJ sigma = 3.64 Å, radius = sigma/2 ≈ 1.82 Å = 1.82e-10 m
const double RADIUS_N2_M = 1.82e-10;

// Define instrument constants
const double EJECTION_SLIT_LENGTH_Z = 20e-3;
const double EJECTION_SLIT_LENGTH_Y = 2e-3;

#endif  // CONSTANTS_H
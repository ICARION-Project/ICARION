// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   -------------------------------------
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        defineFields.h
 *   @brief       Calculates the electric fields.
 *
 *   @details
 *   Includes quadrupolar RF fields, axial DC fields, AC fields in
 *   arbitrary directions, and Orbitrap-like potentials.
 *
 *
 *   @date        2025-10-08
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#pragma once
#include "core/types/IonState.h"
#include "core/physics/collisions/collisionHelpers.h"
#include "core/param/paramUtils.h"

// -----------------------------
// Primitive electric field functions
// -----------------------------

/**
 * @brief Quadrupolar RF + DC field for linear ion traps
 * 
 * @param ion Ion state (position used for field calculation)
 * @param RF_voltage Peak RF voltage [V]
 * @param DC_voltage DC offset voltage [V]
 * @param omega Angular frequency [rad/s] = 2π × frequency
 * @param radius Trap radius (r₀) [m]
 * @param t Current time [s]
 * 
 * @return Electric field vector [V/m]
 * 
 * Quadrupolar field: E_x = (U_DC + U_RF·cos(ωt)) · x/r₀²
 *                     E_y = -(U_DC + U_RF·cos(ωt)) · y/r₀²
 * Used for LQIT, quadrupole ion traps, and similar devices.
 */
Vec3 RFField(const IonState& ion, double RF_voltage, double DC_voltage, double omega, double radius, double t);

/**
 * @brief Uniform axial DC field
 * 
 * @param ion Ion state (z position used)
 * @param voltage Voltage difference [V] across field region
 * @param length Field region length [m]
 * 
 * @return Electric field vector [V/m] in z-direction
 * 
 * Creates constant field E_z = -voltage/length.
 * Used for drift tubes, ion mobility spectrometry, and axial trapping.
 */
Vec3 DCField(const IonState& ion, double voltage, double length);

/**
 * @brief AC field in arbitrary direction
 * 
 * @param ion Ion state (position used for field calculation)
 * @param voltage AC voltage amplitude [V]
 * @param omega Angular frequency [rad/s]
 * @param radius Characteristic length scale [m]
 * @param t Current time [s]
 * @param direction Field direction unit vector
 * 
 * @return Electric field vector [V/m]
 * 
 * Oscillating field: E = (voltage/radius) · cos(ωt) · direction
 * Used for resonant excitation and ion manipulation.
 */
Vec3 ACField(const IonState& ion, double voltage, double omega, double radius, double t,
                    const Vec3& direction);

// -----------------------------
// Special field functions
// -----------------------------

/**
 * @brief Orbitrap hyperlogarithmic field
 * 
 * @param ion Ion state
 * @param k Orbitrap field curvature parameter [V/m²]
 * @param r_char Characteristic radius [m]
 * @param length_m Axial trap length [m]
 * 
 * @return Electric field vector [V/m]
 * 
 * Orbitrap potential: U(r,z) = k/2 · (z² - r²/2 + r_char²·ln(r/r_char))
 * Creates axial harmonic oscillation with radial confinement.
 * Used in Orbitrap mass analyzers for high-resolution mass spectrometry.
 */
Vec3 OrbitrapField(const IonState& ion, double k, double r_char, double length_m);

/**
 * @brief FTICR (Fourier Transform Ion Cyclotron Resonance) Penning trap field
 * 
 * @param ion Ion state
 * @param voltage Trap voltage [V]
 * @param characteristic_length Trap cell size [m]
 * @param instrument_length_axial_m Axial trap dimension [m]
 * 
 * @return Electric field vector [V/m]
 * 
 * Axial harmonic potential for FTICR cells (cylindrical or cubic).
 * Radial confinement from magnetic field (not included here).
 */
Vec3 FTICRField(const IonState& ion, double voltage, double characteristic_length, double instrument_length_axial_m);

// -----------------------------
// Magnetic field function
// -----------------------------

/**
 * @brief Magnetic field vector at ion position
 * 
 * @param ion Ion state
 * @param dom Instrument domain (contains magnetic field configuration)
 * 
 * @return Magnetic field vector [T]
 * 
 * Currently supports uniform magnetic fields specified in domain configuration.
 * Used for FTICR, ion cyclotron resonance, and Lorentz force calculations.
 */
Vec3 MagneticFieldVec(const IonState& ion, const InstrumentDomain& dom);

/**
 * @brief Total electric field at ion position
 * 
 * @param ion Ion state
 * @param dom Instrument domain (field configuration, voltages)
 * @param t_global Current simulation time [s]
 * 
 * @return Electric field vector [V/m]
 * 
 * Combines all electric field contributions:
 * - RF quadrupole fields (if enabled)
 * - DC axial fields
 * - AC excitation fields
 * - Special fields (Orbitrap, FTICR)
 * - External field arrays (if provided)
 * 
 * This is the main field evaluation function called by the integrator.
 */
Vec3 ElectricFieldVec(const IonState& ion, const InstrumentDomain& dom, double t_global);
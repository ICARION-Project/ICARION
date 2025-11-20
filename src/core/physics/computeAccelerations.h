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
 *   @file        computeAccelerations.h
 *   @brief       Computes ion accelerations and updated velocities.
 *
 *   @details
 *   Provides the function to calculate the time derivatives of an ion’s state
 *   vector for use in numerical integration (e.g., Runge-Kutta).
 *
 *   Acceleration contributions include:
 *   - Electric fields (DC, RF, AC, Orbitrap)
 *   - Collision damping (HardSphere, Langevin, Friction, EHSS, HSMC)
 *   - Background gas flow (adds to velocity)
 *
 *   Supports multiple instruments: LQIT, IMS, SIFDT-MS, Orbitrap, QuadrupoleRF, TOF.
 *
 *   @note
 *   Returned IonState contains updated velocity and acceleration.
 *   All units are SI; fields in V/m. Collision events and boundary checks
 *   are handled elsewhere.
 *
 *
 *   @date        2025-10-06
 *   @version     0.1
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */
#pragma once

#include <vector>

#include "integrator/integrator.h"

#include "utils/constants.h"
#include "core/types/IonState.h"
#include "core/physics/fields/defineFields.h"
#include "core/physics/collisions/defineCollisionForces.h"
#include "fieldsolver/utils/IFieldProvider.h"

namespace ICARION {
namespace physics {

// -----------------------------
// ODE function
// -----------------------------

/**
 * @brief Compute ion acceleration from electric fields and collision forces
 * 
 * @param t Current simulation time [s]
 * @param ion Ion state (position, velocity, mass, charge)
 * @param gParams Global simulation parameters (field configuration, gas properties)
 * @param dom Instrument domain (boundary conditions, field/collision settings)
 * @param field_provider Field data provider (interpolates E-field at ion position)
 * @param current_ions All ion states at current time (for space-charge calculations)
 * 
 * @return IonState with updated velocity (derivative of position) and acceleration
 * 
 * Evaluates right-hand side of Newton's equations: F = ma
 * - Electric force from external fields (DC, RF, space-charge)
 * - Collision damping forces (friction, Langevin, EHSS, etc.)
 * - Optional gas flow drag
 * 
 * All units are SI (positions in m, velocities in m/s, fields in V/m).
 * Collision events and boundary checks handled by integrator.
 */
IonState compute_accelerations(double t, const IonState& ion, const GlobalParams& gParams, const InstrumentDomain& dom,
                               const IFieldProvider* field_provider, const std::vector<IonState>& current_ions);

}  // namespace physics
}  // namespace ICARION

// Bring into global namespace for backward compatibility
using ICARION::physics::compute_accelerations;

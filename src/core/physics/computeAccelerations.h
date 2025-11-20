// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   Ion Collision And Reaction IntegratiON (ICARION)
 *   ------------------------------------------------
 *   Modular framework for simulating ion trajectories in custom
 *   electric fields and background gas environments.
 *
 *   @file       computeAccelerations.h
 *   @brief      Computes ion accelerations and updated velocities
 *
 *   @details
 *   Defines functions to compute ion accelerations based on summed forces
 *  from electric fields, space charge, and background gas interactions.
 *
 *   @date       2025-11-20
 *   @version    1.0.0
 *   @authors    ICARION Development Team
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

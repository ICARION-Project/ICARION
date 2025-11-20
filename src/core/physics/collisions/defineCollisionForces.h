// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * =====================================================================
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A modular C++ framework for simulating ion trajectories 
 *   in user-defined electric fields and background gas environments.
 *
 *   @file        defineCollisionForces.h
 *   @brief       Computes collisional damping forces acting on ions.
 *
 *   @details
 *   Implements multiple collision models including:
 *    - Friction (mobility-based damping)
 *    - Langevin (ion-neutral long-range interactions)
 *    - Hard-sphere collisions
 *    - EHSS / HSMC stochastic models (no explicit damping)
 *
 *   Each function returns a force vector divided by the ion mass (acceleration).
 *
 *   @date        2025-11-10
 *   @version     1.0.0
 *
 * =====================================================================
 */
#pragma once
#include "core/types/IonState.h"
#include "core/param/paramUtils.h"

namespace ICARION {
namespace physics {

// -----------------------------
// General collision force
// -----------------------------
Vec3 CollisionForce(const IonState& ion, const GlobalParams& gParams, const InstrumentDomain& dom);

// -----------------------------
// Specific collision models
// -----------------------------
Vec3 HardSphereCollision(const IonState& ion, const InstrumentDomain& dom);
Vec3 LangevinCollision(const IonState& ion, const InstrumentDomain& dom);
Vec3 FrictionCollision(const IonState& ion, const InstrumentDomain& dom);
Vec3 EHSSCollision(const IonState& ion, const InstrumentDomain& dom);
Vec3 HSMCCollision(const IonState& ion, const InstrumentDomain& dom);

}  // namespace physics
}  // namespace ICARION

// Bring into global namespace for backward compatibility
using namespace ICARION::physics;
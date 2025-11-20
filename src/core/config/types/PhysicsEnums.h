// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_PHYSICS_ENUMS_H
#define ICARION_CONFIG_PHYSICS_ENUMS_H

namespace ICARION::config {

/**
 * @brief Collision model types
 * 
 * Defines available collision models for ion-neutral interactions.
 */
enum class CollisionModel {
    NoCollisions,       ///< No collisions (vacuum)
    HardSphere,         ///< Hard-sphere collision model
    Langevin,           ///< Langevin collision model (polarization-limited)
    Friction,           ///< Friction force model
    EHSS,               ///< Energy-dependent hard-sphere scattering
    HSMC,               ///< Hard-sphere Monte Carlo
    UnknownCollisionModel  ///< Invalid/unknown model
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_PHYSICS_ENUMS_H

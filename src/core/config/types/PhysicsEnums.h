// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
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
    HSD,                ///< Hard-Sphere Deterministic
    Langevin,           ///< Langevin collision model (polarization-limited)
    Friction,           ///< Friction force model
    EHSS,               ///< Energy-dependent hard-sphere scattering
    HSS,                ///< Hard-Sphere Stochastic (spherical scattering)
    InteractionPotentialModel, ///< Offline-sampled ion-neutral interaction potential model
    UnknownCollisionModel  ///< Invalid/unknown model
};

/**
 * @brief Orientation sampling mode for the InteractionPotentialModel/IPM collision model
 */
enum class IPMOrientationMode {
    Random, ///< Draw random orientation each collision event
    Fixed   ///< Use a fixed orientation index (ipm_fixed_orientation_index)
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_PHYSICS_ENUMS_H

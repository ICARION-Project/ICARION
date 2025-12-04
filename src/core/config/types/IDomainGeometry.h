// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "core/types/Vec3.h"

namespace ICARION::config {

/**
 * @brief Interface for domain geometry/boundary handling
 *
 * Provides geometry-specific queries and coordinate transforms.
 */
class IDomainGeometry {
public:
    virtual ~IDomainGeometry() = default;

    /// Check whether a global position lies inside the domain.
    virtual bool contains(const Vec3& global_pos) const = 0;

    /// Transform position from global to local coordinates.
    virtual Vec3 global_to_local_pos(const Vec3& global_pos) const = 0;

    /// Transform velocity from global to local coordinates.
    virtual Vec3 global_to_local_vel(const Vec3& global_vel) const = 0;

    /// Transform position from local to global coordinates.
    virtual Vec3 local_to_global_pos(const Vec3& local_pos) const = 0;

    /// Transform velocity from local to global coordinates.
    virtual Vec3 local_to_global_vel(const Vec3& local_vel) const = 0;
};

} // namespace ICARION::config

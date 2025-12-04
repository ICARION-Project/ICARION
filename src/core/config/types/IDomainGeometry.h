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

    /// Domain length [m] along local z.
    virtual double length() const = 0;

    /// Characteristic radius [m] (cylindrical or outer radius).
    virtual double radius() const = 0;

    /// Exit aperture radius [m] (0 if not applicable).
    virtual double end_aperture() const = 0;

    /// Compute an inward-pointing surface normal at a local point on the boundary.
    virtual Vec3 surface_normal(const Vec3& local_pos) const = 0;

    /// Compute first boundary intersection along segment; returns true if hit.
    virtual bool first_boundary_intersection(const Vec3& start_local,
                                             const Vec3& end_local,
                                             Vec3& intersection_local) const = 0;

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

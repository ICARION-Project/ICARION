// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "core/types/Vec3.h"

namespace ICARION::config {

/**
 * @brief Interface for electric/magnetic field models
 *
 * Provides field sampling at a given position/time. Implementations may
 * represent analytical fields or imported maps.
 */
class IFieldModel {
public:
    virtual ~IFieldModel() = default;

    /// Evaluate electric field at global position/time.
    virtual Vec3 E(const Vec3& global_pos, double t) const = 0;

    /// Evaluate magnetic field at global position/time (default zero).
    virtual Vec3 B(const Vec3& global_pos, double t) const { (void)global_pos; (void)t; return Vec3{0.0, 0.0, 0.0}; }

    /// Whether a magnetic field is provided.
    virtual bool has_B() const { return false; }
};

} // namespace ICARION::config

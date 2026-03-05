// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include <vector>
#include <chrono>

namespace ICARION {
namespace fields {

/**
 * @brief Simple 3D vector for field representation
 */
struct Vec3d {
    double x{0.0};
    double y{0.0};
    double z{0.0};
    
    Vec3d() = default;
    Vec3d(double x_, double y_, double z_) : x(x_), y(y_), z(z_) {}
};

/**
 * @brief Structured 3D grid for field arrays
 */
struct Grid3DSnapshot {
    int nx{0}, ny{0}, nz{0};
    double dx{1.0}, dy{1.0}, dz{1.0};
    Vec3d origin{0.0, 0.0, 0.0};
    
    // Flattened field arrays of size nx*ny*nz
    std::vector<double> Ex, Ey, Ez;
};

/**
 * @brief Field snapshot with metadata
 */
struct FieldSnapshot {
    Grid3DSnapshot grid;
    std::chrono::steady_clock::time_point timestamp;
    int version{0};
};

/**
 * @brief Sample electric field using trilinear interpolation (host-only helper)
 * 
 * @param snapshot Field snapshot containing grid data
 * @param pos Position to sample at (in meters)
 * @return Interpolated field E(pos) in V/m
 * 
 * Uses trilinear interpolation between the 8 nearest grid points. Clamps to grid
 * bounds if position is outside. Intended for simple host-side sampling of stored
 * snapshots; not wired into the main simulation loop.
 */
Vec3d sample_field(const FieldSnapshot& snapshot, const Vec3d& pos);

}  // namespace fields
}  // namespace ICARION

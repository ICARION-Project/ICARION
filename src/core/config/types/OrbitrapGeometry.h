// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once

#include "IDomainGeometry.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/Grid3D.h"
#include <cmath>
#include <array>
#include <limits>
#include <algorithm>

namespace ICARION::config {

/**
 * @brief Orbitrap-like radial potential geometry
 *
 * Uses Orbitrap parameters from DomainConfig. Provides global/local transforms
 * and radial/axial containment checks (outer/inner radius). The detailed
 * potential shape is handled in field models.
 */
class OrbitrapGeometry : public IDomainGeometry {
public:
    explicit OrbitrapGeometry(const DomainConfig& cfg)
        : origin_(cfg.geometry.origin_m),
          length_(cfg.geometry.length_m),
          radius_out_(cfg.geometry.radius_out_m),
          radius_in_(cfg.geometry.radius_in_m),
          radius_char_(cfg.geometry.radius_char_m),
          R_g2l_(cfg.rotation_global_to_local),
          R_l2g_(cfg.rotation_local_to_global),
          inner_potential_(cfg.fields.dc.axial_V.constant_value.value_or(0.0)),
          outer_potential_(0.0) {}

    bool contains(const Vec3& global_pos) const override {
        Vec3 local = global_to_local_pos(global_pos);
        const double z_abs = std::fabs(local.z);
        const double z_max = 0.5 * length_;
        if (z_abs > z_max + EPSILON) return false;
        double r2 = local.x * local.x + local.y * local.y;
        const double r = std::sqrt(r2);
        // Radial corridor from hyperlogarithmic electrodes
        const double r_in_allowed = orbitrap_r_for_z(z_abs, radius_in_, radius_char_);
        const double r_out_allowed = orbitrap_r_for_z(z_abs, radius_out_, radius_char_);

        if (!(r_in_allowed > 0.0 && r_out_allowed > r_in_allowed)) {
            return false;
        }
        return (r >= r_in_allowed - EPSILON) && (r <= r_out_allowed + EPSILON);
    }

    Vec3 global_to_local_pos(const Vec3& global_pos) const override {
        Vec3 shifted = global_pos - origin_;
        return {
            R_g2l_.m[0][0] * shifted.x + R_g2l_.m[0][1] * shifted.y + R_g2l_.m[0][2] * shifted.z,
            R_g2l_.m[1][0] * shifted.x + R_g2l_.m[1][1] * shifted.y + R_g2l_.m[1][2] * shifted.z,
            R_g2l_.m[2][0] * shifted.x + R_g2l_.m[2][1] * shifted.y + R_g2l_.m[2][2] * shifted.z
        };
    }

    Vec3 global_to_local_vel(const Vec3& global_vel) const override {
        return {
            R_g2l_.m[0][0] * global_vel.x + R_g2l_.m[0][1] * global_vel.y + R_g2l_.m[0][2] * global_vel.z,
            R_g2l_.m[1][0] * global_vel.x + R_g2l_.m[1][1] * global_vel.y + R_g2l_.m[1][2] * global_vel.z,
            R_g2l_.m[2][0] * global_vel.x + R_g2l_.m[2][1] * global_vel.y + R_g2l_.m[2][2] * global_vel.z
        };
    }

    Vec3 local_to_global_pos(const Vec3& local_pos) const override {
        Vec3 g = {
            R_l2g_.m[0][0] * local_pos.x + R_l2g_.m[0][1] * local_pos.y + R_l2g_.m[0][2] * local_pos.z,
            R_l2g_.m[1][0] * local_pos.x + R_l2g_.m[1][1] * local_pos.y + R_l2g_.m[1][2] * local_pos.z,
            R_l2g_.m[2][0] * local_pos.x + R_l2g_.m[2][1] * local_pos.y + R_l2g_.m[2][2] * local_pos.z
        };
        return g + origin_;
    }

    Vec3 local_to_global_vel(const Vec3& local_vel) const override {
        return {
            R_l2g_.m[0][0] * local_vel.x + R_l2g_.m[0][1] * local_vel.y + R_l2g_.m[0][2] * local_vel.z,
            R_l2g_.m[1][0] * local_vel.x + R_l2g_.m[1][1] * local_vel.y + R_l2g_.m[1][2] * local_vel.z,
            R_l2g_.m[2][0] * local_vel.x + R_l2g_.m[2][1] * local_vel.y + R_l2g_.m[2][2] * local_vel.z
        };
    }

    double length() const override { return length_; }
    double radius() const override { return radius_out_; }
    double end_aperture() const override { return 0.0; } // Not used for Orbitrap
    
    Vec3 surface_normal(const Vec3& /*local_pos*/) const override {
        // Orbitrap uses midpoint termination in DomainManager; provide axial fallback.
        return Vec3{0.0, 0.0, 1.0};
    }

    bool first_boundary_intersection(const Vec3& start_local,
                                     const Vec3& end_local,
                                     Vec3& intersection_local) const override {
        intersection_local = (start_local + end_local) * 0.5;
        return true;
    }

    BoundingBox bounding_box(double margin) const override {
        const double r = radius_out_ + std::max(0.0, margin);
        const double half_len = 0.5 * length_;
        const double z_min = -half_len - std::max(0.0, margin);
        const double z_max = half_len + std::max(0.0, margin);

        BoundingBox box;
        box.min = Vec3{ std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max(),
                        std::numeric_limits<double>::max() };
        box.max = Vec3{ std::numeric_limits<double>::lowest(),
                        std::numeric_limits<double>::lowest(),
                        std::numeric_limits<double>::lowest() };

        for (double x : { -r, r }) {
            for (double y : { -r, r }) {
                for (double z : { z_min, z_max }) {
                    Vec3 global = local_to_global_pos(Vec3{x, y, z});
                    box.min.x = std::min(box.min.x, global.x);
                    box.min.y = std::min(box.min.y, global.y);
                    box.min.z = std::min(box.min.z, global.z);
                    box.max.x = std::max(box.max.x, global.x);
                    box.max.y = std::max(box.max.y, global.y);
                    box.max.z = std::max(box.max.z, global.z);
                }
            }
        }
        return box;
    }

    void apply_spacecharge_dirichlet(const ::Grid3D& grid,
                                     std::vector<char>& mask,
                                     std::vector<double>& values) const override {
        mask.assign(grid.size(), 0);
        values.assign(grid.size(), 0.0);

        const double tol_r = 0.5 * std::max(grid.dx, grid.dy) + 1e-9;
        const double tol_z = 0.5 * grid.dz + 1e-9;
        const double half_len = 0.5 * length_;

        for (int k = 0; k < grid.Nz; ++k) {
            const double z = grid.origin_m.z + static_cast<double>(k) * grid.dz;
            for (int j = 0; j < grid.Ny; ++j) {
                const double y = grid.origin_m.y + static_cast<double>(j) * grid.dy;
                for (int i = 0; i < grid.Nx; ++i) {
                    const double x = grid.origin_m.x + static_cast<double>(i) * grid.dx;
                    const int idx = grid.index(i, j, k);
                    Vec3 local = global_to_local_pos(Vec3{x, y, z});
                    const double r = std::sqrt(local.x * local.x + local.y * local.y);
                    const bool near_inner = r <= radius_in_ + tol_r;
                    const bool near_outer = r >= radius_out_ - tol_r;
                    const bool beyond_inner = r < radius_in_ - tol_r;
                    const bool beyond_outer = r > radius_out_ + tol_r;
                    const bool axial_cap = std::fabs(local.z) >= (half_len - tol_z);

                    if (near_inner || beyond_inner || axial_cap) {
                        mask[idx] = 1;
                        values[idx] = inner_potential_;
                        continue;
                    }
                    if (near_outer || beyond_outer) {
                        mask[idx] = 1;
                        values[idx] = outer_potential_;
                    }
                }
            }
        }
    }

private:
    // Hyperlogarithmic surface solver (matches DomainManager logic)
    static double orbitrap_surface_residual(double r, double z, double R, double R_m) {
        const double term1 = z * z;
        const double term2 = 0.5 * (r * r - R * R);
        const double term3 = R_m * R_m * std::log(R / r);
        return term1 - term2 - term3;
    }

    static double orbitrap_r_for_z(double z, double R, double R_m) {
        const double z_abs = std::fabs(z);
        const double eps = 1e-10;
        const int max_iter = 80;

        if (z_abs < eps) {
            return R;
        }

        double r_lo = 0.1 * R;
        double r_hi = R;

        double f_lo = orbitrap_surface_residual(r_lo, z_abs, R, R_m);
        double f_hi = orbitrap_surface_residual(r_hi, z_abs, R, R_m);

        int expand_iter = 0;
        while (f_lo * f_hi > 0.0 && expand_iter < 10) {
            r_lo *= 0.5;
            r_hi *= 1.5;
            f_lo = orbitrap_surface_residual(r_lo, z_abs, R, R_m);
            f_hi = orbitrap_surface_residual(r_hi, z_abs, R, R_m);
            expand_iter++;
        }

        double r_mid = R;
        for (int i = 0; i < max_iter; ++i) {
            r_mid = 0.5 * (r_lo + r_hi);
            double f_mid = orbitrap_surface_residual(r_mid, z_abs, R, R_m);
            if (f_lo * f_mid <= 0.0) {
                r_hi = r_mid;
                f_hi = f_mid;
            } else {
                r_lo = r_mid;
                f_lo = f_mid;
            }
            if (std::fabs(f_mid) < eps) break;
        }
        return r_mid;
    }

    Vec3 origin_;
    double length_;
    double radius_out_;
    double radius_in_;
    double radius_char_;
    Mat3 R_g2l_;
    Mat3 R_l2g_;
    double inner_potential_;
    double outer_potential_;

    static constexpr double EPSILON = 1e-12;
};

} // namespace ICARION::config

// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IDomainGeometry.h"
#include "core/config/types/DomainConfig.h"
#include <cmath>
#include <array>
#include <limits>
#include <algorithm>

namespace ICARION::config {

/**
 * @brief Cylindrical domain geometry (radius/length, with rotation/origin)
 *
 * Uses parameters from DomainConfig and handles global/local transforms and
 * basic radial/axial containment checks with a small epsilon tolerance.
 */
class CylindricalGeometry : public IDomainGeometry {
public:
    explicit CylindricalGeometry(const DomainConfig& cfg)
        : origin_(cfg.geometry.origin_m),
          radius_(cfg.geometry.radius_m),
          length_(cfg.geometry.length_m),
          end_aperture_(cfg.geometry.end_aperture_m),
          R_g2l_(cfg.rotation_global_to_local),
          R_l2g_(cfg.rotation_local_to_global) {}

    bool contains(const Vec3& global_pos) const override {
        Vec3 local = global_to_local_pos(global_pos);
        if (local.z < -EPSILON || local.z >= length_ + EPSILON) return false;
        double r2 = local.x * local.x + local.y * local.y;
        return r2 <= (radius_ + EPSILON) * (radius_ + EPSILON);
    }

    double length() const override { return length_; }
    double radius() const override { return radius_; }
    double end_aperture() const override { return end_aperture_; }
    
    Vec3 surface_normal(const Vec3& local_pos) const override {
        const double r = std::sqrt(local_pos.x * local_pos.x + local_pos.y * local_pos.y);
        const double EPS = 1e-9;
        if (std::abs(r - radius_) < EPS && r > NUMERICAL_ZERO) {
            return Vec3{-local_pos.x / r, -local_pos.y / r, 0.0};
        }
        if (std::abs(local_pos.z) < EPS) {
            return Vec3{0.0, 0.0, 1.0};
        }
        if (std::abs(local_pos.z - length_) < EPS) {
            return Vec3{0.0, 0.0, -1.0};
        }
        // Fallback
        return Vec3{0.0, 0.0, 1.0};
    }

    bool first_boundary_intersection(const Vec3& start_local,
                                     const Vec3& end_local,
                                     Vec3& intersection_local) const override {
        Vec3 d = end_local - start_local;
        double t_min = 1.0;
        bool hit = false;

        // Radial cylinder
        const double a = d.x*d.x + d.y*d.y;
        const double b = 2.0*(start_local.x*d.x + start_local.y*d.y);
        const double c = start_local.x*start_local.x + start_local.y*start_local.y - radius_*radius_;

        if (a > NUMERICAL_ZERO) {
            double disc = b*b - 4.0*a*c;
            if (disc >= 0.0) {
                double sqrt_disc = std::sqrt(disc);
                double t1 = (-b + sqrt_disc) / (2.0*a);
                double t2 = (-b - sqrt_disc) / (2.0*a);
                for (double t : {t1, t2}) {
                    if (t > 0.0 && t <= 1.0) {
                        Vec3 cand = start_local + d * t;
                        if (cand.z >= -EPSILON && cand.z <= length_ + EPSILON && t < t_min) {
                            t_min = t;
                            intersection_local = cand;
                            hit = true;
                        }
                    }
                }
            }
        }

        // Entrance plane z=0
        if (std::abs(d.z) > NUMERICAL_ZERO) {
            double t = (-start_local.z) / d.z;
            if (t > 0.0 && t <= 1.0 && t < t_min) {
                Vec3 cand = start_local + d * t;
                double r2 = cand.x*cand.x + cand.y*cand.y;
                if (r2 <= radius_*radius_) {
                    t_min = t;
                    intersection_local = cand;
                    hit = true;
                }
            }
        }

        // Exit plane z=length_
        if (std::abs(d.z) > NUMERICAL_ZERO) {
            double t = (length_ - start_local.z) / d.z;
            if (t > 0.0 && t <= 1.0 && t < t_min) {
                Vec3 cand = start_local + d * t;
                double r2 = cand.x*cand.x + cand.y*cand.y;
                double limit = (end_aperture_ > 0.0) ? end_aperture_ : radius_;
                if (r2 <= limit*limit) {
                    t_min = t;
                    intersection_local = cand;
                    hit = true;
                }
            }
        }

        if (!hit) {
            intersection_local = end_local;
        }
        return hit;
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

    BoundingBox bounding_box(double margin) const override {
        const double r = radius_ + std::max(0.0, margin);
        const double z_min = -std::max(0.0, margin);
        const double z_max = length_ + std::max(0.0, margin);

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

private:
    Vec3 origin_;
    double radius_;
    double length_;
    double end_aperture_;
    Mat3 R_g2l_;
    Mat3 R_l2g_;

    static constexpr double EPSILON = 1e-12;
    static constexpr double NUMERICAL_ZERO = 1e-15;
};

} // namespace ICARION::config

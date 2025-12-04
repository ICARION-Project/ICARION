// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IDomainGeometry.h"
#include "core/config/types/DomainConfig.h"
#include <cmath>

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
          epsilon_(cfg.geometry.boundary_epsilon),
          R_g2l_(cfg.geometry.rotation_global_to_local),
          R_l2g_(cfg.geometry.rotation_local_to_global) {}

    bool contains(const Vec3& global_pos) const override {
        Vec3 local = global_to_local_pos(global_pos);
        if (local.z < -epsilon_ || local.z > length_ + epsilon_) return false;
        double r2 = local.x * local.x + local.y * local.y;
        // Basic radial check: inside outer radius and outside inner radius if defined
        if (radius_in_ > 0.0 && r2 < (radius_in_ - epsilon_) * (radius_in_ - epsilon_)) {
            return false;
        }
        return r2 <= (radius_out_ + epsilon_) * (radius_out_ + epsilon_);
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

private:
    Vec3 origin_;
    double length_;
    double radius_out_;
    double radius_in_;
    double radius_char_;
    double epsilon_;
    RotationMatrix R_g2l_;
    RotationMatrix R_l2g_;
};

} // namespace ICARION::config

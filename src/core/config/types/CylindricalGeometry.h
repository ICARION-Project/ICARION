// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "IDomainGeometry.h"
#include "core/config/types/DomainConfig.h"

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
          epsilon_(cfg.geometry.boundary_epsilon),
          R_g2l_(cfg.geometry.rotation_global_to_local),
          R_l2g_(cfg.geometry.rotation_local_to_global) {}

    bool contains(const Vec3& global_pos) const override {
        Vec3 local = global_to_local_pos(global_pos);
        if (local.z < -epsilon_ || local.z > length_ + epsilon_) return false;
        double r2 = local.x * local.x + local.y * local.y;
        return r2 <= (radius_ + epsilon_) * (radius_ + epsilon_);
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
    double radius_;
    double length_;
    double epsilon_;
    RotationMatrix R_g2l_;
    RotationMatrix R_l2g_;
};

} // namespace ICARION::config

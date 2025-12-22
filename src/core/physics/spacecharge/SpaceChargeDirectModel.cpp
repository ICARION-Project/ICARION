// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "SpaceChargeDirectModel.h"
#include "utils/constants.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace ICARION::physics {

SpaceChargeDirectModel::SpaceChargeDirectModel(double softening_m)
    : softening_m_(softening_m) {
    if (softening_m_ < 0.0) {
        throw std::invalid_argument("SpaceChargeDirectModel: softening must be >= 0");
    }
}

void SpaceChargeDirectModel::update_fields(const core::IonEnsemble& ions, double /*time_s*/) {
    const size_t n = ions.size();
    electric_field_.assign(n, core::Vec3{0.0, 0.0, 0.0});
    if (n <= 1) {
        return;
    }

    const auto* pos_x = ions.pos_x_data();
    const auto* pos_y = ions.pos_y_data();
    const auto* pos_z = ions.pos_z_data();
    const auto* charge = ions.charge_data();
    const auto* active = ions.active_data();

    const double epsilon_sq = softening_m_ * softening_m_;
    const double min_distance_sq = std::max(softening_m_ * softening_m_ * 1e-4, 1e-30);

    for (size_t i = 0; i < n; ++i) {
        if (active[i] == 0) continue;

        double ex = 0.0, ey = 0.0, ez = 0.0;
        const double xi = pos_x[i];
        const double yi = pos_y[i];
        const double zi = pos_z[i];

        for (size_t j = 0; j < n; ++j) {
            if (j == i || active[j] == 0) continue;

            const double dx = xi - pos_x[j];
            const double dy = yi - pos_y[j];
            const double dz = zi - pos_z[j];
            const double r_sq = dx * dx + dy * dy + dz * dz;
            if (r_sq < min_distance_sq) continue;

            const double r = std::sqrt(r_sq);
            const double r_eff_sq = r_sq + epsilon_sq;
            const double inv_r = 1.0 / r;

            // Electric field from charge j on i: E = k * q_j / r_eff^2 * r_hat
            const double scale = COULOMB_CONST * charge[j] / r_eff_sq * inv_r;
            ex += dx * scale;
            ey += dy * scale;
            ez += dz * scale;
        }

        electric_field_[i] = core::Vec3{ex, ey, ez};
    }
}

core::Vec3 SpaceChargeDirectModel::sample_electric_field(std::size_t ion_idx) const {
    if (ion_idx >= electric_field_.size()) {
        return core::Vec3{0.0, 0.0, 0.0};
    }
    return electric_field_[ion_idx];
}

} // namespace ICARION::physics

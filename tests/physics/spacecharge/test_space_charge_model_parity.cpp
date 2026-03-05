// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <vector>
#include <cmath>
#include <limits>
#include <memory>

#include "core/physics/spacecharge/SpaceChargeDirectModel.h"
#include "core/physics/spacecharge/SpaceChargeGridModel.h"
#include "core/config/types/CylindricalGeometry.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonEnsemble.h"
#include "core/utils/mathUtils.h"
#include "utils/constants.h"

using ICARION::core::IonEnsemble;
using ICARION::core::IonState;
using ICARION::core::Vec3;
using ICARION::config::CylindricalGeometry;
using ICARION::config::DomainConfig;
using ICARION::physics::SpaceChargeDirectModel;
using ICARION::physics::SpaceChargeGridModel;

namespace {

std::vector<IonState> make_test_ions() {
    std::vector<IonState> ions;
    ions.reserve(64);

    const double base = -0.0015;
    const double step = 0.0010;
    double jitter = 0.0;

    for (int ix = 0; ix < 4; ++ix) {
        for (int iy = 0; iy < 4; ++iy) {
            for (int iz = 0; iz < 4; ++iz) {
                IonState ion;
                ion.pos = Vec3{
                    base + ix * step + jitter * 1e-4,
                    base + iy * step - jitter * 5e-5,
                    base + iz * step + jitter * 2e-4};
                ion.vel = Vec3{0.0, 0.0, 0.0};
                ion.ion_charge_C = ELEM_CHARGE_C;
                ion.mass_kg = AMU_TO_KG;
                ion.active = true;
                ion.born = true;
                ions.push_back(ion);
                jitter += 1.0;
            }
        }
    }
    return ions;
}

DomainConfig make_cylindrical_domain() {
    DomainConfig dom;
    dom.name = "parity_domain";
    dom.geometry.length_m = 0.05;
    dom.geometry.radius_m = 0.01;
    dom.geometry.origin_m = Vec3{0.0, 0.0, -0.025};
    dom.geometry.end_aperture_m = 0.01;
    dom.fields.dc.axial_V.constant_value = 0.0;
    dom.fields.dc.radial_V.constant_value = 0.0;
    dom.rotation_global_to_local = Mat3::identity();
    dom.rotation_local_to_global = Mat3::identity();
    return dom;
}

double vec_norm(const Vec3& v) {
    return std::sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

} // namespace

TEST_CASE("SpaceChargeGridModel approximates SpaceChargeDirectModel", "[spacecharge][grid][parity]") {
    auto ions_legacy = make_test_ions();
    IonEnsemble ensemble = IonEnsemble::from_legacy(ions_legacy);

    SpaceChargeDirectModel direct_model(1e-10);
    direct_model.update_fields(ensemble, 0.0);

    DomainConfig domain = make_cylindrical_domain();
    auto geometry = std::make_unique<CylindricalGeometry>(domain);
    SpaceChargeGridModel grid_model(domain, std::move(geometry), 0.001, 96);
    grid_model.update_fields(ensemble, 0.0);

    double max_relative_error = 0.0;
    double max_relative_significant = 0.0;
    double max_absolute_small = 0.0;
    double min_direct_mag = std::numeric_limits<double>::max();
    double max_direct_mag = 0.0;
    double max_grid_mag = 0.0;
    double sum_direct_mag = 0.0;
    double sum_diff_sq = 0.0;
    constexpr double SIGNIFICANCE = 0.002;
    size_t worst_idx = 0;
    Vec3 worst_direct{0.0, 0.0, 0.0};
    Vec3 worst_grid{0.0, 0.0, 0.0};
    for (size_t i = 0; i < ensemble.size(); ++i) {
        Vec3 E_direct = direct_model.sample_electric_field(i);
        Vec3 E_grid = grid_model.sample_electric_field(i);
        Vec3 diff = E_grid - E_direct;

        const double direct_mag = vec_norm(E_direct);
        const double diff_mag = vec_norm(diff);
        const double grid_mag = vec_norm(E_grid);
        const double denom = std::max(1e-12, direct_mag);
        const double rel_error = diff_mag / denom;
        max_relative_error = std::max(max_relative_error, rel_error);
        if (direct_mag > SIGNIFICANCE) {
            max_relative_significant = std::max(max_relative_significant, rel_error);
            if (rel_error >= max_relative_significant) {
                worst_idx = i;
                worst_direct = E_direct;
                worst_grid = E_grid;
            }
        } else {
            max_absolute_small = std::max(max_absolute_small, diff_mag);
        }
        sum_direct_mag += direct_mag;
        sum_diff_sq += diff_mag * diff_mag;

        min_direct_mag = std::min(min_direct_mag, direct_mag);
        max_direct_mag = std::max(max_direct_mag, direct_mag);
        max_grid_mag = std::max(max_grid_mag, grid_mag);
    }
    INFO("Direct field magnitude range: [" << min_direct_mag << ", " << max_direct_mag
         << "], max grid=" << max_grid_mag);

    const double avg_direct_mag = sum_direct_mag / static_cast<double>(ensemble.size());
    const double rms_diff = std::sqrt(sum_diff_sq / static_cast<double>(ensemble.size()));
    const double rms_relative = rms_diff / std::max(1e-12, avg_direct_mag);

    INFO("Maximum relative error grid vs direct: " << max_relative_error);
    INFO("RMS relative error grid vs direct: " << rms_relative);
    INFO("Max relative error for |E|>2e-3: " << max_relative_significant);
    INFO("Max absolute error for |E|<=2e-3: " << max_absolute_small);
    INFO("Worst ion index: " << worst_idx << " direct=(" << worst_direct.x << "," << worst_direct.y << "," << worst_direct.z
         << ") grid=(" << worst_grid.x << "," << worst_grid.y << "," << worst_grid.z << ")");

    REQUIRE(max_relative_significant < 0.70);
    REQUIRE(max_absolute_small < 1.5e-3);
    REQUIRE(rms_relative < 0.30);
}

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "SpaceChargeGridModel.h"
#include "core/config/types/DomainConfig.h"
#include <algorithm>
#include <cmath>

namespace ICARION::physics {

namespace {

struct GridParams {
    int nx;
    int ny;
    int nz;
    double dx;
    double dy;
    double dz;
    core::Vec3 origin;
};

GridParams make_grid_params(const config::BoundingBox& bbox,
                            int target_resolution,
                            double padding) {
    core::Vec3 min = bbox.min;
    core::Vec3 max = bbox.max;
    min.x -= padding;
    min.y -= padding;
    min.z -= padding;
    max.x += padding;
    max.y += padding;
    max.z += padding;
    core::Vec3 extent = {std::max(max.x - min.x, 1e-6),
                         std::max(max.y - min.y, 1e-6),
                         std::max(max.z - min.z, 1e-6)};

    const double largest = std::max({extent.x, extent.y, extent.z});
    const double base_cell = largest / static_cast<double>(std::max(target_resolution, 8));

    auto compute_cells = [&](double size) -> int {
        int cells = static_cast<int>(std::ceil(size / base_cell));
        return std::max(cells, 8);
    };

    GridParams params;
    params.nx = compute_cells(extent.x);
    params.ny = compute_cells(extent.y);
    params.nz = compute_cells(extent.z);
    params.dx = extent.x / static_cast<double>(params.nx);
    params.dy = extent.y / static_cast<double>(params.ny);
    params.dz = extent.z / static_cast<double>(params.nz);
    params.origin = min;
    return params;
}

} // namespace

SpaceChargeGridModel::SpaceChargeGridModel(const config::DomainConfig& domain,
                                           std::unique_ptr<config::IDomainGeometry> geometry,
                                           double padding,
                                           int target_resolution)
    : domain_(domain),
      geometry_(std::move(geometry)),
      padding_m_(padding),
      resolution_hint_(target_resolution) {
    auto bbox = geometry_->bounding_box(padding_m_);
    GridParams params = make_grid_params(bbox, resolution_hint_, padding_m_);
    ICARION::log::Logger::main()->info(
        "Space charge grid ({}) → {}x{}x{} cells, cell ≈ {:.3e} m",
        domain.name, params.nx, params.ny, params.nz, params.dx);
    solver_ = std::make_unique<SpaceChargeSolver>(
        params.nx, params.ny, params.nz,
        params.dx, params.dy, params.dz,
        params.origin);
    solver_->setGeometryMask(geometry_.get());

    std::vector<char> mask(solver_->grid().size(), 0);
    std::vector<double> values(mask.size(), 0.0);
    geometry_->apply_spacecharge_dirichlet(solver_->grid(), mask, values);
    if (!mask.empty() && mask.size() == values.size()) {
        solver_->setDirichletMask(mask);
        solver_->setDirichletValues(values);
    }
}

void SpaceChargeGridModel::update_fields(const core::IonEnsemble& ions, double /*time_s*/) {
    solver_->update(ions);
    const size_t n = ions.size();
    cached_field_.assign(n, core::Vec3{0.0, 0.0, 0.0});
    const auto* pos_x = ions.pos_x_data();
    const auto* pos_y = ions.pos_y_data();
    const auto* pos_z = ions.pos_z_data();

    for (size_t i = 0; i < n; ++i) {
        core::Vec3 pos(pos_x[i], pos_y[i], pos_z[i]);
        cached_field_[i] = solver_->fieldAt(pos);
    }
}

core::Vec3 SpaceChargeGridModel::sample_electric_field(std::size_t ion_idx) const {
    if (ion_idx >= cached_field_.size()) {
        return core::Vec3{0.0, 0.0, 0.0};
    }
    return cached_field_[ion_idx];
}

} // namespace ICARION::physics

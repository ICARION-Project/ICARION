// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "SpaceChargeModelFactory.h"
#include "SpaceChargeDirectModel.h"
#include "SpaceChargeGridModel.h"
#include "SpaceChargeGPUModel.h"
#include "core/config/types/CylindricalGeometry.h"
#include "core/config/types/OrbitrapGeometry.h"
#include "core/log/Logger.h"
#include <algorithm>

#ifdef ICARION_USE_GPU
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/spacecharge/GPUSpaceChargeP3M.h"
#endif

namespace ICARION::physics {

namespace {

std::unique_ptr<config::IDomainGeometry> make_geometry(const config::DomainConfig& domain) {
    if (domain.instrument == config::Instrument::Orbitrap) {
        return std::make_unique<config::OrbitrapGeometry>(domain);
    }
    return std::make_unique<config::CylindricalGeometry>(domain);
}

#ifdef ICARION_USE_GPU
SpaceChargeModelPtr try_create_gpu_model(const config::FullConfig& config,
                                         const config::DomainConfig& domain,
                                         std::size_t ion_count) {
    constexpr std::size_t GPU_THRESHOLD = 2000;
    constexpr double TARGET_CELL_SIZE = 30e-6;
    constexpr int MIN_RES = 32;
    constexpr int MAX_RES = 256;

    if (!config.physics.enable_space_charge_gpu || ion_count < GPU_THRESHOLD) {
        return nullptr;
    }

    auto gpu_ctx_unique = icarion::gpu::GPUContext::create(0);
    if (!gpu_ctx_unique || !gpu_ctx_unique->is_valid()) {
        ICARION::log::Logger::main()->info(
            "Space charge: GPU context unavailable, falling back to CPU models");
        return nullptr;
    }
    auto gpu_ctx = std::shared_ptr<icarion::gpu::GPUContext>(std::move(gpu_ctx_unique));

    auto geometry = make_geometry(domain);
    auto bbox = geometry->bounding_box(SpaceChargeModelFactory::GRID_PADDING_M);
    Vec3 extent = bbox.max - bbox.min;
    auto ensure_extent = [](double value) {
        constexpr double MIN_EXTENT = 0.01;
        if (value <= 0.0) {
            return MIN_EXTENT;
        }
        return std::max(value, MIN_EXTENT);
    };
    extent.x = ensure_extent(extent.x);
    extent.y = ensure_extent(extent.y);
    extent.z = ensure_extent(extent.z);

    auto cells_from_extent = [&](double axis_length) -> int {
        const double desired = axis_length / TARGET_CELL_SIZE;
        return std::clamp(static_cast<int>(std::ceil(desired)), MIN_RES, MAX_RES);
    };

    icarion::gpu::GPUSpaceChargeP3M::Config cfg;
    cfg.grid_nx = cells_from_extent(extent.x);
    cfg.grid_ny = cells_from_extent(extent.y);
    cfg.grid_nz = cells_from_extent(extent.z);
    cfg.domain_min = bbox.min;
    cfg.domain_max = bbox.max;

    auto solver_unique = icarion::gpu::GPUSpaceChargeP3M::create(*gpu_ctx, cfg);
    if (!solver_unique) {
        ICARION::log::Logger::main()->warn(
            "Space charge: GPU solver creation failed for {}", domain.name);
        return nullptr;
    }

    auto solver = std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M>(std::move(solver_unique));
    auto model = std::make_shared<SpaceChargeGPUModel>(gpu_ctx, solver, domain.name);
    if (!model->is_available()) {
        ICARION::log::Logger::main()->warn(
            "Space charge: GPU solver unavailable for {}", domain.name);
        return nullptr;
    }

    ICARION::log::Logger::main()->info(
        "Space charge: GPU P³M model enabled for {} ({}×{}×{} grid)",
        domain.name, cfg.grid_nx, cfg.grid_ny, cfg.grid_nz);
    return model;
}
#endif

} // namespace

SpaceChargeModelPtr SpaceChargeModelFactory::create(const config::FullConfig& config,
                                                    const config::DomainConfig& domain,
                                                    std::size_t ion_count) {
    if (!config.physics.enable_space_charge || ion_count == 0) {
        return nullptr;
    }

    const auto model_type = config.physics.space_charge_model_type;

    if (model_type == config::SpaceChargeModel::Direct) {
        return std::make_shared<SpaceChargeDirectModel>(DIRECT_SOFTENING_M);
    }

    if (model_type == config::SpaceChargeModel::Grid) {
        auto geometry = make_geometry(domain);
        return std::make_shared<SpaceChargeGridModel>(
            domain, std::move(geometry), GRID_PADDING_M, GRID_TARGET_RESOLUTION);
    }

#ifdef ICARION_USE_GPU
    if (model_type == config::SpaceChargeModel::GPU || config.physics.enable_space_charge_gpu) {
        if (auto gpu_model = try_create_gpu_model(config, domain, ion_count)) {
            return gpu_model;
        }
        if (model_type == config::SpaceChargeModel::GPU) {
            return nullptr;
        }
    }
#else
    if (model_type == config::SpaceChargeModel::GPU) {
        return nullptr;
    }
#endif

    if (ion_count < DIRECT_THRESHOLD) {
        return std::make_shared<SpaceChargeDirectModel>(DIRECT_SOFTENING_M);
    }

    auto geometry = make_geometry(domain);
    return std::make_shared<SpaceChargeGridModel>(
        domain, std::move(geometry), GRID_PADDING_M, GRID_TARGET_RESOLUTION);
}

} // namespace ICARION::physics

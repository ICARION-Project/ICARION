// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once

#include <memory>
#include "core/config/types/FullConfig.h"
#include "core/physics/spacecharge/ISpaceChargeModel.h"

namespace ICARION::physics {

/**
 * @brief Factory for creating per-domain space-charge models (Direct/Grid/GPU).
 *
 * Selection order:
 *  1. GPU P³M (`SpaceChargeGPUModel`) when CUDA build + `enable_space_charge_gpu`.
 *  2. Grid-based Poisson solver (`SpaceChargeGridModel`) using IDomainGeometry bounds.
 *  3. Direct Coulomb model (`SpaceChargeDirectModel`) for low ion counts.
 */
class SpaceChargeModelFactory {
public:
    static constexpr std::size_t DIRECT_THRESHOLD = 1000;
    static constexpr std::size_t GRID_THRESHOLD = 1000;
    static constexpr int GRID_TARGET_RESOLUTION = 64;
    static constexpr double GRID_PADDING_M = 0.0;

    /**
     * @param config Full configuration (for global flags).
     * @param domain Domain config owning the model.
     * @param ion_count Number of ions assigned to the domain.
     * @return SpaceChargeModelPtr (nullptr only when space charge is disabled or model creation fails).
     */
    static SpaceChargeModelPtr create(const config::FullConfig& config,
                                      const config::DomainConfig& domain,
                                      std::size_t ion_count);
};

} // namespace ICARION::physics

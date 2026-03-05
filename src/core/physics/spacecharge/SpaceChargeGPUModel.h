// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#include "core/physics/spacecharge/ISpaceChargeModel.h"
#include "core/log/Logger.h"
#include <memory>
#include <string>
#include <vector>

#ifdef ICARION_USE_GPU
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/spacecharge/GPUSpaceChargeP3M.h"
#endif

namespace ICARION::physics {

/**
 * @brief ISpaceChargeModel wrapper around the experimental GPU P³M solver.
 *
 * When ICARION_USE_GPU is not defined, the model acts as a stub that returns zero
 * fields so CPU-only builds remain functional.
 */
class SpaceChargeGPUModel : public ISpaceChargeModel {
public:
#ifdef ICARION_USE_GPU
    SpaceChargeGPUModel(std::shared_ptr<icarion::gpu::GPUContext> context,
                        std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> solver,
                        std::string label);
#else
    explicit SpaceChargeGPUModel(std::string label = "SpaceChargeGPUModel (stub)");
#endif

    void update_fields(const core::IonEnsemble& ions, double time_s) override;
    core::Vec3 sample_electric_field(std::size_t ion_idx) const override;
    std::string name() const override { return "SpaceChargeGPUModel"; }

    bool is_available() const { return available_; }

private:
#ifdef ICARION_USE_GPU
    std::shared_ptr<icarion::gpu::GPUContext> context_;
    std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> solver_;
#endif
    std::vector<core::Vec3> cached_field_;
    bool available_ = false;
    std::string label_;
};

} // namespace ICARION::physics


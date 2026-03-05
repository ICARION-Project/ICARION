// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "SpaceChargeGPUModel.h"

namespace ICARION::physics {

#ifdef ICARION_USE_GPU
SpaceChargeGPUModel::SpaceChargeGPUModel(
    std::shared_ptr<icarion::gpu::GPUContext> context,
    std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> solver,
    std::string label)
    : context_(std::move(context)),
      solver_(std::move(solver)),
      label_(std::move(label)) {
    available_ = static_cast<bool>(context_) && static_cast<bool>(solver_);
}
#else
SpaceChargeGPUModel::SpaceChargeGPUModel(std::string label)
    : label_(std::move(label)) {
    available_ = false;
}
#endif

void SpaceChargeGPUModel::update_fields(const core::IonEnsemble& ions, double /*time_s*/) {
    const size_t n = ions.size();
    cached_field_.assign(n, core::Vec3{0.0, 0.0, 0.0});

#ifdef ICARION_USE_GPU
    if (!available_ || !solver_) {
        return;
    }
    if (!solver_->compute_space_charge_field(ions, cached_field_)) {
        available_ = false;
        ICARION::log::Logger::main()->warn(
            "SpaceChargeGPUModel: compute_space_charge_field failed for {}", label_);
        std::fill(cached_field_.begin(), cached_field_.end(), core::Vec3{0.0, 0.0, 0.0});
    }
#else
    (void)ions;
#endif
}

core::Vec3 SpaceChargeGPUModel::sample_electric_field(std::size_t ion_idx) const {
    if (ion_idx >= cached_field_.size()) {
        return core::Vec3{0.0, 0.0, 0.0};
    }
    return cached_field_[ion_idx];
}

} // namespace ICARION::physics


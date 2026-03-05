// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once

#include "core/physics/spacecharge/ISpaceChargeModel.h"
#include <vector>

namespace ICARION::physics {

/**
 * @brief Direct N-body space-charge model (O(N²)).
 *
 * Computes electric fields by summing Coulomb contributions from every other ion.
 * Intended for small ion counts or validation of grid/GPU implementations.
 */
class SpaceChargeDirectModel : public ISpaceChargeModel {
public:
    explicit SpaceChargeDirectModel(double softening_m = 1e-10);

    void update_fields(const core::IonEnsemble& ions, double time_s) override;
    core::Vec3 sample_electric_field(std::size_t ion_idx) const override;
    std::string name() const override { return "SpaceChargeDirectModel"; }

    double softening_length() const { return softening_m_; }

private:
    double softening_m_;
    std::vector<core::Vec3> electric_field_;  // cached per-ion field
};

} // namespace ICARION::physics

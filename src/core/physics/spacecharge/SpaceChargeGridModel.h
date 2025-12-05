// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include "core/physics/spacecharge/ISpaceChargeModel.h"
#include "core/physics/spacecharge/spaceChargeSolver.h"
#include "core/config/types/IDomainGeometry.h"
#include "core/config/types/DomainConfig.h"
#include "core/log/Logger.h"
#include <memory>
#include <vector>

namespace ICARION::physics {

class SpaceChargeGridModel : public ISpaceChargeModel {
public:
    SpaceChargeGridModel(const config::DomainConfig& domain,
                         std::unique_ptr<config::IDomainGeometry> geometry,
                         double padding = 0.0,
                         int target_resolution = 64);

    void update_fields(const core::IonEnsemble& ions, double time_s) override;
    core::Vec3 sample_electric_field(std::size_t ion_idx) const override;
    std::string name() const override { return "SpaceChargeGridModel"; }

private:
    config::DomainConfig domain_;
    std::unique_ptr<config::IDomainGeometry> geometry_;
    double padding_m_;
    int resolution_hint_;
    std::unique_ptr<SpaceChargeSolver> solver_;
    std::vector<core::Vec3> cached_field_;
};

} // namespace ICARION::physics

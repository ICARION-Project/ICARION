// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "SpaceChargeModelFactory.h"
#include "SpaceChargeDirectModel.h"

namespace ICARION::physics {

SpaceChargeModelPtr SpaceChargeModelFactory::create(const config::FullConfig& config,
                                                    const config::DomainConfig&,
                                                    std::size_t ion_count) {
    if (!config.physics.enable_space_charge || ion_count == 0) {
        return nullptr;
    }

    constexpr std::size_t DIRECT_THRESHOLD = 1000;
    constexpr double DEFAULT_SOFTENING = 1e-10;

    if (ion_count < DIRECT_THRESHOLD) {
        return std::make_shared<SpaceChargeDirectModel>(DEFAULT_SOFTENING);
    }

    return nullptr;
}

} // namespace ICARION::physics

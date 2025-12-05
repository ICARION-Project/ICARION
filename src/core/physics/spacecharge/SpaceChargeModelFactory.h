// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <memory>
#include "core/config/types/FullConfig.h"
#include "core/physics/spacecharge/ISpaceChargeModel.h"

namespace ICARION::physics {

/**
 * @brief Placeholder factory (vNext) for creating space-charge models per domain.
 *
 * Currently returns nullptr; wiring will land in later commits once the new models
 * are implemented. Documented here so other components can start depending on the
 * API without enabling runtime behavior yet.
 */
class SpaceChargeModelFactory {
public:
    /**
     * @param config Full configuration (for global flags).
     * @param domain Domain config owning the model.
     * @param ion_count Number of ions assigned to the domain.
     * @return Unique pointer to ISpaceChargeModel (nullptr for now).
     */
    static SpaceChargeModelPtr create(const config::FullConfig& config,
                                      const config::DomainConfig& domain,
                                      std::size_t ion_count);
};

} // namespace ICARION::physics

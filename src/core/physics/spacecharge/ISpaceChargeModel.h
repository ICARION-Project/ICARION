// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <cstddef>
#include <memory>
#include <string>

#include "core/types/IonEnsemble.h"
#include "core/types/Vec3.h"

namespace ICARION::physics {

/**
 * @brief Strategy interface for space-charge solvers/models.
 *
 * Implementations encapsulate how charge deposition / Poisson solve / GPU kernels are run.
 * SimulationEngine or a dedicated force calls `update_fields()` once per timestep and then
 * samples per-ion electric fields via `sample_electric_field()`.
 */
class ISpaceChargeModel {
public:
    virtual ~ISpaceChargeModel() = default;

    /**
     * @brief Update internal field state based on current ion ensemble.
     *
     * @param ions  SoA ensemble (positions, charges) for the active domain.
     * @param time_s Current simulation time in seconds (for adaptive tolerances etc.).
     */
    virtual void update_fields(const core::IonEnsemble& ions, double time_s) = 0;

    /**
     * @brief Sample electric field for ion @p ion_idx (SoA index).
     *
     * @param ion_idx Index inside the provided IonEnsemble.
     * @return Electric field vector [V/m].
     */
    virtual core::Vec3 sample_electric_field(std::size_t ion_idx) const = 0;

    /**
     * @brief Model identifier (for logging / diagnostics).
     */
    virtual std::string name() const = 0;
};

using SpaceChargeModelPtr = std::shared_ptr<ISpaceChargeModel>;

} // namespace ICARION::physics

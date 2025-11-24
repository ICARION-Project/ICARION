// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "SpaceChargeGrid.h"
#include "ForceContext.h"
#include <stdexcept>
#include <cmath>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ICARION::physics {

SpaceChargeGrid::SpaceChargeGrid(std::shared_ptr<SpaceChargeSolver> solver)
    : solver_(solver),
      last_update_time_(-1.0),
      solver_updated_this_step_(false)
{
    if (!solver_) {
        throw std::invalid_argument(
            "SpaceChargeGrid: Solver cannot be null. "
            "Pass a valid SpaceChargeSolver instance."
        );
    }
}

Vec3 SpaceChargeGrid::compute(
    const IonState& ion,
    double t,
    const ForceContext& ctx
) const {
    // Validate solver
    if (!solver_) {
        return {0, 0, 0};
    }
    
    // Validate context (need all_ions for grid update)
    if (!ctx.all_ions || ctx.all_ions->empty()) {
        return {0, 0, 0};
    }
    
    // Update solver once per timestep (time-based deduplication)
    // Use epsilon comparison for floating-point time
    constexpr double TIME_EPSILON = 1e-12;  // 1 ps tolerance
    
    if (std::abs(t - last_update_time_) > TIME_EPSILON) {
        // New timestep detected - update solver
        #ifdef _OPENMP
        #pragma omp critical(space_charge_solver_update)
        #endif
        {
            // Double-check inside critical section (other thread might have updated)
            if (std::abs(t - last_update_time_) > TIME_EPSILON) {
                solver_->update(*ctx.all_ions);
                last_update_time_ = t;
                solver_updated_this_step_ = true;
            }
        }
    }
    
    // Interpolate electric field at ion position
    Vec3 E_sc = solver_->fieldAt(ion.pos);
    
    // Compute force: F = q·E
    return E_sc * ion.ion_charge_C;
}

bool SpaceChargeGrid::applies_to(const IonState& ion) const {
    // Space charge applies to all charged particles
    // (Actual computation in compute() checks for active/born status)
    return true;
}

std::string SpaceChargeGrid::name() const {
    return "SpaceChargeGrid";
}

}  // namespace ICARION::physics

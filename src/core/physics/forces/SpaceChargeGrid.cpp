// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

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
    
    // Validate context (need ion ensemble for grid update)
    if (!ctx.all_ions && !ctx.ion_ensemble) {
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
                if (ctx.ion_ensemble) {
                    solver_->update(*ctx.ion_ensemble);
                } else {
                    solver_->update(*ctx.all_ions);
                }
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

Vec3 SpaceChargeGrid::compute_soa(
    const core::IonEnsemble& ensemble,
    size_t ion_idx,
    double t,
    const ForceContext& ctx
) const {
    ForceContext ctx_with = ctx;
    ctx_with.ion_ensemble = &ensemble;
    ctx_with.ion_index = ion_idx;

    // Use SoA update path
    if (!ctx_with.all_ions && !ctx_with.ion_ensemble) {
        ctx_with.ion_ensemble = &ensemble;
    }

    // Build minimal IonState for interpolation point
    IonState ion;
    ion.pos = ensemble.get_pos(ion_idx);
    ion.vel = ensemble.get_vel(ion_idx);
    ion.ion_charge_C = ensemble.charge_data()[ion_idx];
    ion.mass_kg = ensemble.mass_data()[ion_idx];
    ion.active = ensemble.active_data()[ion_idx] != 0;
    ion.born = ensemble.born_data()[ion_idx] != 0;
    ion.current_domain_index = ensemble.domain_index(ion_idx);
    ion.species_id = ensemble.species_id(ion_idx);
    ion.CCS_m2 = ensemble.CCS(ion_idx);
    ion.reduced_mobility_cm2_Vs = ensemble.mobility(ion_idx);
    ion.birth_time_s = ensemble.birth_time(ion_idx);

    return compute(ion, t, ctx_with);
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

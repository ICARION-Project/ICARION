// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifdef ICARION_USE_GPU

#include "SpaceChargeGPU.h"
#include "ForceContext.h"
#include <stdexcept>
#include <cmath>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace ICARION::physics {

SpaceChargeGPU::SpaceChargeGPU(std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> solver)
    : solver_(solver),
      last_update_time_(-1.0),
      solver_updated_this_step_(false)
{
    if (!solver_) {
        throw std::invalid_argument(
            "SpaceChargeGPU: Solver cannot be null. "
            "Pass a valid GPUSpaceChargeP3M instance."
        );
    }
}

Vec3 SpaceChargeGPU::compute(
    const IonState& ion,
    double t,
    const ForceContext& ctx
) const {
    // Validate solver
    if (!solver_) {
        return {0, 0, 0};
    }
    
    // Validate context (need ion ensemble for GPU update)
    if ((!ctx.all_ions || ctx.all_ions->empty()) && !ctx.ion_ensemble) {
        return {0, 0, 0};
    }
    
    const size_t N = ctx.ion_ensemble ? ctx.ion_ensemble->size() : ctx.all_ions->size();
    
    // Update solver once per timestep (time-based deduplication)
    // Use epsilon comparison for floating-point time
    constexpr double TIME_EPSILON = 1e-12;  // 1 ps tolerance
    
    if (std::abs(t - last_update_time_) > TIME_EPSILON) {
        // New timestep detected - update solver
        #ifdef _OPENMP
        #pragma omp critical(space_charge_gpu_solver_update)
        #endif
        {
            // Double-check inside critical section (other thread might have updated)
            if (std::abs(t - last_update_time_) > TIME_EPSILON) {
                // Resize cache if needed
                if (cached_E_fields_.size() != N) {
                    cached_E_fields_.resize(N);
                }
                
                // Call GPU solver: Computes E-field for ALL ions at once
                bool success = ctx.ion_ensemble
                    ? solver_->compute_space_charge_field(*ctx.ion_ensemble, cached_E_fields_)
                    : solver_->compute_space_charge_field(*ctx.all_ions, cached_E_fields_);
                
                if (!success) {
                    // GPU solver failed - zero out fields
                    std::fill(cached_E_fields_.begin(), cached_E_fields_.end(), Vec3{0, 0, 0});
                }
                
                last_update_time_ = t;
                solver_updated_this_step_ = true;
            }
        }
    }
    
    // Use provided ion_index when available (SoA path)
    size_t ion_idx = (ctx.ion_index != static_cast<size_t>(-1)) ? ctx.ion_index : 0;

    if (ctx.ion_index == static_cast<size_t>(-1) && ctx.all_ions) {
        // Linear search fallback (AoS path)
        for (; ion_idx < N; ++ion_idx) {
            if (&(*ctx.all_ions)[ion_idx] == &ion) {
                break;
            }
        }
    }

    if (ion_idx >= N || ion_idx >= cached_E_fields_.size()) {
        return {0, 0, 0};
    }

    // Retrieve precomputed E-field from cache
    Vec3 E_sc = cached_E_fields_[ion_idx];
    
    // Compute force: F = q·E
    return E_sc * ion.ion_charge_C;
}

Vec3 SpaceChargeGPU::compute_soa(
    const core::IonEnsemble& ensemble,
    size_t ion_idx,
    double t,
    const ForceContext& ctx
) const {
    if (!solver_) {
        return {0, 0, 0};
    }

    // Fast path: skip inactive or unborn ions
    if (ensemble.active_data()[ion_idx] == 0 || ensemble.born_data()[ion_idx] == 0) {
        return {0, 0, 0};
    }

    // Ensure we have fields for this timestep
    update_fields_if_needed(ensemble, t);

    if (ion_idx >= cached_E_fields_.size()) {
        return {0, 0, 0};
    }

    const Vec3 E_sc = cached_E_fields_[ion_idx];
    const double q = ensemble.charge_data()[ion_idx];
    return E_sc * q;
}

bool SpaceChargeGPU::applies_to(const IonState& ion) const {
    // Space charge applies to all charged particles
    // (Actual computation in compute() checks for active/born status)
    return true;
}

std::string SpaceChargeGPU::name() const {
    return "SpaceChargeGPU";
}

void SpaceChargeGPU::update_fields_if_needed(
    const core::IonEnsemble& ensemble,
    double t
) const {
    constexpr double TIME_EPSILON = 1e-12;

    if (std::abs(t - last_update_time_) <= TIME_EPSILON) {
        return;
    }

    const size_t N = ensemble.size();

    #ifdef _OPENMP
    #pragma omp critical(space_charge_gpu_solver_update)
    #endif
    {
        if (std::abs(t - last_update_time_) > TIME_EPSILON) {
            if (cached_E_fields_.size() != N) {
                cached_E_fields_.resize(N);
            }

            const bool success = solver_->compute_space_charge_field(ensemble, cached_E_fields_);
            if (!success) {
                std::fill(cached_E_fields_.begin(), cached_E_fields_.end(), Vec3{0, 0, 0});
            }

            last_update_time_ = t;
            solver_updated_this_step_ = true;
        }
    }
}

}  // namespace ICARION::physics

#endif  // ICARION_USE_GPU

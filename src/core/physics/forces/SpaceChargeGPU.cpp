// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifdef ICARION_USE_GPU

#include "SpaceChargeGPU.h"
#include "ForceContext.h"
#include <stdexcept>
#include <cmath>

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
    
    // Validate context (need all_ions for GPU update)
    if (!ctx.all_ions || ctx.all_ions->empty()) {
        return {0, 0, 0};
    }
    
    const auto& all_ions = *ctx.all_ions;
    const size_t N = all_ions.size();
    
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
                bool success = solver_->compute_space_charge_field(all_ions, cached_E_fields_);
                
                if (!success) {
                    // GPU solver failed - zero out fields
                    std::fill(cached_E_fields_.begin(), cached_E_fields_.end(), Vec3{0, 0, 0});
                }
                
                last_update_time_ = t;
                solver_updated_this_step_ = true;
            }
        }
    }
    
    // Find ion index in all_ions (linear search - assume small overhead vs GPU computation)
    // Alternative: Pass ion_index via ForceContext (requires API change)
    size_t ion_idx = 0;
    for (; ion_idx < N; ++ion_idx) {
        if (&all_ions[ion_idx] == &ion) {
            break;
        }
    }
    
    // Fallback: If ion not found, use zero field (should never happen)
    if (ion_idx >= N || ion_idx >= cached_E_fields_.size()) {
        return {0, 0, 0};
    }
    
    // Retrieve precomputed E-field from cache
    Vec3 E_sc = cached_E_fields_[ion_idx];
    
    // Compute force: F = q·E
    return E_sc * ion.ion_charge_C;
}

bool SpaceChargeGPU::applies_to(const IonState& ion) const {
    // Space charge applies to all charged particles
    // (Actual computation in compute() checks for active/born status)
    return true;
}

std::string SpaceChargeGPU::name() const {
    return "SpaceChargeGPU";
}

}  // namespace ICARION::physics

#endif  // ICARION_USE_GPU

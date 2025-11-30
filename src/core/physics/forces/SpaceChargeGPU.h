// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#pragma once

#ifdef ICARION_USE_GPU

#include "IForce.h"
#include "core/gpu/spacecharge/GPUSpaceChargeP3M.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include <memory>
#include <vector>

namespace ICARION::physics {

/**
 * @brief GPU-accelerated space charge force using P³M algorithm (O(N log N))
 * 
 * Integrates GPUSpaceChargeP3M into the IForce framework for automatic
 * GPU-accelerated space charge computation. Updates P³M solver once per 
 * timestep for all ions, then returns precomputed field for each ion.
 * 
 * **Performance vs CPU:**
 * - N = 1,000:   2 ms/timestep (GPU) vs 20 ms (CPU Grid) → **10× speedup**
 * - N = 10,000:  15 ms/timestep (GPU) vs 200 ms (CPU Grid) → **13× speedup**
 * - N = 100,000: 50 ms/timestep (GPU) vs 2000 ms (CPU Grid) → **40× speedup**
 * 
 * **When to use:**
 * - N ≥ 1000 ions (below this, CPU is competitive)
 * - GPU context available
 * - Space charge effects significant
 * 
 * **Algorithm:**
 * 1. Particle-to-Grid (P2G): Scatter charges to 3D grid using CIC interpolation
 * 2. FFT: Forward transform charge density to k-space (cuFFT)
 * 3. Poisson Solve: φ̂(k) = ρ̂(k) / (ε₀ k²) in spectral domain
 * 4. IFFT: Inverse transform potential to real space
 * 5. Gradient: Compute E = -∇φ on grid
 * 6. Grid-to-Particle (G2P): Interpolate E-field to ion positions
 * 
 * **Limitations:**
 * - Rectangular domain only (no cylindrical grids yet)
 * - Fixed grid resolution (configured at construction)
 * - ~20% accuracy error for near-field (ion separation < 50 cells)
 * - Assumes periodic boundaries (FFT limitation)
 * 
 * **Usage:**
 * @code
 * // Create GPU context
 * auto gpu_ctx = icarion::gpu::GPUContext::create(0);
 * 
 * // Configure P³M solver
 * icarion::gpu::GPUSpaceChargeP3M::Config config;
 * config.grid_nx = 128;
 * config.grid_ny = 128;
 * config.grid_nz = 128;
 * config.domain_min = Vec3{-0.05, -0.05, -0.05};  // 10cm box
 * config.domain_max = Vec3{0.05, 0.05, 0.05};
 * 
 * // Create solver
 * auto gpu_solver = icarion::gpu::GPUSpaceChargeP3M::create(gpu_ctx, config);
 * 
 * // Wrap in IForce interface
 * auto force = std::make_unique<SpaceChargeGPU>(gpu_solver);
 * force_registry.add_force(std::move(force));
 * @endcode
 * 
 * @note Requires ctx.all_ions to be populated in ForceContext
 * @note Updates solver only once per unique timestep (time-based deduplication)
 * @note For N < 1000, prefer SpaceChargeDirect (CPU is faster + exact)
 * 
 * @see GPUSpaceChargeP3M for algorithm details
 * @see SpaceChargeGrid for CPU fallback
 * @see SpaceChargeDirect for small N exact solution
 */
class SpaceChargeGPU : public IForce {
public:
    /**
     * @brief Construct GPU space charge force from P³M solver
     * 
     * @param solver Shared pointer to GPUSpaceChargeP3M (must not be null)
     * 
     * @throws std::invalid_argument if solver is null
     */
    explicit SpaceChargeGPU(std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> solver);
    
    // =========================================================================
    // IForce Interface
    // =========================================================================
    
    /**
     * @brief Compute space charge force from GPU P³M field
     * 
     * Updates P³M solver once per timestep (for all ions), then
     * returns precomputed electric field at the given ion position.
     * 
     * @param ion Current ion (compute force ON this ion)
     * @param t Simulation time [s] (used for update deduplication)
     * @param ctx Force context - MUST contain all_ions vector
     * @return Space charge force [N] = q·E_sc(r)
     * 
     * @note Returns {0,0,0} if ctx.all_ions is null or solver fails
     * @note Thread-safe: Only first thread updates solver (OpenMP critical section)
     * @note E-field is cached after update, so subsequent calls are O(1)
     */
    Vec3 compute(
        const IonState& ion,
        double t,
        const ForceContext& ctx
    ) const override;
    
    /**
     * @brief Check if force applies to this ion
     * 
     * Space charge applies to all ions (Coulomb interaction is universal).
     * Individual ion activity/birth checks are handled in compute().
     * 
     * @param ion Ion to check
     * @return Always true (all charged particles feel space charge)
     */
    bool applies_to(const IonState& ion) const override;
    
    /**
     * @brief Get force name for logging/debugging
     * @return "SpaceChargeGPU"
     */
    std::string name() const override;
    
    // =========================================================================
    // GPU-specific API
    // =========================================================================
    
    /**
     * @brief Get underlying P³M solver (for diagnostics/tuning)
     * @return Shared pointer to GPU solver
     */
    std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> get_solver() const { return solver_; }
    
    /**
     * @brief Check if solver was updated this timestep
     * @return True if solver.compute_space_charge_field() was called at current time
     */
    bool was_updated_this_step() const { return solver_updated_this_step_; }

private:
    std::shared_ptr<icarion::gpu::GPUSpaceChargeP3M> solver_;  ///< GPU P³M solver
    
    // Timestep tracking (for once-per-timestep updates)
    mutable double last_update_time_;        ///< Last time solver was updated [s]
    mutable bool solver_updated_this_step_;  ///< Flag: was solver updated this timestep?
    
    // Cached E-fields (updated once per timestep, indexed by ion position in all_ions)
    mutable std::vector<Vec3> cached_E_fields_;  ///< E-field at each ion position [V/m]
};

}  // namespace ICARION::physics

#endif  // ICARION_USE_GPU

// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file GPUCollisionHelper.h
 * @brief GPU batch collision processing
 * 
 * High-level interface for GPU collision handling.
 * Manages memory transfers, cuRAND states, and kernel launches.
 * 
 * **Design:**
 * - Automatic dispatch: GPU if N >= threshold, else return false
 * - Async pipeline: upload → compute → download
 * - Buffer management via GPUMemoryPool
 * - cuRAND state persistence across timesteps
 * - Performance statistics tracking
 * 
 * **Integration with SimulationEngine:**
 * ```cpp
 * if (gpu_collision_helper_ && ions.size() >= threshold) {
 *     if (gpu_collision_helper_->process_collisions_batch(ions, dt, env)) {
 *         return;  // GPU success
 *     }
 * }
 * // Fallback to CPU
 * for (auto& ion : ions) {
 *     collision_handler->handle_collision(ion, dt, rng, env);
 * }
 * ```
 */

#pragma once

#ifdef ICARION_USE_GPU

#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include "core/config/types/EnvironmentConfig.h"
#include <vector>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>

// Forward declarations
namespace ICARION { namespace gpu {
    class GPUContext;
    class GPUMemoryPool;
}}

// cuRAND forward declarations
struct curandStateXORWOW;
typedef struct curandStateXORWOW curandState;

namespace ICARION {
namespace gpu {

// Forward declarations from collision_kernels_gpu.cuh
struct EnvironmentParams_GPU;
struct GeometryData_GPU;

/**
 * @brief Performance statistics for GPU collision processing
 */
struct CollisionStats_GPU {
    size_t total_batches = 0;           ///< Total batches processed
    size_t total_ions_processed = 0;    ///< Total ions processed
    double total_gpu_time_ms = 0.0;     ///< Total GPU kernel time [ms]
    double total_transfer_time_ms = 0.0;///< Total H↔D transfer time [ms]
    double average_speedup = 0.0;       ///< vs CPU (estimated)
};

/**
 * @brief GPU collision batch processor
 * 
 * Manages GPU resources for batch collision processing.
 * Supports both HSS (isotropic) and EHSS (geometry-resolved) models.
 */
class GPUCollisionHelper {
public:
    /**
     * @brief Factory method for creating GPU collision helper
     * 
     * @param context GPU context (must outlive helper)
     * @param threshold Minimum ion count for GPU dispatch (default: 5000)
     * @param collision_model "HSS" or "EHSS"
     * @param rng_seed Random seed for cuRAND initialization
     * 
     * @return Unique pointer to helper, or nullptr if GPU unavailable
     */
    static std::unique_ptr<GPUCollisionHelper> create(
        const GPUContext& context,
        size_t threshold = 5000,
        const std::string& collision_model = "HSS",
        unsigned long long rng_seed = 42
    );
    
    /**
     * @brief Process collisions for batch of ions on GPU
     * 
     * @param ions Ion ensemble (modified in-place)
     * @param dt Timestep [s]
     * @param env Environment configuration
     * 
     * @return true if GPU processing succeeded, false if should fallback to CPU
     * 
     * **Algorithm:**
     * 1. Check if N >= threshold
     * 2. Upload ion velocities + masses + CCS to GPU
     * 3. Launch collision kernel (HSS or EHSS)
     * 4. Download modified velocities
     * 5. Update statistics
     * 
     * **Thread Safety:**
     * - Not thread-safe (uses internal GPU stream)
     * - Do not call from multiple threads
     */
    bool process_collisions_batch(
        std::vector<IonState>& ions,
        double dt,
        const config::EnvironmentConfig& env
    );
    
    /**
     * @brief Type alias for geometry map
     */
    using GeometryMap = std::unordered_map<std::string, std::pair<std::vector<Vec3>, std::vector<double>>>;
    
    /**
     * @brief Set molecular geometry for EHSS model
     * 
     * Must be called before processing if using EHSS model.
     * 
     * @param geometry_map Species → (atom_centers, atom_radii)
     */
    void set_geometry(const GeometryMap& geometry_map);
    
    /**
     * @brief Get performance statistics
     */
    const CollisionStats_GPU& get_stats() const { return stats_; }
    
    /**
     * @brief Reset performance statistics
     */
    void reset_stats() { stats_ = {}; }
    
    /**
     * @brief Get dispatch threshold
     */
    size_t get_threshold() const { return threshold_; }
    
    ~GPUCollisionHelper();

private:
    GPUCollisionHelper(
        const GPUContext& context,
        size_t threshold,
        const std::string& collision_model,
        unsigned long long rng_seed
    );
    
    // Disable copy/move
    GPUCollisionHelper(const GPUCollisionHelper&) = delete;
    GPUCollisionHelper& operator=(const GPUCollisionHelper&) = delete;
    
    /**
     * @brief Initialize cuRAND states (called once)
     */
    void initialize_curand_states(size_t n_states);
    
    /**
     * @brief Upload geometry to GPU (for EHSS)
     */
    void upload_geometry_to_gpu();
    
    /**
     * @brief Convert EnvironmentConfig to GPU format
     */
    EnvironmentParams_GPU convert_environment_params(
        const config::EnvironmentConfig& env
    ) const;
    
    const GPUContext& context_;
    size_t threshold_;
    std::string collision_model_;  // "HSS" or "EHSS"
    unsigned long long rng_seed_;
    
    // cuRAND states (persistent across timesteps)
    curandState* d_curand_states_ = nullptr;
    size_t n_curand_states_ = 0;
    bool curand_initialized_ = false;
    
    // Geometry data (for EHSS)
    GeometryData_GPU geometry_gpu_;
    bool geometry_uploaded_ = false;
    const GeometryMap* geometry_map_host_ = nullptr;
    
    // Performance statistics
    CollisionStats_GPU stats_;
};

} // namespace gpu
} // namespace ICARION

#endif // ICARION_USE_GPU

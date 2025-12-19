// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file GPUCollisionHelper.h
 * @brief GPU batch collision processing
 * 
 * High-level interface for GPU collision handling.
 * Manages memory transfers, cuRAND states, and kernel launches. Experimental
 * and currently only lightly exercised; EHSS geometry mapping assumes a single
 * species index unless the caller populates `set_geometry` appropriately.
 * 
 * **Design:**
 * - Automatic dispatch: GPU if N >= threshold, else return false
 * - Async pipeline: upload -> compute -> download
 * - Persistent buffers and cuRAND states reused across timesteps
 * - Minimal statistics tracking (batch/ion counters, timing placeholders)
 * 
 * **Integration with SimulationEngine (not wired by default):**
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
#include <unordered_set>
#include <utility>
#include <cstdint>  // For uint8_t

// Forward declarations (or full includes for non-CUDA files)
#ifndef __CUDACC__
    // Regular C++ compilation - include full definition
    #include "core/gpu/core/GPUContext.h"
#else
    // CUDA compilation - forward declare only (to avoid header issues)
    namespace icarion {
    namespace gpu {
        class GPUContext;
    }
    }
#endif

namespace icarion { 
namespace gpu {
    struct EnvironmentParams_GPU;  // Forward declaration (defined in collision_kernels_gpu.cuh)
    struct GeometryData_GPU;  // Forward declaration (defined in collision_kernels_gpu.cuh)
}
}

// cuRAND forward declarations
struct curandStateXORWOW;
typedef struct curandStateXORWOW curandState;

namespace icarion {
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
 * Supports both HSS (isotropic) and EHSS (geometry-resolved) models; EHSS
 * requires geometry upload and species-index mapping by the caller.
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
        const gpu::GPUContext& context,
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
        const ICARION::config::EnvironmentConfig& env
    );
    
    /**
     * @brief Type alias for geometry map
     */
    using GeometryMap = std::unordered_map<std::string, std::pair<std::vector<Vec3>, std::vector<double>>>;
    
    /**
     * @brief Set molecular geometry for EHSS model
     * 
     * Must be called before processing if using EHSS model. Geometry order
     * defines species indices expected in `IonState` entries supplied by the
     * caller; helper does not resolve IDs to indices.
     * 
     * @param geometry_map Species -> (atom_centers, atom_radii)
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
        const gpu::GPUContext& context,
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
    
    // Internal helper (implemented in .cu)
    // EnvironmentParams_GPU convert_environment_params(const ICARION::config::EnvironmentConfig&) const;
    
    void* cuda_stream_;  // cudaStream_t stored as void* to avoid nvcc header issues
    size_t threshold_;
    std::string collision_model_;  // "HSS" or "EHSS"
    unsigned long long rng_seed_;
    bool warned_mixture_limit_ = false;
    
    // cuRAND states (persistent across timesteps)
    curandState* d_curand_states_ = nullptr;
    size_t n_curand_states_ = 0;
    bool curand_initialized_ = false;
    
    // Persistent GPU buffers (reused across timesteps to avoid malloc/free overhead)
    double *d_vx_ = nullptr;
    double *d_vy_ = nullptr;
    double *d_vz_ = nullptr;
    double *d_mass_ = nullptr;
    double *d_ccs_ = nullptr;
    uint8_t *d_active_ = nullptr;
    int *d_species_indices_ = nullptr;
    size_t buffer_capacity_ = 0;  // Current buffer size in # of ions
    
    // Geometry data (for EHSS)
    GeometryData_GPU* geometry_gpu_ = nullptr;
    bool geometry_uploaded_ = false;
    const GeometryMap* geometry_map_host_ = nullptr;
    std::unordered_map<std::string, int> species_index_map_;  ///< species_id -> geometry index
    std::unordered_set<std::string> missing_species_warned_;
    
    // Performance statistics
    CollisionStats_GPU stats_;
};

} // namespace gpu
} // namespace icarion

#endif // ICARION_USE_GPU

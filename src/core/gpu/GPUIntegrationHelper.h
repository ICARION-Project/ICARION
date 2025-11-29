#ifndef ICARION_GPU_INTEGRATION_HELPER_H
#define ICARION_GPU_INTEGRATION_HELPER_H

#include "core/gpu/GPUContext.h"
#include "core/gpu/GPUMemoryPool.h"
#include "core/types/IonState.h"
#include "utils/IonState_GPU.h"
#include <vector>
#include <memory>

namespace icarion {
namespace gpu {

/**
 * @brief GPU batch integration helper for SimulationEngine
 * 
 * Provides high-performance batch RK4 integration on GPU when ion count
 * exceeds threshold. Designed to be called from SimulationEngine's parallel
 * ion loop as an optimization path.
 * 
 * Usage in SimulationEngine:
 *   if (gpu_helper && ions.size() >= GPU_THRESHOLD) {
 *       gpu_helper->integrate_batch_rk4(ions, dt, t);
 *   } else {
 *       // CPU integration per-ion
 *   }
 */
class GPUIntegrationHelper {
public:
    /**
     * @brief Create GPU integration helper
     * 
     * @param context GPU context (must outlive this object)
     * @param threshold Minimum ion count for GPU dispatch (default: 5000)
     * @return Unique pointer, or nullptr if GPU unavailable
     */
    static std::unique_ptr<GPUIntegrationHelper> create(
        const GPUContext& context,
        size_t threshold = 5000
    );
    
    ~GPUIntegrationHelper();
    
    /**
     * @brief Batch RK4 integration on GPU
     * 
     * Integrates all ions from t to t+dt using GPU acceleration.
     * Falls back to returning false if GPU error occurs.
     * 
     * @param ions Ion states (modified in-place)
     * @param dt Timestep [s]
     * @param t Current time [s]
     * @return true if GPU integration succeeded, false if fallback needed
     * 
     * @note Currently uses zero fields (Phase 3 will add field interpolation)
     */
    bool integrate_batch_rk4(
        std::vector<ICARION::core::IonState>& ions,
        double dt,
        double t
    );
    
    /**
     * @brief Check if GPU is enabled and operational
     */
    bool is_enabled() const { return enabled_; }
    
    /**
     * @brief Get GPU threshold
     */
    size_t get_threshold() const { return threshold_; }
    
    /**
     * @brief Set GPU threshold dynamically
     */
    void set_threshold(size_t new_threshold) { threshold_ = new_threshold; }
    
    /**
     * @brief Get performance statistics
     */
    struct Stats {
        size_t gpu_integrations;     ///< Number of GPU batch integrations
        size_t total_ions_gpu;        ///< Total ions processed on GPU
        double total_time_ms;         ///< Total GPU time [ms]
    };
    
    const Stats& get_stats() const { return stats_; }
    void reset_stats() { stats_ = {}; }

private:
    GPUIntegrationHelper(
        const GPUContext& context,
        std::unique_ptr<GPUMemoryPool> memory_pool,
        size_t threshold
    );
    
    const GPUContext& context_;
    std::unique_ptr<GPUMemoryPool> memory_pool_;
    
    // GPU buffers (reused)
    IonStateGPU ions_gpu_in_;
    IonStateGPU ions_gpu_out_;
    size_t allocated_capacity_;
    
    size_t threshold_;
    bool enabled_;
    Stats stats_;
};

} // namespace gpu
} // namespace icarion

#endif // ICARION_GPU_INTEGRATION_HELPER_H

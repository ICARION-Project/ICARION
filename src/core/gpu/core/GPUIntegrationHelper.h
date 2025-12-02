// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifndef ICARION_GPU_INTEGRATION_HELPER_H
#define ICARION_GPU_INTEGRATION_HELPER_H

#include "core/gpu/core/GPUContext.h"
#include "core/gpu/core/GPUMemoryPool.h"
#include "core/gpu/fields/FieldArrayGPU.h"
#include "core/types/IonState.h"
#include "utils/IonState_GPU.h"
#include <vector>
#include <memory>

// Forward declarations
class IFieldProvider;
struct FieldArray;

namespace icarion {
namespace gpu {

/**
 * @brief GPU batch integration helper for SimulationEngine
 * 
 * Provides batch RK4/RK45/Boris integration on GPU when ion count exceeds
 * threshold. Experimental and not wired into SimulationEngine by default;
 * callers must handle GPU/CPU selection and field upload.
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
     * Falls back to returning false if GPU error occurs. Not validated against
     * the CPU path and not enabled in production by default.
     * 
     * @param ions Ion states (modified in-place)
     * @param dt Timestep [s]
     * @param t Current time [s]
     * @param field_provider Optional field provider for position-dependent fields
     * @return true if GPU integration succeeded, false if fallback needed
     * 
     * If field_provider is provided and implements GridFieldProvider,
     * fields will be uploaded to GPU for hardware-accelerated interpolation.
     * Otherwise uses zero fields (free-particle motion).
     */
    bool integrate_batch_rk4(
        std::vector<ICARION::core::IonState>& ions,
        double dt,
        double t,
        const IFieldProvider* field_provider = nullptr
    );
    
    /**
     * @brief Batch RK45 (Dormand-Prince) integration on GPU
     * 
     * Adaptive 4th/5th order Runge-Kutta with automatic step size control.
     * Each ion adapts its substep size independently based on local error.
     * Experimental; expect to fall back to CPU if instability is observed.
     * 
     * @param ions Ion states (modified in-place)
     * @param dt Maximum timestep [s]
     * @param t Current time [s]
     * @param field_provider Optional field provider for position-dependent fields
     * @param atol Absolute tolerance (default: 1e-12)
     * @param rtol Relative tolerance (default: 1e-6)
     * @return true if GPU integration succeeded, false if fallback needed
     */
    bool integrate_batch_rk45(
        std::vector<ICARION::core::IonState>& ions,
        double dt,
        double t,
        const IFieldProvider* field_provider = nullptr,
        double atol = 1e-12,
        double rtol = 1e-6
    );
    
    /**
     * @brief Batch Boris pusher integration on GPU
     * 
     * Symplectic method optimal for magnetic field-dominated systems.
     * Experimental and not connected to the main integration flow.
     * 
     * @param ions Ion states (modified in-place)
     * @param dt Timestep [s]
     * @param t Current time [s]
     * @param field_provider Optional field provider for position-dependent fields
     * @return true if GPU integration succeeded, false if fallback needed
     */
    bool integrate_batch_boris(
        std::vector<ICARION::core::IonState>& ions,
        double dt,
        double t,
        const IFieldProvider* field_provider = nullptr
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
     * @brief Batch boundary checking on GPU (cylindrical geometry)
     * 
     * Checks all ions against cylindrical domain boundaries in parallel.
     * Deactivates ions that are outside bounds.
     * 
     * @param ions Ion states (active flags modified in-place)
     * @param length_m Domain length [m]
     * @param radius_m Domain radius [m]
     * @param is_last_domain True if last domain (strict z check)
     * @return true if GPU check succeeded, false if fallback needed
     * 
     * Boundary checks:
     * - Axial: -ε ≤ z < length (last) or z ≤ length+ε (transition)
     * - Radial: sqrt(x² + y²) ≤ radius + ε
     * 
     * @note Orbitrap hyperlogarithmic boundaries not supported (handled on CPU)
     */
    bool check_boundaries_batch(
        std::vector<ICARION::core::IonState>& ions,
        double length_m,
        double radius_m,
        bool is_last_domain
    );
    
    /**
     * @brief Get performance statistics
     */
    struct Stats {
        size_t gpu_integrations;     ///< Number of GPU batch integrations
        size_t gpu_boundary_checks;  ///< Number of GPU boundary checks
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
    
    // GPU field storage (managed)
    FieldArrayGPU field_array_gpu_;
    bool has_gpu_fields_;
    
    size_t threshold_;
    bool enabled_;
    Stats stats_;
    
    // Helper: Try to extract FieldArray from IFieldProvider
    const FieldArray* try_extract_field_array(const IFieldProvider* provider) const;
};

} // namespace gpu
} // namespace icarion

#endif // ICARION_GPU_INTEGRATION_HELPER_H

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * GPU-Accelerated Fast Multipole Method for ICARION
 * 
 * Reduces space charge field calculation complexity from O(N²) to O(N log N)
 * using the Fast Multipole Method with GPU acceleration.
 * 
 */

#pragma once

#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif

#include <vector>
#include <memory>
#include <chrono>
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

namespace ICARION {
namespace gpu {

// Forward declaration
struct FMMNode;

/**
 * @brief GPU-accelerated Fast Multipole Method for space charge calculations
 * 
 * This class implements the Fast Multipole Method (FMM) to compute electrostatic
 * interactions between charged particles with O(N log N) complexity instead of
 * the traditional O(N²) direct summation.
 * 
 * Key features:
 * - O(N log N) complexity vs O(N²) direct summation
 * - GPU acceleration for tree construction and field evaluation
 * - Adaptive precision control for speed vs accuracy trade-offs
 * - Memory pool integration for efficient GPU memory management
 */
class GPUFastMultipoleMethod {
public:
    /**
     * @brief Constructor
     * @param max_ions Maximum number of ions the system will handle
     * @param precision FMM precision parameter (higher = more accurate, slower)
     * @param use_gpu Enable GPU acceleration for FMM operations
     */
    GPUFastMultipoleMethod(int max_ions = 100000, int precision = 8, bool use_gpu = true);
    
    /**
     * @brief Destructor - cleanup GPU resources
     */
    ~GPUFastMultipoleMethod();

    /**
     * @brief Initialize the FMM system with current simulation parameters
     * @param ions Reference ions to establish spatial bounds
     */
    void initialize(const std::vector<core::IonState>& ions);

    /**
     * @brief Compute space charge electric field for all ions using FMM
     * @param ions Vector of all ions in the simulation
     * @param eps0 Permittivity of free space
     * @return Vector of electric field contributions for each ion [V/m]
     */
    std::vector<core::Vec3> computeSpaceChargeFields(const std::vector<core::IonState>& ions, 
                                                     double eps0);

    /**
     * @brief Compute space charge field for a single ion using FMM
     * @param ion Target ion for field calculation
     * @param all_ions Vector of all ions in the simulation
     * @param eps0 Permittivity of free space
     * @return Electric field at ion position [V/m]
     */
    core::Vec3 computeFieldAtIon(const core::IonState& ion,
                                 const std::vector<core::IonState>& all_ions,
                                 double eps0);

    /**
     * @brief Get performance statistics from last FMM calculation
     * @return Struct containing timing and efficiency metrics
     */
    struct PerformanceStats {
        double tree_build_time_ms;      // Time to build octree [ms]
        double fmm_evaluation_time_ms;  // Time for FMM evaluation [ms]
        double total_time_ms;           // Total computation time [ms]
        double speedup_vs_direct;       // Speedup compared to O(N²) method
        int num_particles;              // Number of particles processed
        double flops_per_second;        // Computational performance [GFLOPS]
        bool using_gpu_acceleration;   // Whether GPU was used
    };
    
    PerformanceStats getPerformanceStats() const { return perf_stats_; }

    /**
     * @brief Enable/disable GPU acceleration during runtime
     * @param enable True to enable GPU acceleration
     */
    void setGPUAcceleration(bool enable);

    /**
     * @brief Adjust FMM precision for speed vs accuracy trade-off
     * @param precision Higher values = more accurate but slower (typical: 6-12)
     */
    void setPrecision(int precision);

    /**
     * @brief Get current FMM configuration
     */
    struct Configuration {
        int precision;           // FMM expansion order
        int max_particles;       // Maximum particles supported
        bool gpu_acceleration;   // Whether GPU acceleration is enabled
        int ncrit;              // Max particles per leaf node
        double relative_error;   // Target relative accuracy
    };
    
    Configuration getConfiguration() const { return config_; }

private:
    Configuration config_;
    PerformanceStats perf_stats_;
    
    // FMM internal state
    bool initialized_;
    double domain_size_;     // Size of computational domain
    core::Vec3 domain_center_;  // Center of computational domain
    
#ifdef USE_CUDA
    // GPU resources
    cudaStream_t fmm_stream_;
    void* d_positions_;      // GPU positions buffer
    void* d_charges_;        // GPU charges buffer
    void* d_fields_;         // GPU fields buffer
    size_t gpu_buffer_size_;
#endif

    // FMM library interface
    void* fmm_instance_;     // Opaque pointer to FMM library instance
    
    // Performance monitoring
    std::chrono::steady_clock::time_point timer_start_;
    void startTimer();
    double stopTimer(); // Returns elapsed time in ms
    
    // Helper methods
    void updateDomainBounds(const std::vector<core::IonState>& ions);
    void prepareFMMData(const std::vector<core::IonState>& ions);
    std::vector<core::Vec3> extractResults(const std::vector<core::IonState>& ions);
    std::vector<core::Vec3> computeDirectSummation(const std::vector<core::IonState>& ions, double eps0);
    
    // FMM tree methods
    void buildOctree(const std::vector<core::IonState>& ions);
    void subdivideNode(FMMNode* node, const std::vector<core::IonState>& ions, int depth);
    void computeMultipoleMoments(const std::vector<core::IonState>& ions);
    void computeMultipoleMomentsRecursive(FMMNode* node, const std::vector<core::IonState>& ions);
    std::vector<core::Vec3> computeFieldContributions(const std::vector<core::IonState>& ions, double eps0);
    core::Vec3 computeFieldFromTree(FMMNode* node, const core::IonState& target_ion, 
                                   const std::vector<core::IonState>& ions, double k, double theta);
    
    // GPU-specific helpers
#ifdef USE_CUDA
    void allocateGPUBuffers(int num_ions);
    void freeGPUBuffers();
    void uploadToGPU(const std::vector<core::IonState>& ions);
    void downloadFromGPU(std::vector<core::Vec3>& fields);
#endif
};

/**
 * @brief Factory function to create optimal FMM instance based on problem size
 * @param expected_ion_count Expected number of ions in simulation
 * @param target_accuracy Target relative accuracy (e.g., 1e-6)
 * @param prefer_speed If true, optimize for speed over accuracy
 * @return Configured FMM instance
 */
std::unique_ptr<GPUFastMultipoleMethod> createOptimalFMM(int expected_ion_count,
                                                         double target_accuracy = 1e-6,
                                                         bool prefer_speed = false);

} // namespace gpu
} // namespace ICARION
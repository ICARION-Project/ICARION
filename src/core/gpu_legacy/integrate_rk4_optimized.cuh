/**
 * =====================================================================
 *   ICARION: Optimized RK4 Kernel Header
 * =====================================================================
 *   @file        integrate_rk4_optimized.cuh
 *   @brief       Header for optimized RK4 integration kernel
 *
 *   @date        2025-11-11
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 * =====================================================================
 */

#pragma once

#ifdef __cplusplus
extern "C" {
#endif

#include "utils/IonState_GPU.h"
#include "paramUtils/GPUParams.h"
#include <curand_kernel.h>

/**
 * @brief Optimized RK4 integration step with shared memory caching
 * 
 * Performance improvements over baseline:
 * - Shared memory for domain parameters
 * - Fast math operations
 * - Reduced global memory traffic
 * - Loop unrolling
 * 
 * @param d_ions Device pointer to ion states
 * @param n Number of ions
 * @param g Global parameters
 * @param doms Device pointer to domain array
 * @param n_domains Number of domains
 * @param t_start Starting time
 * @param dt_local Time step
 * @param t_target Target time
 * @param d_rng_states Device pointer to RNG states (can be nullptr)
 */
void integrate_rk4_step_gpu_optimized(
    IonStateGPU* d_ions,
    int n,
    const GlobalParamsGPU& g,
    const DomainGPU* doms,
    int n_domains,
    double t_start,
    double dt_local,
    double t_target,
    curandState* d_rng_states
);

/**
 * @brief Get theoretical occupancy of optimized kernel
 * @return Occupancy percentage (0.0 - 1.0)
 */
float get_optimized_kernel_occupancy();

#ifdef __cplusplus
}
#endif

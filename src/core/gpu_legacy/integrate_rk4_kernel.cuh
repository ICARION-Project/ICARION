/**
 * =====================================================================
 *  _______ _____            _____ ______
 * |__   __|  __ \     /\   / ____|  ____|
 *    | |  | |__) |   /  \ | |    | |__
 *    | |  |  _  /   / /\ \| |    |  __|
 *    | |  | | \ \  / ____ \ |____| |____
 *    |_|  |_|  \_\/_/    \_\_____|______|
 *
 *   ICARION: A Modular Framework for Ion Collision and Reaction Integration
 * =====================================================================
 *   A C++ framework for simulating ion trajectories in electric fields
 *   under the influence of gas.
 *
 *   @file        integrate_rk4_kernel.cuh
 *   @brief       GPU Kernel Header: RK4 Integration Step
 *
 *   @details
 *   Header file declaring the GPU RK4 integration kernel functions.
 *   This file can be included in C++ source files.
 *
 *   @date        <2025-10-17>
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#pragma once

#include "utils/IonState_GPU.h"
#include "paramUtils/GPUParams.h"
#include "gpuUtils/gpu_geometry.h"
#include <curand_kernel.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief GPU wrapper function for RK4 integration step
 * 
 * @param d_ions Pointer to ion states on GPU device memory
 * @param n Number of ions
 * @param g Global parameters for GPU
 * @param dom Domain parameters for GPU
 * @param t_local Current time
 * @param dt_local Time step
 * @param t_target Target time
 * @param d_rng_states Device pointer to per-ion cuRAND states (can be nullptr if no stochastic collisions)
 * @param d_multi_geom Device pointer to multi-species geometry data (can be nullptr if not using EHSS)
 */
void integrate_rk4_step_gpu(IonStateGPU* d_ions, int n,
                           const GlobalParamsGPU& g, const DomainGPU* doms, int n_domains,
                           double t_local, double dt_local, double t_target,
                           curandState* d_rng_states,
                           const MultiSpeciesGeometryGPU* d_multi_geom);

/**
 * @brief Count number of active ions on GPU
 * 
 * @param d_ions Pointer to ion states on GPU device memory
 * @param n Number of ions
 * @return Number of active ions
 */
int count_active_ions_gpu(IonStateGPU* d_ions, int n);

#ifdef __cplusplus
}
#endif
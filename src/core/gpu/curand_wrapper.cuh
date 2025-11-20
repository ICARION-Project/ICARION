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
 *   @file        curand_wrapper.cuh
 *   @brief       GPU Random Number Generation Wrapper
 *
 *   @details
 *   Wrapper for cuRAND library to provide RNG state management for
 *   stochastic collision models (EHSS, HSMC) on GPU.
 *
 *   Each ion gets its own persistent RNG state that is maintained
 *   across timesteps, allowing for reproducible simulations with
 *   fixed seeds.
 *
 *   @date        October 21, 2025
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#pragma once

#include <curand_kernel.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * @brief Initialize cuRAND RNG states for all ions
 * 
 * Allocates and initializes one RNG state per ion. Each state gets a
 * unique sequence number based on its ion index, ensuring different
 * random sequences while maintaining reproducibility with fixed seed.
 * 
 * @param d_rng_states Pointer to device memory for RNG states (must be allocated)
 * @param n_ions Number of ions (and RNG states)
 * @param seed Random seed for reproducibility (0 = time-based)
 * 
 * @note RNG states should persist for entire simulation (don't free until end)
 * @note Each state is ~48 bytes, so memory overhead is minimal
 */
void init_rng_states_gpu(curandState* d_rng_states, int n_ions, unsigned long long seed);

/**
 * @brief Test RNG generation (for validation)
 * 
 * Generates random numbers and computes basic statistics to validate
 * RNG is working correctly. Used for unit testing.
 * 
 * @param d_rng_states Pointer to initialized RNG states
 * @param n_ions Number of ions
 * @param d_output Pointer to device memory for output (size = n_ions)
 */
void test_rng_generation_gpu(curandState* d_rng_states, int n_ions, float* d_output);

#ifdef __cplusplus
}
#endif

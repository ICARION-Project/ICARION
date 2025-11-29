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
 *   @file        curand_wrapper.cu
 *   @brief       GPU Random Number Generation Implementation
 *
 *   @details
 *   Implementation of cuRAND wrapper for stochastic collision models.
 *   Uses XORWOW algorithm (cuRAND default) for high-quality RNG.
 *
 *   @date        October 21, 2025
 *   @version     1.0
 *   @author      Christoph Schäfer
 *   @license     MIT License
 *
 * =====================================================================
 */

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <stdio.h>
#include <time.h>

/**
 * @brief Kernel to initialize RNG states
 * 
 * Each thread initializes one RNG state with:
 * - Common seed (for reproducibility)
 * - Unique sequence number (based on ion index)
 * - Offset = 0 (starting point in sequence)
 */
__global__ void init_rng_kernel(curandState* states, int n, unsigned long long seed) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    // Initialize RNG state for this ion
    // seed: same for all (reproducibility)
    // sequence: unique per ion (different random sequences)
    // offset: 0 (start of sequence)
    curand_init(seed, idx, 0, &states[idx]);
}

/**
 * @brief Kernel to test RNG generation
 * 
 * Generates one random number per ion for validation.
 */
__global__ void test_rng_kernel(curandState* states, int n, float* output) {
    int idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n) return;
    
    // Generate uniform random number [0, 1)
    curandState local_state = states[idx];
    float rand = curand_uniform(&local_state);
    output[idx] = rand;
    
    // Save state for next call (important!)
    states[idx] = local_state;
}

// ============================================================================
// Host wrapper functions
// ============================================================================

extern "C" void init_rng_states_gpu(curandState* d_rng_states, int n_ions, unsigned long long seed) {
    // If seed is 0, use current time for non-deterministic behavior
    if (seed == 0) {
        seed = (unsigned long long)time(NULL);
    }
    
    const int threads = 256;
    const int blocks = (n_ions + threads - 1) / threads;
    
    init_rng_kernel<<<blocks, threads>>>(d_rng_states, n_ions, seed);
    
    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        fprintf(stderr, "[cuRAND] RNG initialization failed: %s\n", 
                cudaGetErrorString(err));
    } else {
        printf("cuRAND: Initialized %d RNG states with seed %llu\n", n_ions, seed);
    }
}

extern "C" void test_rng_generation_gpu(curandState* d_rng_states, int n_ions, float* d_output) {
    const int threads = 256;
    const int blocks = (n_ions + threads - 1) / threads;
    
    test_rng_kernel<<<blocks, threads>>>(d_rng_states, n_ions, d_output);
    
    cudaError_t err = cudaDeviceSynchronize();
    if (err != cudaSuccess) {
        fprintf(stderr, "[cuRAND] RNG test generation failed: %s\n", 
                cudaGetErrorString(err));
    }
}

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once

#ifdef ICARION_USE_GPU

#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <cstdint>

namespace icarion {
namespace gpu {

struct DeviceReaction {
    int reactant_idx;
    int product_idx;
    double k_eff;  // [1/s] precomputed for the domain
};

/**
 * @brief Launch GPU reaction sampling.
 *
 * @param species_inout Device array of species indices (updated in-place)
 * @param domain_indices Device array of per-ion domain indices
 * @param active Device array of active flags
 * @param born Device array of born flags
 * @param n_ions Number of ions
 * @param reactions Device array of reactions (flattened across domains)
 * @param reaction_offsets Device array size (n_domains+1) with prefix offsets
 * @param n_domains Number of domains
 * @param rng_states Device array of curandStateXORWOW (one per ion)
 * @param dt Timestep [s]
 * @param stream CUDA stream to launch on
 */
void launch_reaction_kernel(uint32_t* species_inout,
                            const int32_t* domain_indices,
                            const uint8_t* active,
                            const uint8_t* born,
                            size_t n_ions,
                            const DeviceReaction* reactions,
                            const int32_t* reaction_offsets,
                            int n_domains,
                            curandStateXORWOW* rng_states,
                            double dt,
                            cudaStream_t stream);

/**
 * @brief Initialize XORWOW RNG states.
 */
void init_rng_states_xorwow(curandStateXORWOW* states,
                            size_t n_states,
                            unsigned long long seed,
                            cudaStream_t stream);

}  // namespace gpu
}  // namespace icarion

#endif  // ICARION_USE_GPU

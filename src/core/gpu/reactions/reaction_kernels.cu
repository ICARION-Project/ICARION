// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifdef ICARION_USE_GPU

#include "reaction_kernels.cuh"
#include <cuda_runtime.h>
#include <curand_kernel.h>
#include <math.h>

namespace icarion {
namespace gpu {

__global__ void reaction_kernel(uint32_t* species_inout,
                                const int32_t* domain_indices,
                                const uint8_t* active,
                                const uint8_t* born,
                                size_t n_ions,
                                const DeviceReaction* reactions,
                                const int32_t* reaction_offsets,
                                int n_domains,
                                curandStateXORWOW* rng_states,
                                double dt) {
    const size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx >= n_ions) return;
    if (!active[idx] || !born[idx]) return;

    const int dom = domain_indices[idx];
    if (dom < 0 || dom >= n_domains) return;

    const int32_t start = reaction_offsets[dom];
    const int32_t end = reaction_offsets[dom + 1];
    if (start == end) return;

    uint32_t species_idx = species_inout[idx];
    curandStateXORWOW local_state = rng_states[idx];

    for (int32_t r = start; r < end; ++r) {
        const DeviceReaction rxn = reactions[r];
        if (rxn.reactant_idx != static_cast<int>(species_idx)) {
            continue;
        }

        double p = 1.0 - exp(-rxn.k_eff * dt);
        float u = curand_uniform(&local_state);
        if (u < p) {
            species_idx = static_cast<uint32_t>(rxn.product_idx);
            break;  // apply first matching reaction
        }
    }

    rng_states[idx] = local_state;
    species_inout[idx] = species_idx;
}

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
                            cudaStream_t stream) {
    if (n_ions == 0 || n_domains == 0) return;
    constexpr int THREADS = 256;
    int blocks = static_cast<int>((n_ions + THREADS - 1) / THREADS);
    reaction_kernel<<<blocks, THREADS, 0, stream>>>(species_inout,
                                                    domain_indices,
                                                    active,
                                                    born,
                                                    n_ions,
                                                    reactions,
                                                    reaction_offsets,
                                                    n_domains,
                                                    rng_states,
                                                    dt);
}

__global__ void init_rng_kernel(curandStateXORWOW* states, unsigned long long seed, size_t n) {
    size_t idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (idx < n) {
        curand_init(seed, idx, 0, &states[idx]);
    }
}

void init_rng_states_xorwow(curandStateXORWOW* states,
                            size_t n_states,
                            unsigned long long seed,
                            cudaStream_t stream) {
    if (!states || n_states == 0) return;
    constexpr int THREADS = 256;
    int blocks = static_cast<int>((n_states + THREADS - 1) / THREADS);
    init_rng_kernel<<<blocks, THREADS, 0, stream>>>(states, seed, n_states);
}

}  // namespace gpu
}  // namespace icarion

#endif  // ICARION_USE_GPU

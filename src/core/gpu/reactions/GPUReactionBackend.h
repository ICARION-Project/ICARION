// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#pragma once

#ifdef ICARION_USE_GPU

#include "core/gpu/core/GPUContext.h"
#include "core/config/types/ReactionConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/config/types/DomainConfig.h"
#include "core/physics/reactions/IReactionHandler.h"
#include "reaction_kernels.cuh"
#include "core/types/IonEnsemble.h"
#include <vector>
#include <memory>
#include <curand_kernel.h>

namespace icarion::gpu {

/**
 * @brief Placeholder backend for GPU reaction handling.
 *
 * Current implementation is a stub that always falls back to CPU.
 * The class exists to provide a consistent architecture with other
 * GPU helpers (integration, collisions, space charge).
 */
class GPUReactionBackend {
public:
    struct Config {
        bool enable_logging = false;
        unsigned long long rng_seed = 0;
        size_t gpu_threshold = 2000;
    };

    GPUReactionBackend(std::shared_ptr<GPUContext> context, Config cfg);
    ~GPUReactionBackend();

    /**
     * @brief Attempt to process reactions on the GPU.
     * @return true if GPU path handled the batch, false to trigger CPU fallback.
     *
     * @note Currently unimplemented – always returns false.
     */
    bool process_batch(ICARION::core::IonEnsemble& ensemble,
                       const std::vector<int>& domain_indices,
                       double dt,
                       const ICARION::config::ReactionDatabase& reaction_db,
                       const ICARION::config::SpeciesDatabase& species_db,
                       const std::vector<ICARION::config::DomainConfig>& domains,
                       std::vector<ICARION::physics::PhysicsRng>& rng_pool);

    bool is_initialized() const { return context_ != nullptr; }

private:
    struct FlattenedReactions {
        std::vector<icarion::gpu::DeviceReaction> reactions;
        std::vector<int32_t> reaction_offsets;  // size = n_domains + 1
        int n_domains = 0;
    };

    FlattenedReactions flatten_reactions(
        const ICARION::config::ReactionDatabase& reaction_db,
        const ICARION::config::SpeciesDatabase& species_db,
        const ICARION::core::IonEnsemble& ensemble,
        const std::vector<ICARION::config::DomainConfig>& domains);

    void ensure_rng_states(size_t n);
    void ensure_ion_buffers(size_t n);
    void ensure_reaction_buffers(size_t reaction_count, size_t offset_count);

    std::shared_ptr<GPUContext> context_;
    Config config_;
    bool warned_once_ = false;
    std::unique_ptr<ICARION::physics::IReactionHandler> cpu_fallback_;

    // Device buffers
    curandStateXORWOW* d_rng_states_ = nullptr;
    uint32_t* d_species_idx_ = nullptr;
    int32_t* d_domain_idx_ = nullptr;
    uint8_t* d_active_ = nullptr;
    uint8_t* d_born_ = nullptr;
    icarion::gpu::DeviceReaction* d_reactions_ = nullptr;
    int32_t* d_offsets_ = nullptr;
    size_t ion_capacity_ = 0;
    size_t rng_capacity_ = 0;
    size_t reaction_capacity_ = 0;
    size_t offset_capacity_ = 0;
};

}  // namespace icarion::gpu

#endif  // ICARION_USE_GPU

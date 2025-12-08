// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#ifdef ICARION_USE_GPU

#include "core/gpu/core/GPUContext.h"
#include "core/config/types/ReactionConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonEnsemble.h"
#include <vector>

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
                       const std::vector<ICARION::config::DomainConfig>& domains);

    bool is_initialized() const { return context_ != nullptr; }

private:
    std::shared_ptr<GPUContext> context_;
    Config config_;
    bool warned_once_ = false;
};

}  // namespace icarion::gpu

#endif  // ICARION_USE_GPU

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once

#include "IReactionHandler.h"
#include <memory>

namespace icarion::gpu {
class GPUContext;
class GPUReactionBackend;
}  // namespace icarion::gpu

namespace ICARION::physics {

/**
 * @brief GPU-aware reaction handler wrapper.
 *
 * Delegates to the CPU handler by default, but exposes a batch interface
 * so a GPU backend can be plugged in without changing SimulationEngine.
 */
class GPUReactionHandler : public IReactionHandler {
public:
    GPUReactionHandler(std::unique_ptr<IReactionHandler> cpu_handler,
                       bool enable_logging,
                       unsigned long long rng_seed,
                       size_t gpu_threshold);
    ~GPUReactionHandler() override;

    bool handle_reaction(core::IonReactionData& view,
                         double dt,
                         PhysicsRng& rng,
                         const config::ReactionDatabase& reaction_db,
                         const config::SpeciesDatabase& species_db,
                         const config::EnvironmentConfig& env) override;

    std::string name() const override;
    ReactionStats get_stats() const override;
    void reset_stats() override;

    bool supports_batch() const override;
    bool handle_batch(core::IonEnsemble& ensemble,
                      const std::vector<int>& domain_indices,
                      double dt,
                      const config::ReactionDatabase& reaction_db,
                      const config::SpeciesDatabase& species_db,
                      const std::vector<config::DomainConfig>& domains,
                      std::vector<PhysicsRng>& rng_pool) override;

private:
    std::unique_ptr<IReactionHandler> cpu_handler_;
    bool logging_enabled_;
    size_t threshold_;

#ifdef ICARION_USE_GPU
    std::shared_ptr<icarion::gpu::GPUContext> gpu_context_;
    std::unique_ptr<icarion::gpu::GPUReactionBackend> gpu_backend_;
#endif
};

}  // namespace ICARION::physics

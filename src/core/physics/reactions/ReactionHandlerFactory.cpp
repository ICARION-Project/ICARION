// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include "ReactionHandlerFactory.h"
#include "StochasticReactionHandler.h"
#include "NoReactionHandler.h"
#include "core/log/Logger.h"
#include "GPUReactionHandler.h"

namespace ICARION {
namespace physics {

std::unique_ptr<IReactionHandler> ReactionHandlerFactory::create(
    const config::PhysicsConfig& config,
    bool enable_logging,
    bool enable_gpu,
    unsigned long long gpu_seed,
    size_t gpu_threshold
) {
#ifndef ICARION_USE_GPU
    (void)gpu_seed;
    (void)gpu_threshold;
#endif
    // SSOT: Read enable_reactions directly from config
    if (config.enable_reactions) {
        auto cpu_handler = std::make_unique<StochasticReactionHandler>(enable_logging);

#ifdef ICARION_USE_GPU
        if (enable_gpu) {
            return std::make_unique<GPUReactionHandler>(
                std::move(cpu_handler),
                enable_logging,
                gpu_seed,
                gpu_threshold
            );
        }
#else
        if (enable_gpu) {
            log::Logger::main()->warn(
                "ReactionHandlerFactory: GPU reactions requested but ICARION was built without GPU support");
        }
#endif
        return cpu_handler;
    } else {
        // Null Object Pattern: Return no-op handler (avoids null checks in caller)
        return std::make_unique<NoReactionHandler>();
    }
}

} // namespace physics
} // namespace ICARION

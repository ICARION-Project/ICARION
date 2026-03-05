// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include "GPUReactionHandler.h"
#include "core/log/Logger.h"

#ifdef ICARION_USE_GPU
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/reactions/GPUReactionBackend.h"
#endif

namespace ICARION::physics {

GPUReactionHandler::GPUReactionHandler(std::unique_ptr<IReactionHandler> cpu_handler,
                                       bool enable_logging,
                                       unsigned long long rng_seed,
                                       size_t gpu_threshold)
    : cpu_handler_(std::move(cpu_handler)),
      logging_enabled_(enable_logging),
      threshold_(gpu_threshold) {
#ifdef ICARION_USE_GPU
    try {
        gpu_context_ = icarion::gpu::GPUContext::create(0);
        if (gpu_context_) {
            icarion::gpu::GPUReactionBackend::Config cfg;
            cfg.enable_logging = enable_logging;
            cfg.rng_seed = rng_seed;
            cfg.gpu_threshold = gpu_threshold;
            gpu_backend_ = std::make_unique<icarion::gpu::GPUReactionBackend>(gpu_context_, cfg);
        }
    } catch (const std::exception& e) {
        log::Logger::main()->warn(
            "GPUReactionHandler: Failed to initialize GPU backend: {}", e.what());
        gpu_context_.reset();
        gpu_backend_.reset();
    }
#else
    (void)rng_seed;
#endif
}

GPUReactionHandler::~GPUReactionHandler() = default;

bool GPUReactionHandler::handle_reaction(core::IonReactionData& view,
                                         double dt,
                                         PhysicsRng& rng,
                                         const config::ReactionDatabase& reaction_db,
                                         const config::SpeciesDatabase& species_db,
                                         const config::EnvironmentConfig& env) {
    if (!cpu_handler_) {
        return false;
    }
    return cpu_handler_->handle_reaction(view, dt, rng, reaction_db, species_db, env);
}

std::string GPUReactionHandler::name() const {
    if (cpu_handler_) {
        return cpu_handler_->name();
    }
    return "GPUReactionHandler";
}

ReactionStats GPUReactionHandler::get_stats() const {
    if (cpu_handler_) {
        return cpu_handler_->get_stats();
    }
    return {};
}

void GPUReactionHandler::reset_stats() {
    if (cpu_handler_) {
        cpu_handler_->reset_stats();
    }
}

bool GPUReactionHandler::supports_batch() const {
#ifdef ICARION_USE_GPU
    return gpu_backend_ != nullptr;
#else
    return false;
#endif
}

bool GPUReactionHandler::handle_batch(core::IonEnsemble& ensemble,
                                      const std::vector<int>& domain_indices,
                                      double dt,
                                      const config::ReactionDatabase& reaction_db,
                                      const config::SpeciesDatabase& species_db,
                                      const std::vector<config::DomainConfig>& domains,
                                      std::vector<PhysicsRng>& rng_pool) {

#ifdef ICARION_USE_GPU
    if (!gpu_backend_) {
        return false;
    }
    if (ensemble.size() < threshold_) {
        return false;
    }
    return gpu_backend_->process_batch(
        ensemble, domain_indices, dt, reaction_db, species_db, domains, rng_pool);
#else
    (void)ensemble;
    (void)domain_indices;
    (void)dt;
    (void)reaction_db;
    (void)species_db;
    (void)domains;
    (void)rng_pool;
    return false;
#endif
}

}  // namespace ICARION::physics

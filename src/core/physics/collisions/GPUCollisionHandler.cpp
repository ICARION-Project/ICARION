// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include "GPUCollisionHandler.h"
#include "core/log/Logger.h"

namespace ICARION::physics {

GPUCollisionHandler::GPUCollisionHandler(
    std::unique_ptr<ICollisionHandler> cpu_handler,
    const std::string& collision_model,
    const GeometryMap* geometry_map,
    bool enable_logging,
    const config::SpeciesDatabase* species_db,
    unsigned long long rng_seed,
    size_t gpu_threshold)
    : cpu_handler_(std::move(cpu_handler)),
      model_name_(collision_model),
      logging_enabled_(enable_logging),
      threshold_(gpu_threshold),
      species_db_(species_db) {
    if (geometry_map) {
        geometry_copy_ = std::make_unique<GeometryMap>(*geometry_map);
    }

#ifdef ICARION_USE_GPU
    try {
        gpu_context_ = icarion::gpu::GPUContext::create(0);
        if (gpu_context_) {
            gpu_helper_ = icarion::gpu::GPUCollisionHelper::create(
                *gpu_context_, threshold_, model_name_, rng_seed);
            if (gpu_helper_ && geometry_copy_ && model_name_ == "EHSS") {
                gpu_helper_->set_geometry(*geometry_copy_);
            }
        }
    } catch (const std::exception& e) {
        log::Logger::main()->warn("GPUCollisionHandler: Failed to initialize GPU helper: {}", e.what());
        gpu_context_.reset();
        gpu_helper_.reset();
    }
#else
    (void)rng_seed;
#endif
}

bool GPUCollisionHandler::handle_collision(core::IonCollisionData& view,
                                           double dt,
                                           PhysicsRng& rng,
                                           const config::EnvironmentConfig& env) {
    return cpu_handler_ && cpu_handler_->handle_collision(view, dt, rng, env);
}

std::string GPUCollisionHandler::name() const {
    if (cpu_handler_) {
        return cpu_handler_->name();
    }
    return model_name_;
}

CollisionStats GPUCollisionHandler::get_stats() const {
    if (cpu_handler_) {
        return cpu_handler_->get_stats();
    }
    return {};
}

void GPUCollisionHandler::reset_stats() {
    if (cpu_handler_) {
        cpu_handler_->reset_stats();
    }
#ifdef ICARION_USE_GPU
    if (gpu_helper_) {
        gpu_helper_->reset_stats();
    }
#endif
}

bool GPUCollisionHandler::supports_batch() const {
#ifdef ICARION_USE_GPU
    return static_cast<bool>(gpu_helper_);
#else
    return false;
#endif
}

bool GPUCollisionHandler::handle_batch(core::IonEnsemble& ensemble,
                                       const std::vector<size_t>& ion_indices,
                                       double dt,
                                       const config::EnvironmentConfig& env,
                                       std::vector<physics::PhysicsRng>& rng_pool) {
    (void)rng_pool;
    if (ion_indices.empty()) {
        return true;
    }

#ifdef ICARION_USE_GPU
    if (!supports_batch() || ion_indices.size() < threshold_) {
        return false;
    }

    std::vector<IonState> ions;
    ions.reserve(ion_indices.size());
    for (size_t idx : ion_indices) {
        ions.push_back(ensemble.ion_state(idx));
    }

    if (!gpu_helper_->process_collisions_batch(ions, dt, env)) {
        if (!logging_enabled_) {
            log::Logger::main()->warn(
                "GPUCollisionHandler: GPU batch failed for {} – fallback to CPU",
                model_name_);
        }
        return false;
    }

    for (size_t i = 0; i < ion_indices.size(); ++i) {
        ensemble.apply_ion_state(ion_indices[i], ions[i]);
    }
    return true;
#else
    (void)ensemble;
    (void)ion_indices;
    (void)dt;
    (void)env;
    return false;
#endif
}

}  // namespace ICARION::physics

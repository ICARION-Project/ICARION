// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once

#include "ICollisionHandler.h"
#include "EHSSCollisionHandler.h"
#include <memory>

#ifdef ICARION_USE_GPU
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#endif

namespace ICARION::physics {

class GPUCollisionHandler : public ICollisionHandler {
public:
    GPUCollisionHandler(std::unique_ptr<ICollisionHandler> cpu_handler,
                        const std::string& collision_model,
                        const GeometryMap* geometry_map,
                        bool enable_logging,
                        const config::SpeciesDatabase* species_db,
                        unsigned long long rng_seed,
                        size_t gpu_threshold = 5000);

    bool handle_collision(core::IonCollisionData& view,
                          double dt,
                          PhysicsRng& rng,
                          const config::EnvironmentConfig& env) override;

    std::string name() const override;
    CollisionStats get_stats() const override;
    void reset_stats() override;

    bool supports_batch() const override;
    bool handle_batch(core::IonEnsemble& ensemble,
                      const std::vector<size_t>& ion_indices,
                      double dt,
                      const config::EnvironmentConfig& env,
                      std::vector<physics::PhysicsRng>& rng_pool) override;

private:
    std::unique_ptr<ICollisionHandler> cpu_handler_;
    std::string model_name_;
    bool logging_enabled_;
    size_t threshold_;
    std::unique_ptr<GeometryMap> geometry_copy_;
    const config::SpeciesDatabase* species_db_;

#ifdef ICARION_USE_GPU
    std::shared_ptr<icarion::gpu::GPUContext> gpu_context_;
    std::shared_ptr<icarion::gpu::GPUCollisionHelper> gpu_helper_;
#endif
};

}  // namespace ICARION::physics

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#pragma once

#include <vector>
#include <memory>
#include <optional>
#include <cmath>

#include "core/config/types/FullConfig.h"
#include "core/config/types/PhysicsEnums.h"
#include "core/config/conversion/EnumMapper.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/physics/forces/DampingForce.h"
#include "core/physics/collisions/CollisionHandlerFactory.h"
#include "core/types/CollisionTypes.h"
#include "core/integrator/SimulationEngine.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/types/IonEnsemble.h"

namespace ICARION::tests {

using physics::PhysicsRng;

struct SimpleSimTrace {
    std::vector<double> times;
    std::vector<double> z_positions;  // First ion z-trace (optional)
};

struct SimpleSimResult {
    std::vector<core::IonState> ions;
    SimpleSimTrace trace;
};

inline config::DomainConfig make_single_domain(config::Instrument inst) {
    config::DomainConfig dom;
    dom.instrument = inst;
    dom.name = "test-domain";
    return dom;
}

inline std::shared_ptr<physics::ForceRegistry> build_force_registry(
    const config::DomainConfig& dom,
    config::CollisionModel collision_model,
    const config::SpeciesDatabase* species_db,
    bool enable_ou
) {
    (void)enable_ou;
    auto registry = std::make_shared<physics::ForceRegistry>(dom);
    registry->add_force(std::make_unique<physics::ElectricFieldForce>(dom));

    // Deterministic damping models: Friction, Langevin, HSD
    //
    // These models use DampingForce (continuous friction in RK4)
    // Optional: Add OU handler for thermal kicks (stochastic)
    //
    // Implementation:
    // - DampingForce: F = -γ·m·v (continuous, deterministic)
    // - OU handler (if enabled): Adds thermal diffusion Δv ~ √(2γkBT dt)
    //
    // The OU handler uses apply_damping=false to avoid double-damping!
    physics::DampingModel damping = physics::DampingModel::None;
    switch (collision_model) {
        case config::CollisionModel::Friction:  damping = physics::DampingModel::Friction; break;
        case config::CollisionModel::Langevin:  damping = physics::DampingModel::Langevin; break;
        case config::CollisionModel::HSD:       damping = physics::DampingModel::HardSphere; break;
        default: break;
    }
    if (damping != physics::DampingModel::None) {
        // Always add DampingForce for deterministic models
        // (Provides continuous friction in RK4 integration)
        registry->add_force(std::make_unique<physics::DampingForce>(
            dom.environment, damping, species_db));
    }

    return registry;
}

inline SimpleSimResult run_simple_simulation(
    config::FullConfig cfg,
    std::vector<core::IonState> ions,
    bool capture_trace = false,
    const physics::GeometryMap* geometry_map = nullptr,
    double gamma_for_ou = 0.0
) {
    // Finalize config (derived quantities)
    cfg.finalize_all();

    // Build force registries (one per domain)
    std::vector<std::shared_ptr<physics::ForceRegistry>> registries;
    registries.reserve(cfg.domains.size());
    for (const auto& dom : cfg.domains) {
        registries.push_back(build_force_registry(dom, cfg.physics.collision_model, &cfg.species_db, cfg.physics.enable_ou_thermalization));
    }

    // Collision handler (stochastic models only)
    std::shared_ptr<physics::ICollisionHandler> collision_handler;
    {
        auto handler = physics::CollisionHandlerFactory::create(
            cfg.physics,
            geometry_map,     // geometry_map (needed for EHSS)
            gamma_for_ou,     // gamma_for_ou (for OU with deterministic models)
            false,            // logging
            &cfg.species_db
        );
        if (handler) {
            collision_handler = std::move(handler);
        }
    }

    // Integrator
    auto integrator = std::make_shared<integrator::RK4Strategy>();

    SimpleSimResult result;

    if (!capture_trace) {
        // Use full SimulationEngine path for stability and SSOT alignment
        integrator::SimulationEngine engine(cfg, registries, integrator, collision_handler, nullptr);
        auto ensemble = core::IonEnsemble::from_legacy(ions);
        auto final_ensemble = engine.run(ensemble);
        result.ions = final_ensemble.to_legacy();
        return result;
    }

    // Manual loop for trace capture (single domain assumption, SoA)
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    const auto& dom = cfg.domains.front();
    const double dt = cfg.simulation.dt_s;
    const int total_steps = cfg.simulation.total_steps;
    PhysicsRng rng(1337);

    if (total_steps > 0 && total_steps < 10'000'000) {
        result.trace.times.reserve(total_steps + 1);
        result.trace.z_positions.reserve(total_steps + 1);
    }
    result.trace.times.push_back(0.0);
    result.trace.z_positions.push_back(ensemble.get_pos(0).z);

    for (int step = 0; step < total_steps; ++step) {
        // Collisions (stochastic)
        if (collision_handler) {
            for (size_t i = 0; i < ensemble.size(); ++i) {
                auto view = ensemble.collision_data(i);
                collision_handler->handle_collision(view, dt, rng, dom.environment);
            }
        }

        // Deterministic forces via RK4 (SoA)
        for (size_t i = 0; i < ensemble.size(); ++i) {
            integrator->step(ensemble, i, step * dt, dt, *registries.front());
        }

        result.trace.times.push_back((step + 1) * dt);
        result.trace.z_positions.push_back(ensemble.get_pos(0).z);
    }

    result.ions = ensemble.to_legacy();
    return result;
}

} // namespace ICARION::tests

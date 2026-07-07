// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#define private public
#include "core/integrator/SimulationEngine.h"
#undef private
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/collisions/OUCollisionHandler.h"
#include "core/config/types/FullConfig.h"
#include "core/types/IonEnsemble.h"

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;
using Catch::Matchers::WithinAbs;

namespace {
class ZeroForce : public IForce {
public:
    Vec3 compute(const core::IonEnsemble&, size_t, double, const ForceContext&) const override {
        return Vec3{0.0, 0.0, 0.0};
    }
    Vec3 compute_soa(const ForceState&, double, const ForceContext&) const override {
        return Vec3{0.0, 0.0, 0.0};
    }
    bool applies_to(const IonState&) const override { return true; }
    std::string name() const override { return "ZeroForce"; }
};
}

TEST_CASE("SimulationEngine RNG determinism after compaction with per-ion dt", "[engine][rng][compaction]") {
    FullConfig cfg;
    cfg.simulation.total_time_s = 1e-5;
    cfg.simulation.dt_s = 1e-5;
    cfg.simulation.integrator = "RK4";
    cfg.simulation.rng_seed = 12345;
    cfg.physics.collision_model = config::CollisionModel::Friction;
    cfg.physics.enable_ou_thermalization = true;
    cfg.physics.enable_reactions = false;
    cfg.output.folder = "/tmp";
    DomainConfig domain;
    domain.name = "rng_compaction_domain";
    domain.instrument = Instrument::IMS;
    domain.geometry.length_m = 0.1;
    domain.geometry.radius_m = 0.01;
    domain.environment.temperature_K = 300.0;
    domain.environment.pressure_Pa = 101325.0;
    domain.finalize();
    cfg.domains.push_back(domain);

    auto registry = std::make_shared<ForceRegistry>(cfg.domains[0]);
    registry->add_force(std::make_unique<ZeroForce>());
    std::vector<std::shared_ptr<ForceRegistry>> registries = {registry};
    auto integrator = std::make_shared<RK4Strategy>();

    // Two runs: baseline and with compaction after first step
    auto run_engine = [&](bool compact_after_first) {
        auto collision_handler = std::make_shared<OUCollisionHandler>(1.0e6, true);
        SimulationEngine eng(cfg, registries, integrator, collision_handler, nullptr);
        core::IonEnsemble ens;
        ens.resize(3);
        ens.prepopulate_species_pool({"test+"});
        for (size_t i = 0; i < 3; ++i) {
            ens.mass_data()[i] = 1.0;
            ens.charge_data()[i] = 1.0;
            ens.mobility_data()[i] = 1.0;
            ens.CCS_data()[i] = 1e-18;
            ens.born_data()[i] = 1;
            ens.active_data()[i] = 1;
            ens.pos_x_data()[i] = static_cast<double>(i);
            ens.pos_y_data()[i] = 0.0;
            ens.pos_z_data()[i] = 0.05;
            ens.vel_x_data()[i] = 100.0 + static_cast<double>(i);
            ens.vel_y_data()[i] = 50.0;
            ens.vel_z_data()[i] = 0.0;
            ens.time_data()[i] = 0.0;
            ens.reaction_data(i).set_species_id("test+");
        }
        eng.initialize(ens);
        eng.current_time_ = 0.0;
        eng.current_step_ = 0;
        eng.dt_per_ion_.assign(ens.size(), cfg.simulation.dt_s);
        eng.deep_collision_diagnostics_.reset(ens.size());
        eng.collision_runtime_stats_.reset();

        auto do_step = [&]() {
            const double t_new = eng.process_timestep(ens);
            eng.current_time_ = t_new;
            eng.current_step_++;
        };

        // Step 1 with identical initial state in both paths.
        do_step();

        // Deactivate one ion between steps; only one path compacts.
        ens.active_data()[1] = 0;
        if (compact_after_first) {
            ens.compact_inactive();
        }

        eng.dt_per_ion_.assign(ens.size(), cfg.simulation.dt_s);
        do_step();

        std::vector<double> out;
        for (size_t i = 0; i < ens.size(); ++i) {
            if (!ens.is_active(i)) {
                continue;
            }
            out.push_back(ens.pos_x_data()[i]);
            out.push_back(ens.vel_x_data()[i]);
        }
        return out;
    };

    auto pos_baseline = run_engine(false);
    auto pos_compact = run_engine(true);

    REQUIRE(pos_baseline.size() == pos_compact.size());
    for (size_t i = 0; i < pos_baseline.size(); ++i) {
        REQUIRE_THAT(pos_baseline[i], WithinAbs(pos_compact[i], 1e-12));
    }
}

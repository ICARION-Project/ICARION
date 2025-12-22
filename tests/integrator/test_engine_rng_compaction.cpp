// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include "core/integrator/SimulationEngine.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/config/types/FullConfig.h"
#include "core/types/IonEnsemble.h"

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;

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
    cfg.simulation.total_time_s = 1e-4;
    cfg.simulation.dt_s = 1e-5;
    cfg.simulation.integrator = "RK4";
    cfg.simulation.rng_seed = 12345;
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    cfg.physics.enable_reactions = false;
    cfg.output.folder = "/tmp";
    cfg.domains.push_back(DomainConfig{});

    auto registry = std::make_shared<ForceRegistry>(cfg.domains[0]);
    registry->add_force(std::make_unique<ZeroForce>());
    std::vector<std::shared_ptr<ForceRegistry>> registries = {registry};
    auto integrator = std::make_shared<RK4Strategy>();

    // Two runs: baseline and with compaction after first step
    auto run_engine = [&](bool compact_after_first) {
        SimulationEngine eng(cfg, registries, integrator, nullptr, nullptr);
        core::IonEnsemble ens;
        ens.resize(3);
        for (size_t i = 0; i < 3; ++i) {
            ens.mass_data()[i] = 1.0;
            ens.charge_data()[i] = 1.0;
            ens.born_data()[i] = 1;
            ens.active_data()[i] = 1;
            ens.pos_x_data()[i] = static_cast<double>(i);
            ens.time_data()[i] = 0.0;
        }
        if (compact_after_first) {
            eng.run(ens); // first run fills RNGs/dt state
            // deactivate one ion and compact
            ens.active_data()[1] = 0;
            ens.compact_inactive();
            eng.run(ens);
        } else {
            eng.run(ens);
            eng.run(ens);
        }
        return std::vector<double>(ens.pos_x_data(), ens.pos_x_data() + ens.size());
    };

    auto pos_baseline = run_engine(false);
    auto pos_compact = run_engine(true);

    REQUIRE(pos_baseline.size() == pos_compact.size());
    for (size_t i = 0; i < pos_baseline.size(); ++i) {
        REQUIRE(pos_baseline[i] == Approx(pos_compact[i]));
    }
}

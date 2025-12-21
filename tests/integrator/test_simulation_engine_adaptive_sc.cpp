// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/integrator/SimulationEngine.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/reactions/NoReactionHandler.h"
#include "core/physics/spacecharge/ISpaceChargeModel.h"
#include "core/types/IonState.h"
#include <cstdlib>
#include <memory>
#include <optional>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using Catch::Matchers::WithinRel;

namespace {

struct EnvGuard {
    explicit EnvGuard(const char* key, const char* value) : key_(key) {
        const char* existing = std::getenv(key);
        if (existing) prev_ = existing;
        setenv(key, value, 1);
    }
    ~EnvGuard() {
        if (prev_) {
            setenv(key_, prev_->c_str(), 1);
        } else {
            unsetenv(key_);
        }
    }
    const char* key_;
    std::optional<std::string> prev_;
};

class DummySpaceChargeModel : public ISpaceChargeModel {
public:
    void update_fields(const core::IonEnsemble&, double) override { ++updates_; }
    core::Vec3 sample_electric_field(std::size_t) const override { return {0.0, 0.0, 0.0}; }
    std::string name() const override { return "dummy_sc"; }
    int updates() const { return updates_; }
private:
    int updates_ = 0;
};

config::FullConfig make_config() {
    config::FullConfig cfg;
    cfg.simulation.total_time_s = 2e-9;  // two nanoseconds (one RK45 step)
    cfg.simulation.dt_s = 1e-9;
    cfg.simulation.compute_derived();

    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    cfg.physics.enable_reactions = false;

    cfg.output.folder = "/tmp";
    cfg.output.trajectory_file = "adaptive_sc_test.h5";

    config::DomainConfig dom;
    dom.name = "test";
    dom.instrument = config::Instrument::IMS;
    dom.geometry.length_m = 0.01;
    dom.geometry.radius_m = 0.005;
    dom.environment.temperature_K = 300.0;
    dom.environment.pressure_Pa = 101325.0;
    dom.finalize();
    cfg.domains.push_back(dom);

    return cfg;
}

core::IonState make_ion() {
    core::IonState ion;
    ion.pos = {0.0, 0.0, 0.001};
    ion.vel = {0.0, 0.0, 100.0};
    ion.mass_kg = 100.0 * 1.66053906660e-27;
    ion.ion_charge_C = 1.602176634e-19;
    ion.species_id = "test_species";
    ion.current_domain_index = 0;
    ion.active = true;
    ion.born = true;
    ion.birth_time_s = 0.0;
    return ion;
}

} // namespace

TEST_CASE("SimulationEngine adaptive SC: enabled path executes", "[simulation][adaptive_sc]") {
    unsetenv("ICARION_ADAPTIVE_SC");
    auto cfg = make_config();

    auto sc_model = std::make_shared<DummySpaceChargeModel>();
    auto registry = std::make_shared<ForceRegistry>(cfg.domains[0]);
    registry->set_space_charge_model(sc_model);

    auto integrator = std::make_shared<RK45Strategy>();
    auto reaction_handler = std::make_shared<NoReactionHandler>();

    SimulationEngine engine(cfg, {registry}, integrator, nullptr, reaction_handler);

    std::vector<core::IonState> ions{make_ion()};
    auto ensemble = core::IonEnsemble::from_legacy(ions);
    auto result_ens = engine.run(ensemble);
    auto out = result_ens.to_legacy();

    REQUIRE(out.size() == 1);
    REQUIRE(out[0].pos.z > ions[0].pos.z);
    REQUIRE(sc_model->updates() > 0); // staged SC rebuilds occurred
}

TEST_CASE("SimulationEngine adaptive SC: env guard disables path", "[simulation][adaptive_sc]") {
    EnvGuard guard("ICARION_ADAPTIVE_SC", "0");
    auto cfg = make_config();

    auto sc_model = std::make_shared<DummySpaceChargeModel>();
    auto registry = std::make_shared<ForceRegistry>(cfg.domains[0]);
    registry->set_space_charge_model(sc_model);

    auto integrator = std::make_shared<RK45Strategy>();
    auto reaction_handler = std::make_shared<NoReactionHandler>();

    SimulationEngine engine(cfg, {registry}, integrator, nullptr, reaction_handler);

    std::vector<core::IonState> ions{make_ion()};
    auto ensemble = core::IonEnsemble::from_legacy(ions);

    REQUIRE_THROWS_AS(engine.run(ensemble), std::runtime_error);
}

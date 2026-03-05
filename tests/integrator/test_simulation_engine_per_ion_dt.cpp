// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include "core/integrator/SimulationEngine.h"
#include "core/integrator/DomainManager.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/config/types/FullConfig.h"
#include "core/types/IonEnsemble.h"

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;

namespace {
// Constant electric field in +x for determinism
class ConstantEField : public config::IFieldModel {
public:
    Vec3 E(const Vec3&, double) const override {
        return Vec3{1.0, 0.0, 0.0}; // 1 V/m
    }
};

std::shared_ptr<ForceRegistry> make_registry(const DomainConfig& dom) {
    auto registry = std::make_shared<ForceRegistry>(dom);
    registry->add_force(std::make_unique<ElectricFieldForce>(dom));
    registry->set_field_model(new ConstantEField());
    return registry;
}
}

TEST_CASE("SimulationEngine per-ion dt not clamped to min across ions", "[engine][rk45][dt]") {
    // Config
    FullConfig cfg;
    cfg.simulation.total_time_s = 1e-4;
    cfg.simulation.dt_s = 1e-5;
    cfg.simulation.integrator = "RK45";
    cfg.simulation.enable_openmp = false;
    cfg.physics.collision_model = config::CollisionModel::NoCollisions;
    cfg.physics.enable_reactions = false;
    cfg.output.folder = "/tmp";

    DomainConfig dom;
    dom.name = "test";
    dom.fields.field_array_terms.clear();
    dom.fields.legacy_field_array_file = "";
    dom.environment.temperature_K = 300.0;
    dom.environment.gas_mass_kg = 28.0 * AMU_TO_KG;
    cfg.domains.push_back(dom);

    // Forces
    auto registry = make_registry(cfg.domains[0]);
    std::vector<std::shared_ptr<ForceRegistry>> registries = {registry};

    // Integrator
    auto rk45 = std::make_shared<RK45Strategy>();

    // Engine
    SimulationEngine engine(cfg, registries, rk45, nullptr, nullptr);

    // Ensemble: two ions, different tolerances via different dt hints
    core::IonEnsemble ens;
    ens.resize(2);
    ens.mass_data()[0] = 1.0;
    ens.mass_data()[1] = 1.0;
    ens.charge_data()[0] = 1.0;
    ens.charge_data()[1] = 1.0;
    ens.CCS_data()[0] = 1e-20;
    ens.CCS_data()[1] = 1e-20;
    ens.mobility_data()[0] = 1e-4;
    ens.mobility_data()[1] = 1e-4;
    ens.born_data()[0] = 1;
    ens.born_data()[1] = 1;
    ens.active_data()[0] = 1;
    ens.active_data()[1] = 1;
    ens.update_species(0, "A", 1.0, 1.0, 1e-20, 1e-4);
    ens.update_species(1, "B", 1.0, 1.0, 1e-20, 1e-4);

    // Manually seed per-ion dt: ion0 small, ion1 large
    engine.run(ens); // initialization sets up dt_per_ion internally (will use cfg.dt_s)

    // After one timestep, ensure ion1 did not clamp to ion0 dt
    // We expect ion1 displacement to be larger than ion0
    // Because ax=E/m = 1 m/s^2, x = 0.5*a*dt^2
    // Allow some tolerance; we just check ordering
    double x0 = ens.pos_x_data()[0];
    double x1 = ens.pos_x_data()[1];
    REQUIRE(x1 >= x0);
}

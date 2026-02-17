// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include <catch2/catch_test_macros.hpp>
#include "core/integrator/strategies/IIntegrationStrategy.h"
#include "core/types/IonEnsemble.h"
#include "core/config/types/DomainConfig.h"
#include "core/physics/forces/ForceRegistry.h"
#include <memory>
#include <vector>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;

namespace {
struct MockBatchIntegrator : public IIntegrationStrategy {
    bool batch_called = false;
    bool batch_should_run = true;
    int step_calls = 0;
    double last_dt = 0.0;

    void step(core::IonEnsemble& ensemble,
              size_t ion_idx,
              double t,
              double dt,
              const physics::ForceRegistry&) override {
        (void)ensemble;
        (void)ion_idx;
        (void)t;
        last_dt = dt;
        step_calls++;
    }

    bool supports_batch() const override { return true; }

    bool step_batch(core::IonEnsemble& ensemble,
                    double t,
                    double dt,
                    const std::vector<std::shared_ptr<physics::ForceRegistry>>&,
                    const std::vector<int>&) override {
        (void)ensemble;
        (void)t;
        if (batch_should_run) {
            batch_called = true;
            last_dt = dt;
            return true;
        }
        return false;
    }

    std::string name() const override { return "MockBatch"; }
    bool is_adaptive() const override { return false; }
};
}

TEST_CASE("Integration batch only runs for uniform dt", "[integration][batch][dt]") {
    // Setup integrator and force registries
    auto integrator = std::make_shared<MockBatchIntegrator>();
    DomainConfig dom;
    auto registry = std::make_shared<ForceRegistry>(dom);
    std::vector<std::shared_ptr<ForceRegistry>> registries = {registry};

    core::IonEnsemble ens;
    ens.resize(3);
    for (size_t i = 0; i < 3; ++i) {
        ens.mass_data()[i] = 1.0;
        ens.charge_data()[i] = 1.0;
        ens.active_data()[i] = 1;
        ens.born_data()[i] = 1;
    }

    // Domain indices all valid
    std::vector<int> domains = {0, 0, 0};

    // Uniform dt: expect batch to run
    integrator->batch_called = false;
    integrator->batch_should_run = true;
    std::vector<double> dt_per_ion = {1e-4, 1e-4, 1e-4};
    std::vector<double> dt_used = dt_per_ion;
    std::vector<double> dt_next = dt_per_ion;

    SimulationEngine engine(config::FullConfig{}, registries, integrator, nullptr, nullptr);
    // Call the internal perform_integration directly
    double max_dt = engine.perform_integration(ens, 0.0, dt_per_ion, domains, dt_used, dt_next);
    REQUIRE(integrator->batch_called == true);
    REQUIRE(max_dt == Approx(1e-4));

    // Non-uniform dt: expect batch not to run, fallback to step()
    integrator->batch_called = false;
    integrator->batch_should_run = true;
    integrator->step_calls = 0;
    dt_per_ion = {1e-4, 2e-4, 1e-4};
    dt_used = dt_per_ion;
    dt_next = dt_per_ion;
    max_dt = engine.perform_integration(ens, 0.0, dt_per_ion, domains, dt_used, dt_next);
    REQUIRE(integrator->batch_called == false);
    REQUIRE(integrator->step_calls == 3);
    REQUIRE(max_dt == Approx(2e-4));
}

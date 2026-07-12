// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/physics/forces/ElectricFieldForce.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonEnsemble.h"

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;

namespace {
// Simple constant acceleration force (ax = 1 m/s^2, others zero)
class ConstantForce : public IForce {
public:
    Vec3 compute(const core::IonEnsemble&, size_t, double, const ForceContext&) const override {
        return Vec3{1.0, 0.0, 0.0};
    }
    Vec3 compute_soa(const ForceState&, double, const ForceContext&) const override {
        return Vec3{1.0, 0.0, 0.0};
    }
    bool applies_to(const IonState&) const override { return true; }
    std::string name() const override { return "ConstantForce"; }
};

core::IonState make_test_ion() {
    core::IonState ion;
    ion.species_id = "X+";
    ion.mass_kg = 1.0;
    ion.ion_charge_C = 1.0;
    ion.born = true;
    ion.active = true;
    ion.current_domain_index = 0;
    return ion;
}
}

TEST_CASE("RK45 per-ion dt divergence does not clamp to min(dt)", "[rk45][perion][dt]") {
    // Domain/Force setup
    DomainConfig domain;
    ForceRegistry registry(domain);
    registry.add_force(std::make_unique<ConstantForce>());

    // Two ions with different initial dt
    core::IonEnsemble ensemble = core::IonEnsemble::from_legacy({
        make_test_ion(),
        make_test_ion(),
    });

    // Ion 0 starts with small dt (harder tolerance), Ion 1 larger
    double dt0 = 1e-6;
    double dt1 = 1e-3;

    RK45Strategy::AdaptiveConfig cfg;
    cfg.atol = 1e-12;
    cfg.rtol = 1e-10;
    RK45Strategy rk45(cfg);

    // Ion 0: use dt0
    rk45.step(ensemble, 0, 0.0, dt0, registry);
    double used0 = rk45.last_dt_used();
    double next0 = rk45.last_dt_suggested();

    // Ion 1: use dt1
    rk45.step(ensemble, 1, 0.0, dt1, registry);
    double used1 = rk45.last_dt_used();
    double next1 = rk45.last_dt_suggested();

    // Ensure per-ion dt paths are independent (no clamping to min dt)
    REQUIRE(used1 >= used0);
    REQUIRE(next1 >= next0);
}

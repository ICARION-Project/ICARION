// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>

#include "core/config/types/FullConfig.h"
#include "core/config/types/InstrumentTypes.h"
#include "core/physics/spacecharge/SpaceChargeModelFactory.h"
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"

using ICARION::config::DomainConfig;
using ICARION::config::FullConfig;
using ICARION::config::Instrument;
using ICARION::config::SpaceChargeModel;
using ICARION::physics::SpaceChargeModelFactory;
using namespace ICARION;

namespace {

FullConfig make_config(bool enable_space_charge) {
    FullConfig cfg;
    cfg.physics.enable_space_charge = enable_space_charge;
    return cfg;
}

FullConfig make_config(bool enable_space_charge, SpaceChargeModel model_type) {
    FullConfig cfg = make_config(enable_space_charge);
    cfg.physics.space_charge_model_type = model_type;
    return cfg;
}

DomainConfig make_domain() {
    DomainConfig dom;
    dom.instrument = Instrument::IMS;
    dom.name = "factory_test_domain";
    dom.geometry.length_m = 0.1;
    dom.geometry.radius_m = 0.01;
    dom.geometry.origin_m = Vec3{0.0, 0.0, 0.0};
    dom.rotation_global_to_local = Mat3::identity();
    dom.rotation_local_to_global = Mat3::identity();
    dom.finalize();
    return dom;
}

} // namespace

TEST_CASE("SpaceChargeModelFactory returns nullptr when disabled or empty",
          "[spacecharge][factory]") {
    auto domain = make_domain();

    SECTION("space charge disabled") {
        auto cfg = make_config(false);
        auto model = SpaceChargeModelFactory::create(cfg, domain, 500);
        REQUIRE(model == nullptr);
    }

    SECTION("zero ions") {
        auto cfg = make_config(true);
        auto model = SpaceChargeModelFactory::create(cfg, domain, 0);
        REQUIRE(model == nullptr);
    }
}

TEST_CASE("SpaceChargeModelFactory auto-selects direct or grid models",
          "[spacecharge][factory]") {
    auto cfg = make_config(true);
    auto domain = make_domain();

    SECTION("small ion count creates direct model") {
        constexpr std::size_t ion_count = SpaceChargeModelFactory::DIRECT_THRESHOLD - 1;
        auto model = SpaceChargeModelFactory::create(cfg, domain, ion_count);

        REQUIRE(model != nullptr);
        REQUIRE(model->name() == "SpaceChargeDirectModel");
    }

    SECTION("large ion count creates grid model") {
        constexpr std::size_t ion_count = SpaceChargeModelFactory::GRID_THRESHOLD + 1;
        auto model = SpaceChargeModelFactory::create(cfg, domain, ion_count);

        REQUIRE(model != nullptr);
        REQUIRE(model->name() == "SpaceChargeGridModel");
    }
}

TEST_CASE("SpaceChargeModelFactory honors explicit CPU model selection",
          "[spacecharge][factory]") {
    auto domain = make_domain();

    SECTION("direct model ignores large ion count") {
        auto cfg = make_config(true, SpaceChargeModel::Direct);
        auto model = SpaceChargeModelFactory::create(
            cfg, domain, SpaceChargeModelFactory::GRID_THRESHOLD + 1);

        REQUIRE(model != nullptr);
        REQUIRE(model->name() == "SpaceChargeDirectModel");
    }

    SECTION("grid model ignores small ion count") {
        auto cfg = make_config(true, SpaceChargeModel::Grid);
        auto model = SpaceChargeModelFactory::create(cfg, domain, 1);

        REQUIRE(model != nullptr);
        REQUIRE(model->name() == "SpaceChargeGridModel");
    }
}

TEST_CASE("SpaceChargeModelFactory created model is callable",
          "[spacecharge][factory]") {
    auto cfg = make_config(true);
    auto domain = make_domain();
    auto model = SpaceChargeModelFactory::create(cfg, domain, 2);
    REQUIRE(model != nullptr);

    ICARION::core::IonState ion;
    ion.pos = Vec3{0.0, 0.0, 0.0};
    ion.vel = Vec3{0.0, 0.0, 0.0};
    ion.ion_charge_C = 1.6e-19;
    ion.mass_kg = 1e-26;
    ion.active = true;
    ion.born = true;
    ICARION::core::IonEnsemble ensemble = ICARION::core::IonEnsemble::from_legacy({ion});

    REQUIRE_NOTHROW(model->update_fields(ensemble, 0.0));
    auto E = model->sample_electric_field(0);
    REQUIRE(E.x == 0.0);
    REQUIRE(E.y == 0.0);
    REQUIRE(E.z == 0.0);
}

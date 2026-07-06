// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/config/types/TIMSAxialFieldModel.h"
#include "core/physics/forces/ElectricFieldForce.h"

using Catch::Approx;

namespace {

ICARION::config::DomainConfig make_tims_domain() {
    ICARION::config::DomainConfig dom;
    dom.instrument = ICARION::config::Instrument::TIMS;
    dom.geometry.origin_m = {0.0, 0.0, 0.0};
    dom.geometry.length_m = 0.02;
    dom.geometry.radius_m = 0.005;
    dom.fields.tims.enabled = true;
    return dom;
}

} // namespace

TEST_CASE("TIMSAxialFieldModel evaluates linear ramp", "[field][tims]") {
    auto dom = make_tims_domain();
    dom.fields.tims.axial_field_initial_uniform_V_m = 1000.0;
    dom.fields.tims.axial_field_final_uniform_V_m = 0.0;
    dom.fields.tims.ramp_start_s = 1.0e-3;
    dom.fields.tims.ramp_end_s = 2.0e-3;

    ICARION::config::TIMSAxialFieldModel model(dom);
    const Vec3 p{0.0, 0.0, 0.01};

    REQUIRE(model.E(p, 0.5e-3).z == Approx(1000.0));
    REQUIRE(model.E(p, 1.5e-3).z == Approx(500.0));
    REQUIRE(model.E(p, 2.5e-3).z == Approx(0.0));
}

TEST_CASE("TIMSAxialFieldModel interpolates axial profiles", "[field][tims]") {
    auto dom = make_tims_domain();
    dom.fields.tims.z_positions_m = {0.0, 0.02};
    dom.fields.tims.axial_field_initial_profile_V_m = {0.0, -4000.0};
    dom.fields.tims.axial_field_final_profile_V_m = {0.0, -2000.0};
    dom.fields.tims.ramp_start_s = 0.0;
    dom.fields.tims.ramp_end_s = 1.0;

    ICARION::config::TIMSAxialFieldModel model(dom);
    const Vec3 p{0.0, 0.0, 0.01};

    REQUIRE(model.E(p, 0.0).z == Approx(-2000.0));
    REQUIRE(model.E(p, 0.5).z == Approx(-1500.0));
    REQUIRE(model.E(p, 1.0).z == Approx(-1000.0));
}

TEST_CASE("ElectricFieldForce uses TIMS model for TIMS domains", "[forces][electric][tims]") {
    auto dom = make_tims_domain();
    dom.fields.tims.axial_field_initial_uniform_V_m = 800.0;
    dom.fields.tims.axial_field_final_uniform_V_m = 400.0;
    dom.fields.tims.ramp_start_s = 0.0;
    dom.fields.tims.ramp_end_s = 1.0;

    ICARION::physics::ElectricFieldForce force(dom);
    REQUIRE(force.name() == "ElectricField(TIMS)");
}

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/spacecharge/SpaceChargeGPUModel.h"
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"

using Catch::Approx;
using ICARION::core::IonEnsemble;
using ICARION::core::IonState;
using ICARION::core::Vec3;
using ICARION::physics::SpaceChargeGPUModel;

#ifndef ICARION_USE_GPU

TEST_CASE("SpaceChargeGPUModel stub returns zero field", "[spacecharge][gpu][stub]") {
    std::vector<IonState> ions;
    ions.resize(4);
    for (size_t i = 0; i < ions.size(); ++i) {
        ions[i].pos = Vec3{static_cast<double>(i) * 1e-3,
                           static_cast<double>(i) * 2e-3,
                           static_cast<double>(i) * -1e-3};
        ions[i].vel = Vec3{0.0, 0.0, 0.0};
        ions[i].ion_charge_C = (i % 2 == 0) ? 1.0e-19 : -1.0e-19;
        ions[i].mass_kg = 1.0e-26;
        ions[i].active = true;
        ions[i].born = true;
    }

    IonEnsemble ensemble = IonEnsemble::from_legacy(ions);
    SpaceChargeGPUModel model("gpu-stub-test");

    REQUIRE_FALSE(model.is_available());
    Vec3 first_before = model.sample_electric_field(0);
    REQUIRE(first_before.x == Approx(0.0));
    REQUIRE(first_before.y == Approx(0.0));
    REQUIRE(first_before.z == Approx(0.0));

    model.update_fields(ensemble, 0.0);
    for (size_t i = 0; i < ensemble.size(); ++i) {
        Vec3 E = model.sample_electric_field(i);
        REQUIRE(E.x == Approx(0.0));
        REQUIRE(E.y == Approx(0.0));
        REQUIRE(E.z == Approx(0.0));
    }
}

#else

TEST_CASE("SpaceChargeGPUModel compilation placeholder", "[spacecharge][gpu]") {
    SUCCEED("GPU-enabled builds validate runtime behaviour via integration tests.");
}

#endif

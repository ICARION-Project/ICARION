// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/physics/spacecharge/SpaceChargeDirectModel.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include "utils/constants.h"

#include <stdexcept>

using namespace ICARION;
using namespace ICARION::physics;
using namespace ICARION::config;
using ICARION::core::IonEnsemble;
using Catch::Matchers::WithinAbs;

namespace {

IonState make_ion(const core::Vec3& pos, double charge_c) {
    IonState ion;
    ion.pos = pos;
    ion.vel = core::Vec3{0.0, 0.0, 0.0};
    ion.ion_charge_C = charge_c;
    ion.mass_kg = 50.0 * AMU_TO_KG;
    ion.active = true;
    ion.born = true;
    return ion;
}

DomainConfig make_domain() {
    DomainConfig domain;
    domain.instrument = Instrument::IMS;
    domain.geometry.length_m = 0.1;
    domain.geometry.radius_m = 0.01;
    domain.environment.temperature_K = 300.0;
    domain.environment.pressure_Pa = 101325.0;
    domain.environment.gas_species = "N2";
    domain.environment.compute_derived_properties();
    return domain;
}

} // namespace

TEST_CASE("SpaceChargeDirectModel matches Coulomb field for two ions", "[spacecharge][direct-model]") {
    IonState ion0 = make_ion(core::Vec3{0.0, 0.0, 0.0}, ELEM_CHARGE_C);
    IonState ion1 = make_ion(core::Vec3{1e-6, 0.0, 0.0}, -ELEM_CHARGE_C);
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion0, ion1});

    SpaceChargeDirectModel model{0.0};
    model.update_fields(ensemble, 0.0);

    const double r = 1e-6;
    const double expected_scalar = COULOMB_CONST * ion1.ion_charge_C / (r * r);
    const double direction = (ion0.pos.x - ion1.pos.x) / r;  // ±1
    const double expected_ex = expected_scalar * direction;

    auto E0 = model.sample_electric_field(0);
    REQUIRE_THAT(E0.x, WithinAbs(expected_ex, 1e-6 * std::abs(expected_ex)));
    REQUIRE_THAT(E0.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(E0.z, WithinAbs(0.0, 1e-12));

    const double expected_ex1 = COULOMB_CONST * ion0.ion_charge_C / (r * r) * ((ion1.pos.x - ion0.pos.x) / r);
    auto E1 = model.sample_electric_field(1);
    REQUIRE_THAT(E1.x, WithinAbs(expected_ex1, 1e-6 * std::abs(expected_ex1)));
    REQUIRE_THAT(E1.y, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(E1.z, WithinAbs(0.0, 1e-12));
}

TEST_CASE("ForceRegistry uses space-charge model contribution", "[spacecharge][force-registry]") {
    auto domain = make_domain();
    auto registry = std::make_shared<ForceRegistry>(domain);
    auto model = std::make_shared<SpaceChargeDirectModel>(0.0);
    registry->set_space_charge_model(model);

    IonState ion0 = make_ion(core::Vec3{0.0, 0.0, 0.0}, ELEM_CHARGE_C);
    IonState ion1 = make_ion(core::Vec3{0.0, 1e-6, 0.0}, 2 * ELEM_CHARGE_C);
    IonEnsemble ensemble = IonEnsemble::from_legacy({ion0, ion1});

    model->update_fields(ensemble, 0.0);

    ForceContext ctx;
    ctx.ion_ensemble = &ensemble;
    ctx.ion_index = 0;

    auto total_force = registry->compute_total_force(ensemble, 0, 0.0, ctx);

    const double r = 1e-6;
    const double expected_force_mag = COULOMB_CONST * ion0.ion_charge_C * ion1.ion_charge_C / (r * r);
    REQUIRE_THAT(total_force.y, WithinAbs(-expected_force_mag, 1e-6 * expected_force_mag));
    REQUIRE_THAT(total_force.x, WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(total_force.z, WithinAbs(0.0, 1e-12));
}

TEST_CASE("SpaceChargeDirectModel handles edge cases", "[spacecharge][direct-model][edge]") {
    SECTION("Zero ions") {
        IonEnsemble empty = IonEnsemble::from_legacy({});
        SpaceChargeDirectModel model{0.0};

        REQUIRE_NOTHROW(model.update_fields(empty, 0.0));
        auto E = model.sample_electric_field(0);

        REQUIRE_THAT(E.x, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(E.y, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(E.z, WithinAbs(0.0, 1e-30));
    }

    SECTION("Single ion has no self-field") {
        IonState ion = make_ion(core::Vec3{0.0, 0.0, 0.0}, ELEM_CHARGE_C);
        IonEnsemble ensemble = IonEnsemble::from_legacy({ion});
        SpaceChargeDirectModel model{0.0};

        model.update_fields(ensemble, 0.0);
        auto E = model.sample_electric_field(0);

        REQUIRE_THAT(E.x, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(E.y, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(E.z, WithinAbs(0.0, 1e-30));
    }

    SECTION("Inactive ion is excluded") {
        IonState active_ion = make_ion(core::Vec3{0.0, 0.0, 0.0}, ELEM_CHARGE_C);
        IonState inactive_ion = make_ion(core::Vec3{1e-6, 0.0, 0.0}, ELEM_CHARGE_C);
        inactive_ion.active = false;
        IonEnsemble ensemble = IonEnsemble::from_legacy({active_ion, inactive_ion});
        SpaceChargeDirectModel model{0.0};

        model.update_fields(ensemble, 0.0);
        auto E = model.sample_electric_field(0);

        REQUIRE_THAT(E.x, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(E.y, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(E.z, WithinAbs(0.0, 1e-30));
    }

    SECTION("N-body force sum cancels") {
        std::vector<core::Vec3> positions = {
            {0.0, 0.0, 0.0},
            {3e-6, 0.0, 0.0},
            {0.0, 4e-6, 0.0},
            {1e-6, 1e-6, 2e-6}};
        std::vector<IonState> ions;
        for (const auto& pos : positions) {
            ions.push_back(make_ion(pos, ELEM_CHARGE_C));
        }

        IonEnsemble ensemble = IonEnsemble::from_legacy(ions);
        SpaceChargeDirectModel model{0.0};
        model.update_fields(ensemble, 0.0);

        core::Vec3 force_sum{0.0, 0.0, 0.0};
        for (size_t i = 0; i < ensemble.size(); ++i) {
            force_sum += model.sample_electric_field(i) * ELEM_CHARGE_C;
        }

        REQUIRE_THAT(force_sum.x, WithinAbs(0.0, 1e-20));
        REQUIRE_THAT(force_sum.y, WithinAbs(0.0, 1e-20));
        REQUIRE_THAT(force_sum.z, WithinAbs(0.0, 1e-20));
    }

    SECTION("Negative softening is rejected") {
        REQUIRE_THROWS_AS(SpaceChargeDirectModel(-1e-12), std::invalid_argument);
    }

    SECTION("Out-of-bounds sample returns zero") {
        SpaceChargeDirectModel model{0.0};
        auto E = model.sample_electric_field(999);

        REQUIRE_THAT(E.x, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(E.y, WithinAbs(0.0, 1e-30));
        REQUIRE_THAT(E.z, WithinAbs(0.0, 1e-30));
    }
}

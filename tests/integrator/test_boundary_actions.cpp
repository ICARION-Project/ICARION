// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/config/types/BoundaryConfig.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/InstrumentTypes.h"
#include "core/integrator/DomainManager.h"
#include "core/integrator/boundaries/BoundaryActionFactory.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

#include <cmath>
#include <memory>
#include <random>
#include <vector>

using namespace ICARION;
using namespace ICARION::config;
using namespace ICARION::integrator;

namespace {

IonState make_ion() {
    IonState ion;
    ion.pos = Vec3{0.0, 0.0, 0.0};
    ion.vel = Vec3{1.0, -2.0, 3.0};
    ion.mass_kg = 100.0 * 1.66053906660e-27;
    ion.active = true;
    ion.death_time_s = -1.0;
    return ion;
}

DomainConfig make_domain(const std::string& name, BoundaryActionType type) {
    DomainConfig d;
    d.name = name;
    d.instrument = Instrument::IMS;
    d.geometry.origin_m = Vec3{0.0, 0.0, 0.0};
    d.geometry.length_m = 0.1;
    d.geometry.radius_m = 0.01;
    d.boundary.type = type;
    d.finalize();
    return d;
}

}  // namespace

TEST_CASE("BoundaryActionFactory: absorption deactivates ion at boundary", "[integrator][boundary]") {
    BoundaryConfig cfg;
    cfg.type = BoundaryActionType::Absorption;

    auto action = BoundaryActionFactory::create(cfg, nullptr);
    REQUIRE(action != nullptr);

    IonState ion = make_ion();
    const Vec3 normal{0.0, 0.0, 1.0};
    const Vec3 hit{0.001, -0.002, 0.003};
    constexpr double t = 1.25e-6;

    action->apply(ion, normal, hit, 300.0, t);

    REQUIRE_FALSE(ion.active);
    REQUIRE_THAT(ion.pos.x, Catch::Matchers::WithinAbs(hit.x, 1e-15));
    REQUIRE_THAT(ion.pos.y, Catch::Matchers::WithinAbs(hit.y, 1e-15));
    REQUIRE_THAT(ion.pos.z, Catch::Matchers::WithinAbs(hit.z, 1e-15));
    REQUIRE_THAT(ion.vel.x, Catch::Matchers::WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ion.vel.y, Catch::Matchers::WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ion.vel.z, Catch::Matchers::WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(ion.death_time_s, Catch::Matchers::WithinRel(t, 1e-12));
}

TEST_CASE("BoundaryActionFactory: specular reflection mirrors normal component", "[integrator][boundary]") {
    BoundaryConfig cfg;
    cfg.type = BoundaryActionType::SpecularReflection;

    auto action = BoundaryActionFactory::create(cfg, nullptr);
    REQUIRE(action != nullptr);

    IonState ion = make_ion();
    const Vec3 normal{0.0, 0.0, 1.0};
    const Vec3 hit{0.0, 0.0, 0.01};

    action->apply(ion, normal, hit, 300.0, 0.0);

    REQUIRE(ion.active);
    REQUIRE_THAT(ion.pos.z, Catch::Matchers::WithinAbs(hit.z, 1e-15));
    REQUIRE_THAT(ion.vel.x, Catch::Matchers::WithinAbs(1.0, 1e-15));
    REQUIRE_THAT(ion.vel.y, Catch::Matchers::WithinAbs(-2.0, 1e-15));
    REQUIRE_THAT(ion.vel.z, Catch::Matchers::WithinAbs(-3.0, 1e-15));
    REQUIRE_THAT(ion.death_time_s, Catch::Matchers::WithinAbs(-1.0, 1e-15));
}

TEST_CASE("BoundaryActionFactory: diffuse/thermal require RNG when created directly", "[integrator][boundary]") {
    BoundaryConfig diffuse;
    diffuse.type = BoundaryActionType::DiffuseReflection;
    REQUIRE_THROWS_AS(BoundaryActionFactory::create(diffuse, nullptr), std::runtime_error);

    BoundaryConfig thermal;
    thermal.type = BoundaryActionType::ThermalReflection;
    REQUIRE_THROWS_AS(BoundaryActionFactory::create(thermal, nullptr), std::runtime_error);
}

TEST_CASE("BoundaryActionFactory: diffuse and thermal reflection produce inward finite velocities",
          "[integrator][boundary]") {
    std::mt19937 rng(42);
    const Vec3 normal{0.0, 0.0, 1.0};
    const Vec3 hit{0.0, 0.0, 0.01};

    SECTION("Diffuse") {
        BoundaryConfig cfg;
        cfg.type = BoundaryActionType::DiffuseReflection;
        cfg.accommodation_coeff = 1.0;

        auto action = BoundaryActionFactory::create(cfg, &rng);
        IonState ion = make_ion();
        action->apply(ion, normal, hit, 300.0, 0.0);

        const double v_dot_n = ion.vel.x * normal.x + ion.vel.y * normal.y + ion.vel.z * normal.z;
        REQUIRE(ion.active);
        REQUIRE(v_dot_n > 0.0);
        REQUIRE(std::isfinite(ion.vel.x));
        REQUIRE(std::isfinite(ion.vel.y));
        REQUIRE(std::isfinite(ion.vel.z));
    }

    SECTION("Thermal") {
        BoundaryConfig cfg;
        cfg.type = BoundaryActionType::ThermalReflection;

        auto action = BoundaryActionFactory::create(cfg, &rng);
        IonState ion = make_ion();
        action->apply(ion, normal, hit, 300.0, 0.0);

        const double v_dot_n = ion.vel.x * normal.x + ion.vel.y * normal.y + ion.vel.z * normal.z;
        REQUIRE(ion.active);
        REQUIRE(v_dot_n > 0.0);
        REQUIRE(std::isfinite(ion.vel.x));
        REQUIRE(std::isfinite(ion.vel.y));
        REQUIRE(std::isfinite(ion.vel.z));
    }
}

TEST_CASE("BoundaryActionFactory: diffuse/thermal reflection is deterministic per event signature",
          "[integrator][boundary][rng]") {
    std::mt19937 rng(42);
    const Vec3 normal{0.0, 0.0, 1.0};
    const Vec3 hit{1.5e-3, -2.5e-3, 1.0e-2};
    constexpr double t_event = 2.75e-6;

    SECTION("Diffuse deterministic") {
        BoundaryConfig cfg;
        cfg.type = BoundaryActionType::DiffuseReflection;
        cfg.accommodation_coeff = 1.0;

        auto action = BoundaryActionFactory::create(cfg, &rng);
        IonState ion_a = make_ion();
        IonState ion_b = make_ion();
        ion_a.species_id = "H3O+";
        ion_b.species_id = "H3O+";
        ion_a.history_index = 17;
        ion_b.history_index = 17;

        action->apply(ion_a, normal, hit, 300.0, t_event);
        action->apply(ion_b, normal, hit, 300.0, t_event);

        REQUIRE_THAT(ion_a.vel.x, Catch::Matchers::WithinAbs(ion_b.vel.x, 1e-15));
        REQUIRE_THAT(ion_a.vel.y, Catch::Matchers::WithinAbs(ion_b.vel.y, 1e-15));
        REQUIRE_THAT(ion_a.vel.z, Catch::Matchers::WithinAbs(ion_b.vel.z, 1e-15));
    }

    SECTION("Thermal deterministic") {
        BoundaryConfig cfg;
        cfg.type = BoundaryActionType::ThermalReflection;

        auto action = BoundaryActionFactory::create(cfg, &rng);
        IonState ion_a = make_ion();
        IonState ion_b = make_ion();
        ion_a.species_id = "N2+";
        ion_b.species_id = "N2+";
        ion_a.history_index = 8;
        ion_b.history_index = 8;

        action->apply(ion_a, normal, hit, 315.0, t_event);
        action->apply(ion_b, normal, hit, 315.0, t_event);

        REQUIRE_THAT(ion_a.vel.x, Catch::Matchers::WithinAbs(ion_b.vel.x, 1e-15));
        REQUIRE_THAT(ion_a.vel.y, Catch::Matchers::WithinAbs(ion_b.vel.y, 1e-15));
        REQUIRE_THAT(ion_a.vel.z, Catch::Matchers::WithinAbs(ion_b.vel.z, 1e-15));
    }
}

TEST_CASE("DomainManager builds configured boundary actions for each domain", "[integrator][boundary][domain]") {
    std::vector<DomainConfig> domains;
    domains.push_back(make_domain("absorption", BoundaryActionType::Absorption));
    domains.push_back(make_domain("specular", BoundaryActionType::SpecularReflection));
    domains.push_back(make_domain("diffuse", BoundaryActionType::DiffuseReflection));
    domains.push_back(make_domain("thermal", BoundaryActionType::ThermalReflection));

    DomainManager manager(domains, 1234);

    REQUIRE(manager.boundary_action(0) != nullptr);
    REQUIRE(manager.boundary_action(1) != nullptr);
    REQUIRE(manager.boundary_action(2) != nullptr);
    REQUIRE(manager.boundary_action(3) != nullptr);

    REQUIRE(manager.boundary_action(0)->name() == "Absorption");
    REQUIRE(manager.boundary_action(1)->name() == "Specular Reflection");
    REQUIRE(manager.boundary_action(2)->name() == "Diffuse Reflection");
    REQUIRE(manager.boundary_action(3)->name() == "Thermal Reflection");
}

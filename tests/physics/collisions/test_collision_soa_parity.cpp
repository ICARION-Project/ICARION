// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/physics/collisions/EHSSCollisionHandler.h"
#include "core/physics/collisions/OUCollisionHandler.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/IonEnsemble.h"

using namespace ICARION;
using namespace ICARION::physics;
using Catch::Approx;

namespace {
config::EnvironmentConfig make_env() {
    config::EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.particle_density_m_3 = 1e25;          // dense to ensure collision
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.gas_mass_kg = 28.0 * AMU_TO_KG;       // N2-like
    return env;
}

IonState make_ion(const std::string& species = "X+") {
    IonState ion;
    ion.pos = Vec3{0.0, 0.0, 0.0};
    ion.vel = Vec3{100.0, 50.0, -25.0};
    ion.mass_kg = 28.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1e-18;  // large to make collision probability ~1
    ion.species_id = species;
    ion.active = true;
    ion.born = true;
    return ion;
}

void require_parity(auto&& handler_aos, auto&& handler_soa) {
    auto env = make_env();
    const double dt = 1e-6;

    IonState ion_aos = make_ion();
    IonState ion_base = make_ion();

    core::IonEnsemble ensemble = core::IonEnsemble::from_legacy({ion_base});
    auto view = ensemble.collision_data(0);

    PhysicsRng rng1(12345);
    PhysicsRng rng2(12345);

    bool hit_aos = handler_aos.handle_collision(ion_aos, dt, rng1, env);
    bool hit_soa = handler_soa.handle_collision_soa(view, dt, rng2, env);

    REQUIRE(hit_aos == hit_soa);
    if (hit_aos && hit_soa) {
        REQUIRE(ion_aos.vel.x == Approx(view.kin.vel().x));
        REQUIRE(ion_aos.vel.y == Approx(view.kin.vel().y));
        REQUIRE(ion_aos.vel.z == Approx(view.kin.vel().z));
    }
}
} // namespace

TEST_CASE("Collision handlers: AoS vs SoA parity", "[collision][soa][parity]") {
    SECTION("HSS") {
        HSSCollisionHandler handler(false, nullptr);
        require_parity(handler, handler);
    }

    SECTION("OU") {
        OUCollisionHandler handler(1e6 /*gamma*/, true);
        require_parity(handler, handler);
    }
}

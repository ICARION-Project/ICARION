// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "utils/constants.h"

using Catch::Approx;
using namespace ICARION;

TEST_CASE("HSS uses gas-specific CCS map in mixture", "[collision][multigas]") {
    config::SpeciesDatabase db;
    config::SpeciesProperties sp;
    sp.id = "X+";
    sp.mass_amu = 28.0;
    sp.charge = 1;
    sp.CCS_m2 = 1.0e-18;
    sp.ccs_hss_m2["N2"] = 1.0e-18;
    sp.ccs_hss_m2["O2"] = 2.0e-18;  // larger CCS in O2
    db.species[sp.id] = sp;

    physics::HSSCollisionHandler handler(false, &db);
    EhssRng rng(1234);

    IonState ion;
    ion.species_id = "X+";
    ion.mass_kg = sp.mass_amu * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = sp.CCS_m2;
    ion.vel = Vec3{1000.0, 0.0, 0.0};

    config::EnvironmentConfig env;
    env.pressure_Pa = 1000.0;  // lower pressure to keep rates reasonable
    env.temperature_K = 300.0;
    env.gas_mixture = {
        {"N2", 0.75, 1.0e-18, -1.0},
        {"O2", 0.25, 1.0e-18, -1.0}
    };
    env.compute_derived_properties();

    const int trials = 2000;
    int collisions = 0;
    for (int i = 0; i < trials; ++i) {
        IonState ion_copy = ion;
        bool occurred = handler.handle_collision(ion_copy, 1e-7, rng, env);
        if (occurred) {
            collisions++;
        }
    }

    auto stats_map = handler.collisions_by_species();
    double n2 = static_cast<double>(stats_map["N2"]);
    double o2 = static_cast<double>(stats_map["O2"]);
    double total = n2 + o2;

    REQUIRE(collisions > 0);
    REQUIRE(total == Approx(collisions));

    double frac_o2 = (total > 0.0) ? o2 / total : 0.0;
    // Expect weighting ~ f_i * sigma_i: N2 0.75*1 vs O2 0.25*2 → O2 weight ~0.4
    REQUIRE(frac_o2 == Approx(0.4).margin(0.1));
}

TEST_CASE("HSS throws when mixture has no sigma", "[collision][multigas][safety]") {
    physics::HSSCollisionHandler handler(false, nullptr);
    EhssRng rng(42);

    IonState ion;
    ion.mass_kg = 28.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 0.0;  // force missing sigma
    ion.vel = Vec3{1000.0, 0.0, 0.0};

    config::EnvironmentConfig env;
    env.pressure_Pa = 1000.0;
    env.temperature_K = 300.0;
    env.gas_mixture = {{"N2", 1.0, -1.0, -1.0}};
    env.compute_derived_properties();

    REQUIRE_THROWS(handler.handle_collision(ion, 1e-7, rng, env));
}

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"

using Catch::Approx;
using namespace ICARION;

TEST_CASE("MultiGasCollisionHandler selects gases by mole fraction", "[collision][multigas]") {
    physics::HSSCollisionHandler handler(false);
    EhssRng rng(1234);

    IonState ion;
    ion.mass_kg = 28.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 1.0e-18;
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

    double frac_n2 = (total > 0.0) ? n2 / total : 0.0;
    REQUIRE(frac_n2 == Approx(0.75).margin(0.1));
}

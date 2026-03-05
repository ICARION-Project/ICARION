// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/physics/reactions/StochasticReactionHandler.h"
#include "core/config/types/ReactionConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/IonEnsemble.h"
#include "utils/constants.h"

using namespace ICARION;
using namespace ICARION::physics;
using namespace ICARION::config;
using Catch::Approx;

namespace {
ReactionDatabase make_reaction_db() {
    ReactionDatabase db;
    Reaction rxn;
    rxn.id = "A_to_B";
    rxn.reactant = "A+";
    rxn.product = "B+";
    rxn.rate_constant = 1e9;  // large → P≈1
    ReactionOrderTerm term;
    term.species = "Buffer";
    term.exponent = 1;
    term.concentration_m3 = -1.0;  // use env particle density
    rxn.order_terms.push_back(term);
    db.reactions.push_back(rxn);
    return db;
}

SpeciesDatabase make_species_db() {
    SpeciesDatabase db;
    SpeciesProperties a;
    a.mass_amu = 10.0;
    a.charge = 1;
    a.mass_kg = a.mass_amu * AMU_TO_KG;
    a.charge_C = ELEM_CHARGE_C;
    a.CCS_m2 = 10e-20;
    a.mobility_m2Vs = 2e-4;
    db.species["A+"] = a;

    SpeciesProperties b = a;
    b.mass_amu = 20.0;
    b.mass_kg = b.mass_amu * AMU_TO_KG;
    b.CCS_m2 = 20e-20;
    b.mobility_m2Vs = 1.5e-4;
    db.species["B+"] = b;
    return db;
}

EnvironmentConfig make_env() {
    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.particle_density_m_3 = 1e25;
    env.gas_mass_kg = 28.0 * AMU_TO_KG;
    return env;
}
} // namespace

TEST_CASE("StochasticReactionHandler: deterministic parity with fixed seed", "[reaction][soa][parity]") {
    auto db = make_reaction_db();
    auto species_db = make_species_db();
    auto env = make_env();
    const double dt = 1e-6;

    // Two ions so species pool includes both reactant and product
    IonState ion0;
    ion0.species_id = "A+";
    ion0.mass_kg = species_db.get("A+").mass_kg;
    ion0.ion_charge_C = species_db.get("A+").charge_C;
    ion0.CCS_m2 = species_db.get("A+").CCS_m2;
    ion0.reduced_mobility_cm2_Vs = species_db.get("A+").mobility_m2Vs / CM2_TO_M2;
    ion0.active = true; ion0.born = true;

    IonState ion1 = ion0;
    ion1.species_id = "B+";  // ensure product exists in pool

    StochasticReactionHandler handler;

    core::IonEnsemble ensemble1 = core::IonEnsemble::from_legacy({ion0, ion1});
    core::IonEnsemble ensemble2 = core::IonEnsemble::from_legacy({ion0, ion1});

    auto view1 = ensemble1.reaction_data(0);
    auto view2 = ensemble2.reaction_data(0);
    PhysicsRng rng1(777);
    PhysicsRng rng2(777);

    bool reacted1 = handler.handle_reaction(view1, dt, rng1, db, species_db, env);
    bool reacted2 = handler.handle_reaction(view2, dt, rng2, db, species_db, env);

    REQUIRE(reacted1 == reacted2);
    if (reacted1 && reacted2) {
        REQUIRE(view1.species_id() == "B+");
        REQUIRE(view2.species_id() == "B+");
        REQUIRE(view1.kin.get_mass() == Approx(species_db.get("B+").mass_kg));
        REQUIRE(view1.kin.get_charge() == Approx(species_db.get("B+").charge_C));
        REQUIRE(view1.get_CCS() == Approx(species_db.get("B+").CCS_m2));
        REQUIRE(view1.get_mobility() == Approx(species_db.get("B+").mobility_m2Vs / CM2_TO_M2));
    }
}

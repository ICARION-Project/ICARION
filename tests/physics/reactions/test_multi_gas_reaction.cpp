// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <unordered_map>

#include "core/physics/reactions/StochasticReactionHandler.h"
#include "core/config/types/ReactionConfig.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/types/IonEnsemble.h"
#include "core/types/CollisionTypes.h"
#include "utils/constants.h"

using Catch::Approx;
using namespace ICARION;
using ICARION::physics::PhysicsRng;

TEST_CASE("StochasticReactionHandler uses mixture partial densities", "[reaction][multigas]") {
    // Species DB
    config::SpeciesDatabase db;
    auto make_species = [](const std::string& id, double mass_amu, int charge) {
        config::SpeciesProperties p;
        p.id = id;
        p.mass_amu = mass_amu;
        p.charge = charge;
        p.CCS_A2 = 100.0;
        p.mobility_cm2Vs = 1.0;
        p.convert_to_SI();
        return p;
    };
    db.species["A+"] = make_species("A+", 28.0, 1);
    db.species["B+"] = make_species("B+", 30.0, 1);
    db.species["C+"] = make_species("C+", 32.0, 1);

    // Reactions: A+ + N2 vs A+ + O2
    config::Reaction rxn1;
    rxn1.id = "A_to_B";
    rxn1.reactant = "A+";
    rxn1.product = "B+";
    rxn1.rate_constant = 1.0e-10;  // m3/s
    rxn1.order_terms.push_back({"N2", 1, -1.0});

    config::Reaction rxn2;
    rxn2.id = "A_to_C";
    rxn2.reactant = "A+";
    rxn2.product = "C+";
    rxn2.rate_constant = 5.0e-10;  // m3/s (5x faster)
    rxn2.order_terms.push_back({"O2", 1, -1.0});

    config::ReactionDatabase rxn_db;
    rxn_db.reactions.push_back(rxn1);
    rxn_db.reactions.push_back(rxn2);

    config::EnvironmentConfig env;
    env.pressure_Pa = 1000.0;
    env.temperature_K = 300.0;
    env.gas_mixture = {
        {"N2", 0.75, -1.0, -1.0},
        {"O2", 0.25, -1.0, -1.0}
    };
    env.compute_derived_properties();

    physics::StochasticReactionHandler handler(false);
    const int trials = 500;
    int count_B = 0;
    int count_C = 0;

    IonState ion_A;
    ion_A.species_id = "A+";
    ion_A.mass_kg = db.get("A+").mass_kg;
    ion_A.ion_charge_C = db.get("A+").charge_C;
    ion_A.CCS_m2 = db.get("A+").CCS_m2;
    ion_A.reduced_mobility_cm2_Vs = db.get("A+").mobility_m2Vs / CM2_TO_M2;
    ion_A.active = true; ion_A.born = true;

    IonState ion_B = ion_A; ion_B.species_id = "B+";
    IonState ion_C = ion_A; ion_C.species_id = "C+";

    for (int i = 0; i < trials; ++i) {
        // Ensure species pool contains reactant and products
        core::IonEnsemble ensemble = core::IonEnsemble::from_legacy({ion_A, ion_B, ion_C});
        auto view = ensemble.reaction_data(0);

        PhysicsRng rng(static_cast<uint64_t>(42 + i));
        bool reacted = handler.handle_reaction(view, 1e-7, rng, rxn_db, db, env);
        if (reacted) {
            const std::string product = view.species_id();
            if (product == "B+") count_B++;
            if (product == "C+") count_C++;
        }
    }

    int total = count_B + count_C;
    REQUIRE(total > 0);
    double frac_B = static_cast<double>(count_B) / static_cast<double>(total);
    // Expected weighting ~ (k1*nN2) / (k1*nN2 + k2*nO2) with k1=1, k2=5, n fractions 0.75/0.25
    double nN2 = env.gas_mixture[0].density_m3;
    double nO2 = env.gas_mixture[1].density_m3;
    double expected = (1.0 * nN2) / (1.0 * nN2 + 5.0 * nO2);
    REQUIRE(frac_B == Approx(expected).margin(0.1));
}

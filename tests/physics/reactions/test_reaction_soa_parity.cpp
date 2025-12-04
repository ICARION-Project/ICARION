// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

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

TEST_CASE("StochasticReactionHandler: AoS vs SoA parity", "[reaction][soa][parity]") {
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

    // AoS path
    StochasticReactionHandler handler;
    PhysicsRng rng1(777);
    IonState ion_aos = ion0;
    bool reacted_aos = handler.handle_reaction(ion_aos, dt, rng1, db, species_db, env);

    // SoA path
    core::IonEnsemble ensemble = core::IonEnsemble::from_legacy({ion0, ion1});
    auto view = ensemble.reaction_data(0);
    double* CCS = ensemble.CCS_data();
    double* mobility = ensemble.mobility_data();
    PhysicsRng rng2(777);
    bool reacted_soa = handler.handle_reaction_soa(view, CCS, mobility, dt, rng2, db, species_db, env);

    REQUIRE(reacted_aos == reacted_soa);
    if (reacted_aos && reacted_soa) {
        // Species should be product
        REQUIRE(ion_aos.species_id == "B+");
        auto pool = ensemble.species_pool();
        REQUIRE(pool->end() != std::find(pool->begin(), pool->end(), "B+"));
        // Check mass/CCS/mobility updated
        REQUIRE(ion_aos.mass_kg == Approx(species_db.get("B+").mass_kg));
        REQUIRE(CCS[0] == Approx(species_db.get("B+").CCS_m2));
        REQUIRE(mobility[0] == Approx(species_db.get("B+").mobility_m2Vs / CM2_TO_M2));
    }
}


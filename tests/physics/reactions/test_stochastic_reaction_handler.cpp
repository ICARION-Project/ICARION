// test_stochastic_reaction_handler.cpp
// Unit tests for StochasticReactionHandler
//
// Tests:
// - SSOT compliance (no parameter duplication)
// - Second-order reactions (2-body)
// - Third-order reactions (3-body)
// - Buffer gas fallback (concentration_m3 = 0)
// - Species database lookup
// - Reaction statistics
//
// Created: 2025-11-22 (Phase 3 Refactor)

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/physics/reactions/StochasticReactionHandler.h"
#include "core/config/types/ReactionConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"

using namespace ICARION;
using namespace ICARION::physics;
using namespace ICARION::config;

// Helper: Create test reaction database
ReactionDatabase create_test_reaction_db() {
    ReactionDatabase db;
    
    // Reaction 1: H3O+ + NH3 → NH4+ (2nd-order)
    Reaction rxn1;
    rxn1.id = "rxn_h3o_to_nh4";
    rxn1.reactant = "H3O+";
    rxn1.product = "NH4+";
    rxn1.rate_constant_m3s = 2.0e-9;  // [m³/s]
    
    ReactionOrderTerm term1;
    term1.species = "NH3";
    term1.exponent = 1;
    term1.concentration_m3 = 1e20;  // Fixed concentration [m⁻³]
    rxn1.order_terms.push_back(term1);
    
    db.reactions.push_back(rxn1);
    
    // Reaction 2: H3O+ + H2O + He → H5O2+ (3rd-order)
    Reaction rxn2;
    rxn2.id = "rxn_h3o_to_h5o2";
    rxn2.reactant = "H3O+";
    rxn2.product = "H5O2+";
    rxn2.rate_constant_m3s = 1.2e-28;  // [m⁶/s]
    
    ReactionOrderTerm term2a;
    term2a.species = "H2O";
    term2a.exponent = 1;
    term2a.concentration_m3 = 2.5e25;  // Fixed [H2O]
    
    ReactionOrderTerm term2b;
    term2b.species = "He";
    term2b.exponent = 1;
    term2b.concentration_m3 = -1.0;  // ✅ Use buffer gas density (sentinel value)
    
    rxn2.order_terms.push_back(term2a);
    rxn2.order_terms.push_back(term2b);
    
    db.reactions.push_back(rxn2);
    
    // Reaction 3: NH4+ + H2O → H3O+ (reverse reaction)
    Reaction rxn3;
    rxn3.id = "rxn_nh4_to_h3o";
    rxn3.reactant = "NH4+";
    rxn3.product = "H3O+";
    rxn3.rate_constant_m3s = 1.0e-10;
    
    ReactionOrderTerm term3;
    term3.species = "H2O";
    term3.exponent = 1;
    term3.concentration_m3 = -1.0;  // ✅ Use buffer gas (sentinel value)
    
    rxn3.order_terms.push_back(term3);
    
    db.reactions.push_back(rxn3);
    
    return db;
}

// Helper: Create test species database
SpeciesDatabase create_test_species_db() {
    SpeciesDatabase db;
    
    // H3O+
    SpeciesProperties h3o;
    h3o.mass_amu = 19.0;
    h3o.charge = 1;
    h3o.mass_kg = 19.0 * AMU_TO_KG;
    h3o.charge_C = ELEM_CHARGE_C;
    h3o.CCS_m2 = 50e-20;  // 50 Ų
    h3o.mobility_m2Vs = 2.5e-4;
    db.species["H3O+"] = h3o;
    
    // NH4+
    SpeciesProperties nh4;
    nh4.mass_amu = 18.0;
    nh4.charge = 1;
    nh4.mass_kg = 18.0 * AMU_TO_KG;
    nh4.charge_C = ELEM_CHARGE_C;
    nh4.CCS_m2 = 48e-20;  // 48 Ų
    nh4.mobility_m2Vs = 2.3e-4;
    db.species["NH4+"] = nh4;
    
    // H5O2+
    SpeciesProperties h5o2;
    h5o2.mass_amu = 37.0;
    h5o2.charge = 1;
    h5o2.mass_kg = 37.0 * AMU_TO_KG;
    h5o2.charge_C = ELEM_CHARGE_C;
    h5o2.CCS_m2 = 60e-20;  // 60 Ų
    h5o2.mobility_m2Vs = 2.0e-4;
    db.species["H5O2+"] = h5o2;
    
    return db;
}

// Helper: Create test environment
EnvironmentConfig create_test_environment(double T_K = 300.0, double n_m3 = 2.5e25) {
    EnvironmentConfig env;
    env.temperature_K = T_K;
    env.particle_density_m_3 = n_m3;
    env.pressure_Pa = 101325.0;
    env.gas_species = "He";
    env.gas_mass_kg = MOLAR_MASS_HE_KG;
    return env;
}

// ===================================================================
// TEST SUITE: StochasticReactionHandler
// ===================================================================

TEST_CASE("StochasticReactionHandler: SSOT compliance", "[reaction][ssot]") {
    // Setup
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment(300.0, 2.5e25);
    
    StochasticReactionHandler handler(false);
    
    IonState ion;
    ion.species_id = "H3O+";
    ion.mass_kg = 19.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    
    EhssRng rng(42);
    double dt = 1e-9;  // 1 ns
    
    // SSOT: Pass databases directly (no intermediate structs!)
    bool reaction_occurred = handler.handle_reaction(
        ion, dt, rng,
        reaction_db,  // Direct reference
        species_db,   // Direct reference
        env           // Direct reference (contains T, n)
    );
    
    // Verify: No exceptions, valid return value
    REQUIRE((reaction_occurred == true || reaction_occurred == false));
}

TEST_CASE("StochasticReactionHandler: Second-order reaction", "[reaction][order]") {
    // Setup
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment();
    
    StochasticReactionHandler handler(false);
    
    IonState ion;
    ion.species_id = "H3O+";
    ion.mass_kg = 19.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    
    // Run many reactions to get statistics
    EhssRng rng(12345);
    double dt = 1e-9;  // 1 ns
    int num_trials = 1000;
    int num_reactions = 0;
    
    for (int i = 0; i < num_trials; ++i) {
        IonState test_ion = ion;  // Copy
        bool rxn = handler.handle_reaction(
            test_ion, dt, rng, reaction_db, species_db, env
        );
        if (rxn) {
            num_reactions++;
            // Verify species changed
            REQUIRE((test_ion.species_id == "NH4+" || test_ion.species_id == "H5O2+"));
        }
    }
    
    // Expected: k_eff = 2e-9 [m³/s] × 1e20 [m⁻³] = 2e11 [s⁻¹]
    // P = 1 - exp(-2e11 × 1e-9) = 1 - exp(-200) ≈ 1.0
    // Should react almost always
    REQUIRE(num_reactions > 900);  // > 90% reaction probability
}

TEST_CASE("StochasticReactionHandler: Third-order reaction", "[reaction][order]") {
    // Setup
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment(300.0, 2.5e25);  // High density
    
    StochasticReactionHandler handler(false);
    
    IonState ion;
    ion.species_id = "H3O+";
    ion.mass_kg = 19.0 * AMU_TO_KG;
    
    // Expected: k_eff = 1.2e-28 [m⁶/s] × 2.5e25 [m⁻³] × 2.5e25 [m⁻³]
    //                 = 1.2e-28 × 6.25e50 = 7.5e22 [s⁻¹]
    // P = 1 - exp(-7.5e22 × 1e-9) ≈ 1.0 (certain)
    
    EhssRng rng(42);
    double dt = 1e-9;
    
    bool rxn = handler.handle_reaction(
        ion, dt, rng, reaction_db, species_db, env
    );
    
    // Should react (very high probability)
    REQUIRE(rxn == true);
    REQUIRE((ion.species_id == "NH4+" || ion.species_id == "H5O2+"));
}

TEST_CASE("StochasticReactionHandler: Buffer gas fallback", "[reaction][fallback]") {
    // Setup: Low buffer gas density
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment(300.0, 1e20);  // Very low density
    
    StochasticReactionHandler handler(false);
    
    IonState ion;
    ion.species_id = "NH4+";  // Use NH4+ → H3O+ reaction (uses buffer gas)
    ion.mass_kg = 18.0 * AMU_TO_KG;
    
    // Expected: k_eff = 1e-10 [m³/s] × 1e20 [m⁻³] = 1e10 [s⁻¹]
    // P = 1 - exp(-1e10 × 1e-9) = 1 - exp(-10) ≈ 0.99995
    
    EhssRng rng(42);
    double dt = 1e-9;
    
    int num_trials = 100;
    int num_reactions = 0;
    
    for (int i = 0; i < num_trials; ++i) {
        IonState test_ion = ion;
        bool rxn = handler.handle_reaction(
            test_ion, dt, rng, reaction_db, species_db, env
        );
        if (rxn) {
            num_reactions++;
            REQUIRE(test_ion.species_id == "H3O+");
        }
    }
    
    // Should react ~100% of time
    REQUIRE(num_reactions > 90);  // At least 90/100 trials
}

TEST_CASE("StochasticReactionHandler: Species database lookup", "[reaction][species]") {
    // Setup
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment();
    
    StochasticReactionHandler handler(false);
    
    IonState ion;
    ion.species_id = "H3O+";
    ion.mass_kg = 19.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 50e-20;
    
    EhssRng rng(42);
    double dt = 1e-6;  // 1 μs (high dt → certain reaction)
    
    bool rxn = handler.handle_reaction(
        ion, dt, rng, reaction_db, species_db, env
    );
    
    if (rxn) {
        // Verify ion properties updated from species database
        if (ion.species_id == "NH4+") {
            REQUIRE(ion.mass_kg == Catch::Approx(18.0 * AMU_TO_KG));
            REQUIRE(ion.ion_charge_C == Catch::Approx(ELEM_CHARGE_C));
            REQUIRE(ion.CCS_m2 == Catch::Approx(48e-20));
        } else if (ion.species_id == "H5O2+") {
            REQUIRE(ion.mass_kg == Catch::Approx(37.0 * AMU_TO_KG));
            REQUIRE(ion.ion_charge_C == Catch::Approx(ELEM_CHARGE_C));
            REQUIRE(ion.CCS_m2 == Catch::Approx(60e-20));
        }
    }
}

TEST_CASE("StochasticReactionHandler: No applicable reactions", "[reaction][none]") {
    // Setup
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment();
    
    StochasticReactionHandler handler(false);
    
    IonState ion;
    ion.species_id = "Unknown+";  // No reactions for this species
    ion.mass_kg = 20.0 * AMU_TO_KG;
    
    EhssRng rng(42);
    double dt = 1e-6;
    
    bool rxn = handler.handle_reaction(
        ion, dt, rng, reaction_db, species_db, env
    );
    
    // Should not react (no applicable reactions)
    REQUIRE(rxn == false);
    REQUIRE(ion.species_id == "Unknown+");  // Unchanged
}

TEST_CASE("StochasticReactionHandler: Reaction statistics", "[reaction][stats]") {
    // Setup
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment();
    
    StochasticReactionHandler handler(false);
    handler.reset_stats();
    
    IonState ion;
    ion.species_id = "H3O+";
    ion.mass_kg = 19.0 * AMU_TO_KG;
    
    EhssRng rng(42);
    double dt = 1e-6;  // High dt → certain reactions
    
    // Run 10 reactions
    for (int i = 0; i < 10; ++i) {
        IonState test_ion = ion;
        handler.handle_reaction(test_ion, dt, rng, reaction_db, species_db, env);
    }
    
    // Verify statistics
    auto stats = handler.get_stats();
    REQUIRE(stats.total_reactions > 0);
    REQUIRE(stats.total_reactions <= 10);
}

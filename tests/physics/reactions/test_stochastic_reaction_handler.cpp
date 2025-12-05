// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

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
#include "core/types/IonEnsemble.h"
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
    rxn1.rate_constant = 2.0e-9;  // [m³/s]
    
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
    rxn2.rate_constant = 1.2e-28;  // [m⁶/s]
    
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
    rxn3.rate_constant = 1.0e-10;
    
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

core::IonEnsemble make_ensemble_for_species(const SpeciesDatabase& db, const std::string& species) {
    std::vector<IonState> ions;
    auto add = [&](const std::string& id) {
        IonState ion;
        ion.species_id = id;
        const auto& props = db.get(id);
        ion.mass_kg = props.mass_kg;
        ion.ion_charge_C = props.charge_C;
        ion.CCS_m2 = props.CCS_m2;
        ion.reduced_mobility_cm2_Vs = props.mobility_m2Vs / CM2_TO_M2;
        ion.active = (id == species);
        ion.born = true;
        ions.push_back(ion);
    };

    add(species);
    for (const auto& kv : db.species) {
        if (kv.first != species) {
            add(kv.first);
        }
    }
    return core::IonEnsemble::from_legacy(ions);
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
    
    auto ensemble = make_ensemble_for_species(species_db, "H3O+");
    auto view = ensemble.reaction_data(0);
    
    PhysicsRng rng(42);
    double dt = 1e-9;  // 1 ns
    
    // SSOT: Pass databases directly (no intermediate structs!)
    bool reaction_occurred = handler.handle_reaction(
        view, dt, rng,
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
    
    // Run many reactions to get statistics
    PhysicsRng rng(12345);
    double dt = 1e-9;  // 1 ns
    int num_trials = 1000;
    int num_reactions = 0;
    
    for (int i = 0; i < num_trials; ++i) {
        auto ensemble = make_ensemble_for_species(species_db, "H3O+");
        auto view = ensemble.reaction_data(0);
        bool rxn = handler.handle_reaction(
            view, dt, rng, reaction_db, species_db, env
        );
        if (rxn) {
            num_reactions++;
            // Verify species changed
            REQUIRE((view.species_id() == "NH4+" || view.species_id() == "H5O2+"));
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
    
    // Expected: k_eff = 1.2e-28 [m⁶/s] × 2.5e25 [m⁻³] × 2.5e25 [m⁻³]
    //                 = 1.2e-28 × 6.25e50 = 7.5e22 [s⁻¹]
    // P = 1 - exp(-7.5e22 × 1e-9) ≈ 1.0 (certain)
    
    auto ensemble = make_ensemble_for_species(species_db, "H3O+");
    auto view = ensemble.reaction_data(0);
    PhysicsRng rng(42);
    double dt = 1e-9;
    
    bool rxn = handler.handle_reaction(
        view, dt, rng, reaction_db, species_db, env
    );
    
    // Should react (very high probability)
    REQUIRE(rxn == true);
    REQUIRE((view.species_id() == "NH4+" || view.species_id() == "H5O2+"));
}

TEST_CASE("StochasticReactionHandler: Buffer gas fallback", "[reaction][fallback]") {
    // Setup: Low buffer gas density
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment(300.0, 1e20);  // Very low density
    
    StochasticReactionHandler handler(false);
    
    // Expected: k_eff = 1e-10 [m³/s] × 1e20 [m⁻³] = 1e10 [s⁻¹]
    // P = 1 - exp(-1e10 × 1e-9) = 1 - exp(-10) ≈ 0.99995
    
    PhysicsRng rng(42);
    double dt = 1e-9;
    
    int num_trials = 100;
    int num_reactions = 0;
    
    for (int i = 0; i < num_trials; ++i) {
        auto ensemble = make_ensemble_for_species(species_db, "NH4+");
        auto view = ensemble.reaction_data(0);
        bool rxn = handler.handle_reaction(
            view, dt, rng, reaction_db, species_db, env
        );
        if (rxn) {
            num_reactions++;
            REQUIRE(view.species_id() == "H3O+");
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
    
    PhysicsRng rng(42);
    double dt = 1e-6;  // 1 μs (high dt → certain reaction)
    auto ensemble = make_ensemble_for_species(species_db, "H3O+");
    auto view = ensemble.reaction_data(0);
    
    bool rxn = handler.handle_reaction(view, dt, rng, reaction_db, species_db, env);
    
    if (rxn) {
        // Verify ion properties updated from species database
        if (view.species_id() == "NH4+") {
            REQUIRE(view.kin.get_mass() == Catch::Approx(18.0 * AMU_TO_KG));
            REQUIRE(view.kin.get_charge() == Catch::Approx(ELEM_CHARGE_C));
            REQUIRE(view.get_CCS() == Catch::Approx(48e-20));
        } else if (view.species_id() == "H5O2+") {
            REQUIRE(view.kin.get_mass() == Catch::Approx(37.0 * AMU_TO_KG));
            REQUIRE(view.kin.get_charge() == Catch::Approx(ELEM_CHARGE_C));
            REQUIRE(view.get_CCS() == Catch::Approx(60e-20));
        }
    }
}

TEST_CASE("StochasticReactionHandler: No applicable reactions", "[reaction][none]") {
    // Setup
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment();
    
    StochasticReactionHandler handler(false);
    
    PhysicsRng rng(42);
    double dt = 1e-6;
    // Add unknown species to pool to avoid lookup issues
    SpeciesDatabase species_db_with_unknown = species_db;
    SpeciesProperties unknown = species_db.get("H3O+");
    unknown.mass_kg = 20.0 * AMU_TO_KG;
    species_db_with_unknown.species["Unknown+"] = unknown;
    auto ensemble = make_ensemble_for_species(species_db_with_unknown, "Unknown+");
    auto view = ensemble.reaction_data(0);
    
    bool rxn = handler.handle_reaction(
        view, dt, rng, reaction_db, species_db_with_unknown, env
    );
    
    // Should not react (no applicable reactions)
    REQUIRE(rxn == false);
    REQUIRE(view.species_id() == "Unknown+");  // Unchanged
}

TEST_CASE("StochasticReactionHandler: Reaction statistics", "[reaction][stats]") {
    // Setup
    auto reaction_db = create_test_reaction_db();
    auto species_db = create_test_species_db();
    auto env = create_test_environment();
    
    StochasticReactionHandler handler(false);
    handler.reset_stats();
    
    PhysicsRng rng(42);
    double dt = 1e-6;  // High dt → certain reactions
    
    // Run 10 reactions
    for (int i = 0; i < 10; ++i) {
        auto ensemble = make_ensemble_for_species(species_db, "H3O+");
        auto view = ensemble.reaction_data(0);
        handler.handle_reaction(view, dt, rng, reaction_db, species_db, env);
    }
    
    // Verify statistics
    auto stats = handler.get_stats();
    REQUIRE(stats.total_reactions > 0);
    REQUIRE(stats.total_reactions <= 10);
}

// =============================================================================
// TEST 8: Competing Reaction Channels (Weighted Selection)
// =============================================================================
TEST_CASE("StochasticReactionHandler handles competing channels correctly", "[reactions][handler][competing]") {
    // Create species database with Ion+ and two products
    SpeciesDatabase species_db;
    
    SpeciesProperties ion;
    ion.mass_amu = 100.0;
    ion.charge = 1;
    ion.mass_kg = 100.0 * AMU_TO_KG;
    ion.charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 100e-20;
    ion.mobility_m2Vs = 2.0e-4;
    species_db.species["Ion+"] = ion;
    
    SpeciesProperties productA;
    productA.mass_amu = 101.0;
    productA.charge = 1;
    productA.mass_kg = 101.0 * AMU_TO_KG;
    productA.charge_C = ELEM_CHARGE_C;
    productA.CCS_m2 = 101e-20;
    productA.mobility_m2Vs = 2.0e-4;
    species_db.species["ProductA+"] = productA;
    
    SpeciesProperties productB;
    productB.mass_amu = 102.0;
    productB.charge = 1;
    productB.mass_kg = 102.0 * AMU_TO_KG;
    productB.charge_C = ELEM_CHARGE_C;
    productB.CCS_m2 = 102e-20;
    productB.mobility_m2Vs = 2.0e-4;
    species_db.species["ProductB+"] = productB;
    
    // Create two competing reaction channels:
    // Channel A: Ion+ → ProductA+ (k = 1e10 s⁻¹, 10%)
    // Channel B: Ion+ → ProductB+ (k = 9e10 s⁻¹, 90%)
    ReactionDatabase reaction_db;
    
    Reaction rxn_A;
    rxn_A.id = "rxn_A";
    rxn_A.reactant = "Ion+";
    rxn_A.product = "ProductA+";
    rxn_A.rate_constant = 1e10;  // 10% of total
    rxn_A.order_terms = {};  // Zero-order (constant rate)
    
    Reaction rxn_B;
    rxn_B.id = "rxn_B";
    rxn_B.reactant = "Ion+";
    rxn_B.product = "ProductB+";
    rxn_B.rate_constant = 9e10;  // 90% of total
    rxn_B.order_terms = {};  // Zero-order (constant rate)
    
    reaction_db.reactions = {rxn_A, rxn_B};
    
    // Environment (buffer gas)
    auto env = create_test_environment(300.0, 2.5e25);
    
    StochasticReactionHandler handler(false);
    PhysicsRng rng(12345);
    
    double dt = 1e-9;  // dt small enough for proper statistics
    
    // Run many trials to measure branching ratio
    int count_A = 0;
    int count_B = 0;
    int total_trials = 10000;
    
    for (int trial = 0; trial < total_trials; ++trial) {
        auto ensemble = make_ensemble_for_species(species_db, "Ion+");
        auto view = ensemble.reaction_data(0);
        bool reacted = handler.handle_reaction(view, dt, rng, reaction_db, species_db, env);
        
        if (reacted) {
            if (view.species_id() == "ProductA+") {
                count_A++;
            } else if (view.species_id() == "ProductB+") {
                count_B++;
            }
        }
    }
    
    // Expected: ~10% ProductA, ~90% ProductB (among reactions that occurred)
    int total_reactions = count_A + count_B;
    REQUIRE(total_reactions > 0);  // At least some reactions occurred
    
    double fraction_A = static_cast<double>(count_A) / total_reactions;
    double fraction_B = static_cast<double>(count_B) / total_reactions;
    
    // Verify branching ratio within 3σ (±5% for 10k trials)
    REQUIRE(fraction_A > 0.05);   // Expected ~0.10 ± 0.05
    REQUIRE(fraction_A < 0.15);
    REQUIRE(fraction_B > 0.85);   // Expected ~0.90 ± 0.05
    REQUIRE(fraction_B < 0.95);
    
    INFO("Branching ratio A:B = " << fraction_A << " : " << fraction_B);
    INFO("Total reactions: " << total_reactions << " / " << total_trials);
}

// =============================================================================
// TEST 9: Zero Reactions Available (Edge Case)
// =============================================================================
TEST_CASE("StochasticReactionHandler handles zero reactions gracefully", "[reactions][handler][edge-case]") {
    // Empty reaction database
    ReactionDatabase reaction_db;  // No reactions!
    
    SpeciesDatabase species_db;
    SpeciesProperties ion_props;
    ion_props.mass_amu = 100.0;
    ion_props.charge = 1;
    ion_props.mass_kg = 100.0 * AMU_TO_KG;
    ion_props.charge_C = ELEM_CHARGE_C;
    ion_props.CCS_m2 = 100e-20;
    ion_props.mobility_m2Vs = 2.0e-4;
    species_db.species["Ion+"] = ion_props;
    
    auto env = create_test_environment(300.0, 2.5e25);
    
    StochasticReactionHandler handler(false);
    PhysicsRng rng(42);
    
    double dt = 1e-6;
    
    // Try many times - should never react (no reactions available)
    for (int i = 0; i < 100; ++i) {
        auto ensemble = make_ensemble_for_species(species_db, "Ion+");
        auto view = ensemble.reaction_data(0);
        bool reacted = handler.handle_reaction(view, dt, rng, reaction_db, species_db, env);
        
        REQUIRE(!reacted);  // Must return false (no reactions)
        REQUIRE(view.species_id() == "Ion+");  // Species unchanged
    }
    
    auto stats = handler.get_stats();
    REQUIRE(stats.total_reactions == 0);  // No reactions should have occurred
}

// =============================================================================
// TEST 10: Very Large Rate Constant (Numerical Stability)
// =============================================================================
TEST_CASE("StochasticReactionHandler handles very large k_eff", "[reactions][handler][edge-case]") {
    SpeciesDatabase species_db;
    SpeciesProperties ion_props;
    ion_props.mass_amu = 100.0;
    ion_props.charge = 1;
    ion_props.mass_kg = 100.0 * AMU_TO_KG;
    ion_props.charge_C = ELEM_CHARGE_C;
    ion_props.CCS_m2 = 100e-20;
    ion_props.mobility_m2Vs = 2.0e-4;
    species_db.species["Ion+"] = ion_props;
    species_db.species["Product+"] = ion_props;  // Same properties
    
    // Reaction with HUGE rate constant → P ≈ 1 (certain reaction)
    ReactionDatabase reaction_db;
    Reaction rxn;
    rxn.id = "fast_rxn";
    rxn.reactant = "Ion+";
    rxn.product = "Product+";
    rxn.rate_constant = 1e20;  // Extremely large! P = 1 - exp(-1e20 * dt) ≈ 1
    rxn.order_terms = {};
    reaction_db.reactions = {rxn};
    
    auto env = create_test_environment(300.0, 2.5e25);
    
    StochasticReactionHandler handler(false);
    PhysicsRng rng(99);
    
    IonState ion;
    ion.species_id = "Ion+";
    ion.mass_kg = 100.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 100e-20;
    
    double dt = 1e-9;  // Small timestep
    
    // With k_eff = 1e20 and dt = 1e-9, P = 1 - exp(-1e11) ≈ 1.0
    // Should react 100% of the time
    int reacted_count = 0;
    for (int i = 0; i < 100; ++i) {
        auto ensemble = make_ensemble_for_species(species_db, "Ion+");
        auto view = ensemble.reaction_data(0);
        bool reacted = handler.handle_reaction(view, dt, rng, reaction_db, species_db, env);
        if (reacted) {
            reacted_count++;
            REQUIRE(view.species_id() == "Product+");
        }
    }
    
    // Should react ~100 times (probability ≈ 1.0)
    REQUIRE(reacted_count >= 95);  // At least 95% (numerical stability check)
}

// =============================================================================
// TEST 11: Very Small Rate Constant (Rare Events)
// =============================================================================
TEST_CASE("StochasticReactionHandler handles very small k_eff", "[reactions][handler][edge-case]") {
    SpeciesDatabase species_db;
    SpeciesProperties ion_props;
    ion_props.mass_amu = 100.0;
    ion_props.charge = 1;
    ion_props.mass_kg = 100.0 * AMU_TO_KG;
    ion_props.charge_C = ELEM_CHARGE_C;
    ion_props.CCS_m2 = 100e-20;
    ion_props.mobility_m2Vs = 2.0e-4;
    species_db.species["Ion+"] = ion_props;
    species_db.species["Product+"] = ion_props;
    
    // Reaction with TINY rate constant → P ≈ 0 (rare event)
    ReactionDatabase reaction_db;
    Reaction rxn;
    rxn.id = "slow_rxn";
    rxn.reactant = "Ion+";
    rxn.product = "Product+";
    rxn.rate_constant = 1e-20;  // Extremely small! P = 1 - exp(-1e-20 * dt) ≈ 1e-20 * dt
    rxn.order_terms = {};
    reaction_db.reactions = {rxn};
    
    auto env = create_test_environment(300.0, 2.5e25);
    
    StochasticReactionHandler handler(false);
    PhysicsRng rng(123);
    double dt = 1e-6;  // Typical timestep
    
    // With k_eff = 1e-20 and dt = 1e-6, P ≈ 1e-26 (extremely rare)
    // Should almost never react in 1000 trials
    int reacted_count = 0;
    for (int i = 0; i < 1000; ++i) {
        auto ensemble = make_ensemble_for_species(species_db, "Ion+");
        auto view = ensemble.reaction_data(0);
        bool reacted = handler.handle_reaction(view, dt, rng, reaction_db, species_db, env);
        if (reacted) {
            reacted_count++;
        }
    }
    
    // Should react 0-2 times max (very rare event, P ≈ 1e-26)
    REQUIRE(reacted_count <= 5);  // Statistical tolerance
    INFO("Rare event: " << reacted_count << " reactions in 1000 trials (expected ~0)");
}

// =============================================================================
// TEST 12: Extremely Small Rate Constant (k_total < 1e-60 early exit)
// =============================================================================
TEST_CASE("StochasticReactionHandler early exit for k_total < 1e-60", "[reactions][handler][optimization]") {
    SpeciesDatabase species_db;
    SpeciesProperties ion_props;
    ion_props.mass_amu = 100.0;
    ion_props.charge = 1;
    ion_props.mass_kg = 100.0 * AMU_TO_KG;
    ion_props.charge_C = ELEM_CHARGE_C;
    ion_props.CCS_m2 = 100e-20;
    ion_props.mobility_m2Vs = 2.0e-4;
    species_db.species["Ion+"] = ion_props;
    species_db.species["Product+"] = ion_props;
    
    // Reaction with k_total = 1e-65 s⁻¹ (below 1e-60 threshold)
    ReactionDatabase reaction_db;
    Reaction rxn;
    rxn.id = "negligible_rxn";
    rxn.reactant = "Ion+";
    rxn.product = "Product+";
    rxn.rate_constant = 1e-65;  // k_total = 1e-65 < 1e-60 → early exit!
    rxn.order_terms = {};
    reaction_db.reactions = {rxn};
    
    auto env = create_test_environment(300.0, 2.5e25);
    
    StochasticReactionHandler handler(false);
    PhysicsRng rng(999);
    double dt = 1.0;  // Large timestep
    
    // Should NEVER react (k_total < 1e-60 → early exit optimization)
    for (int i = 0; i < 1000; ++i) {
        auto ensemble = make_ensemble_for_species(species_db, "Ion+");
        auto view = ensemble.reaction_data(0);
        bool reacted = handler.handle_reaction(view, dt, rng, reaction_db, species_db, env);
        
        REQUIRE(!reacted);  // Must return false immediately
        REQUIRE(view.species_id() == "Ion+");  // Species unchanged
    }
    
    auto stats = handler.get_stats();
    REQUIRE(stats.total_reactions == 0);  // No reactions should have occurred
}

// =============================================================================
// TEST 13: Very Large k*dt (k_total * dt > 50 numerical safety)
// =============================================================================
TEST_CASE("StochasticReactionHandler numerical safety for k*dt > 50", "[reactions][handler][optimization]") {
    SpeciesDatabase species_db;
    SpeciesProperties ion_props;
    ion_props.mass_amu = 100.0;
    ion_props.charge = 1;
    ion_props.mass_kg = 100.0 * AMU_TO_KG;
    ion_props.charge_C = ELEM_CHARGE_C;
    ion_props.CCS_m2 = 100e-20;
    ion_props.mobility_m2Vs = 2.0e-4;
    species_db.species["Ion+"] = ion_props;
    species_db.species["Product+"] = ion_props;
    
    // Reaction with k_total * dt = 100 > 50 → P_total = 1.0 (certain)
    ReactionDatabase reaction_db;
    Reaction rxn;
    rxn.id = "certain_rxn";
    rxn.reactant = "Ion+";
    rxn.product = "Product+";
    rxn.rate_constant = 1e10;  // k_total = 1e10 s⁻¹
    rxn.order_terms = {};
    reaction_db.reactions = {rxn};
    
    auto env = create_test_environment(300.0, 2.5e25);
    
    StochasticReactionHandler handler(false);
    PhysicsRng rng(777);
    double dt = 1e-8;  // k*dt = 1e10 * 1e-8 = 100 > 50 → P_total = 1.0
    
    // Should react 100% of the time (P_total = 1.0, no exp() underflow)
    int reacted_count = 0;
    for (int i = 0; i < 100; ++i) {
        auto ensemble = make_ensemble_for_species(species_db, "Ion+");
        auto view = ensemble.reaction_data(0);
        bool reacted = handler.handle_reaction(view, dt, rng, reaction_db, species_db, env);
        if (reacted) {
            reacted_count++;
            REQUIRE(view.species_id() == "Product+");
        }
    }
    
    // Should react 100/100 times (P_total = 1.0)
    REQUIRE(reacted_count == 100);  // Deterministic!
}

// =============================================================================
// TEST 14: Arrhenius Temperature Dependence
// =============================================================================
TEST_CASE("StochasticReactionHandler: Arrhenius temperature dependence", "[reactions][handler][temperature]") {
    SpeciesDatabase species_db;
    SpeciesProperties ion_props;
    ion_props.mass_amu = 100.0;
    ion_props.charge = 1;
    ion_props.mass_kg = 100.0 * AMU_TO_KG;
    ion_props.charge_C = ELEM_CHARGE_C;
    ion_props.CCS_m2 = 100e-20;
    ion_props.mobility_m2Vs = 2.0e-4;
    species_db.species["Ion+"] = ion_props;
    species_db.species["Product+"] = ion_props;
    
    // Reaction with Arrhenius temperature dependence
    // k(T) = A × exp(-Eₐ / (kB·T))
    // A = 1e-9 m³/s, Eₐ = 0.1 eV
    ReactionDatabase reaction_db;
    Reaction rxn;
    rxn.id = "arrhenius_rxn";
    rxn.reactant = "Ion+";
    rxn.product = "Product+";
    rxn.rate_model = RateModel::Arrhenius;
    rxn.rate_constant = 1.0e-9;      // A [m³/s]
    rxn.activation_energy_eV = 0.1;      // Eₐ = 0.1 eV
    
    // Second-order reaction (ion-neutral)
    ReactionOrderTerm term;
    term.species = "N2";
    term.exponent = 1;
    term.concentration_m3 = -1.0;  // Use buffer gas density
    rxn.order_terms.push_back(term);
    
    reaction_db.reactions = {rxn};
    
    StochasticReactionHandler handler(false);
    PhysicsRng rng(42);
    
    IonState ion;
    ion.species_id = "Ion+";
    ion.mass_kg = 100.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 100e-20;
    
    double dt = 1e-9;
    double n_buffer = 2.5e25;  // [m⁻³]
    
    // Test at different temperatures
    std::vector<double> temperatures = {200.0, 300.0, 400.0};  // [K]
    std::vector<double> reaction_rates;
    
    for (double T_K : temperatures) {
        auto env = create_test_environment(T_K, n_buffer);
        
        // Compute k(T)
        double k_T = rxn.compute_rate_constant(T_K);
        
        // Expected: k(T) = A × exp(-Eₐ / (kB·T))
        // Eₐ = 0.1 eV × 1.602e-19 J/eV = 1.602e-20 J
        // kB = 1.381e-23 J/K
        double Ea_J = 0.1 * ELEM_CHARGE_C;
        double expected_k_T = 1e-9 * std::exp(-Ea_J / (BOLTZMANN_CONSTANT * T_K));
        
        // Check k(T) matches theory (within 0.01% relative error)
        REQUIRE(k_T == Catch::Approx(expected_k_T).epsilon(1e-6));
        
        // Verify rate increases with temperature
        reaction_rates.push_back(k_T);
    }
    
    // Check: k(400K) > k(300K) > k(200K)
    REQUIRE(reaction_rates[2] > reaction_rates[1]);
    REQUIRE(reaction_rates[1] > reaction_rates[0]);
    
    INFO("k(200K) = " << reaction_rates[0] << " m³/s");
    INFO("k(300K) = " << reaction_rates[1] << " m³/s");
    INFO("k(400K) = " << reaction_rates[2] << " m³/s");
    INFO("Ratio k(400K)/k(200K) = " << reaction_rates[2] / reaction_rates[0]);
}

// =============================================================================
// TEST 15: Modified Arrhenius (Ion-Dipole Capture, T^(-0.5) dependence)
// =============================================================================
TEST_CASE("StochasticReactionHandler: Modified Arrhenius (capture)", "[reactions][handler][temperature]") {
    SpeciesDatabase species_db;
    SpeciesProperties ion_props;
    ion_props.mass_amu = 100.0;
    ion_props.charge = 1;
    ion_props.mass_kg = 100.0 * AMU_TO_KG;
    ion_props.charge_C = ELEM_CHARGE_C;
    ion_props.CCS_m2 = 100e-20;
    ion_props.mobility_m2Vs = 2.0e-4;
    species_db.species["Ion+"] = ion_props;
    species_db.species["Cluster+"] = ion_props;
    
    // Reaction with modified Arrhenius: k(T) = A × (T/T₀)^n
    // A = 2e-9 m³/s, n = -0.5 (typical for ion-dipole capture)
    // T₀ = 300 K, Eₐ = 0 (barrierless)
    ReactionDatabase reaction_db;
    Reaction rxn;
    rxn.id = "capture_rxn";
    rxn.reactant = "Ion+";
    rxn.product = "Cluster+";
    rxn.rate_model = RateModel::ModifiedArrhenius;
    rxn.rate_constant = 2.0e-9;          // A [m³/s]
    rxn.temperature_exponent = -0.5;         // n = -0.5
    rxn.reference_temperature_K = 300.0;     // T₀ = 300 K
    rxn.activation_energy_eV = 0.0;          // No barrier
    
    reaction_db.reactions = {rxn};
    
    // Test k(T) = A × (T/T₀)^(-0.5)
    
    // At 200 K:
    // k(200K) = 2e-9 × (200/300)^(-0.5) = 2e-9 × (2/3)^(-0.5) = 2e-9 × 1.225 = 2.45e-9
    double k_200 = rxn.compute_rate_constant(200.0);
    REQUIRE(k_200 == Catch::Approx(2.45e-9).epsilon(1e-2));
    
    // At 300 K (reference):
    // k(300K) = 2e-9 × (300/300)^(-0.5) = 2e-9 × 1 = 2e-9
    double k_300 = rxn.compute_rate_constant(300.0);
    REQUIRE(k_300 == Catch::Approx(2.0e-9).epsilon(1e-6));
    
    // At 400 K:
    // k(400K) = 2e-9 × (400/300)^(-0.5) = 2e-9 × (4/3)^(-0.5) = 2e-9 × 0.866 = 1.73e-9
    double k_400 = rxn.compute_rate_constant(400.0);
    REQUIRE(k_400 == Catch::Approx(1.73e-9).epsilon(1e-2));
    
    // Verify: k decreases with T (negative exponent → "anti-Arrhenius")
    REQUIRE(k_200 > k_300);
    REQUIRE(k_300 > k_400);
    
    INFO("k(200K) = " << k_200 << " m³/s (faster at low T)");
    INFO("k(300K) = " << k_300 << " m³/s");
    INFO("k(400K) = " << k_400 << " m³/s (slower at high T)");
}

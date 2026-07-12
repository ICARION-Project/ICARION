// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file test_reaction_validation.cpp
 * @brief Tests for reaction order term validation rules
 * 
 * Validates reaction order term validation rules:
 * 1. exponent ∈ {0, 1, 2}
 * 2. concentration_m3 ≥ -1.0
 * 3. species exists in database (if not a runtime placeholder)
 * 4. No duplicate species
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "core/config/loader/ReactionLoader.h"
#include "core/config/loader/SpeciesLoader.h"
#include <json/json.h>
#include <fstream>
#include <iostream>
#include <sstream>

using namespace ICARION::config;

// Helper: Load species database for validation
SpeciesDatabase load_test_species_db() {
    Json::Value root;
    root["species"] = Json::objectValue;
    
    // Add test species
    Json::Value h3o_plus;
    h3o_plus["mass_amu"] = 19.0;
    h3o_plus["charge"] = 1;
    root["species"]["H3O+"] = h3o_plus;
    
    Json::Value h5o2_plus;
    h5o2_plus["mass_amu"] = 37.0;
    h5o2_plus["charge"] = 1;
    root["species"]["H5O2+"] = h5o2_plus;
    
    Json::Value nh4_plus;
    nh4_plus["mass_amu"] = 18.0;
    nh4_plus["charge"] = 1;
    root["species"]["NH4+"] = nh4_plus;
    
    Json::Value h2o;
    h2o["mass_amu"] = 18.0;
    h2o["charge"] = 0;
    root["species"]["H2O"] = h2o;
    
    Json::Value nh3;
    nh3["mass_amu"] = 17.0;
    nh3["charge"] = 0;
    root["species"]["NH3"] = nh3;
    
    return SpeciesLoader::load_from_json(root);
}

TEST_CASE("Reaction validation rules", "[reaction][validation]") {
    // Note: We test validation WITHOUT species_db to focus on order term rules
    // Species validation is tested separately in test_reaction_loader_unit.cpp
    
    SECTION("Valid: Pseudo-first-order with neutral=-1") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "valid_pseudo";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 3.5e-9;
        
        Json::Value order_term;
        order_term["species"] = "neutral";
        order_term["exponent"] = 1;
        order_term["concentration_m3"] = -1.0;
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);
        
        root["reactions"].append(rxn);
        
        REQUIRE_NOTHROW(ReactionLoader::load_from_json(root, nullptr));
    }

    SECTION("Valid: neutral placeholder with species database for non-equilibrium kinetics") {
        auto species_db = load_test_species_db();
        Json::Value root;
        root["reactions"] = Json::arrayValue;

        Json::Value rxn;
        rxn["id"] = "valid_neutral_placeholder";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 1.0e-15;

        Json::Value order_term;
        order_term["species"] = "neutral";
        order_term["exponent"] = 1;
        order_term["concentration_m3"] = -1.0;
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);

        root["reactions"].append(rxn);

        REQUIRE_NOTHROW(ReactionLoader::load_from_json(root, &species_db));
    }

    SECTION("Valid: M placeholder with species database") {
        auto species_db = load_test_species_db();
        Json::Value root;
        root["reactions"] = Json::arrayValue;

        Json::Value rxn;
        rxn["id"] = "valid_third_body_placeholder";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 1.0e-39;

        Json::Value order_term;
        order_term["species"] = "M";
        order_term["exponent"] = 1;
        order_term["concentration_m3"] = -1.0;
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);

        root["reactions"].append(rxn);

        REQUIRE_NOTHROW(ReactionLoader::load_from_json(root, &species_db));
    }

    SECTION("Invalid: equilibrium cannot use neutral placeholder") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;

        Json::Value rxn;
        rxn["id"] = "bad_equilibrium_neutral";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 1.0e-15;
        rxn["equilibrium"] = true;
        rxn["delta_r_H_J_mol"] = -42000.0;
        rxn["delta_r_S_J_molK"] = -95.0;

        Json::Value order_term;
        order_term["species"] = "neutral";
        order_term["exponent"] = 1;
        order_term["concentration_m3"] = -1.0;
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);

        root["reactions"].append(rxn);

        REQUIRE_THROWS_WITH(
            ReactionLoader::load_from_json(root, nullptr),
            Catch::Matchers::ContainsSubstring("equilibrium=true cannot use order term species 'neutral'")
        );
    }
    
    SECTION("Valid: Explicit concentration") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "valid_explicit";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "NH4+";
        rxn["rate_constant"] = 1.5e-9;
        
        Json::Value order_term;
        order_term["species"] = "NH3";
        order_term["exponent"] = 1;
        order_term["concentration_m3"] = 2.5e25;
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);
        
        root["reactions"].append(rxn);
        
        REQUIRE_NOTHROW(ReactionLoader::load_from_json(root, nullptr));
    }
    
    SECTION("Valid: Termolecular (exponent=2)") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "valid_termolecular";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 1.0e-28;
        
        Json::Value order_term;
        order_term["species"] = "H2O";
        order_term["exponent"] = 2;
        order_term["concentration_m3"] = 2.5e25;
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);
        
        root["reactions"].append(rxn);
        
        REQUIRE_NOTHROW(ReactionLoader::load_from_json(root, nullptr));
    }
    
    SECTION("❌ INVALID: Exponent > 2 (RULE #1)") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "bad_exponent";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 3.5e-9;
        
        Json::Value order_term;
        order_term["species"] = "H2O";
        order_term["exponent"] = 3;  // ❌ INVALID
        order_term["concentration_m3"] = 1e25;
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);
        
        root["reactions"].append(rxn);
        
        REQUIRE_THROWS_WITH(
            ReactionLoader::load_from_json(root, nullptr),
            Catch::Matchers::ContainsSubstring("exponent must be 0, 1, or 2")
        );
    }
    
    SECTION("❌ INVALID: Negative concentration < -1 (RULE #2)") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "bad_concentration";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 3.5e-9;
        
        Json::Value order_term;
        order_term["species"] = "H2O";
        order_term["exponent"] = 1;
        order_term["concentration_m3"] = -5.0;  // ❌ INVALID
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);
        
        root["reactions"].append(rxn);
        
        REQUIRE_THROWS_WITH(
            ReactionLoader::load_from_json(root, nullptr),
            Catch::Matchers::ContainsSubstring("must be ≥ -1.0")
        );
    }
    
    SECTION("❌ INVALID: Ambiguous negative concentration (RULE #2)") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "ambiguous_concentration";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 3.5e-9;
        
        Json::Value order_term;
        order_term["species"] = "H2O";
        order_term["exponent"] = 1;
        order_term["concentration_m3"] = -0.5;  // ❌ INVALID (ambiguous)
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);
        
        root["reactions"].append(rxn);
        
        REQUIRE_THROWS_WITH(
            ReactionLoader::load_from_json(root, nullptr),
            Catch::Matchers::ContainsSubstring("cannot be between -1.0 and 0.0")
        );
    }
    
    // NOTE: RULE #3 (species exists) is tested in test_reaction_loader_unit.cpp
    // because it requires a species_db. Here we test without species_db to focus
    // on order term validation rules #1, #2, #4.
    
    SECTION("❌ INVALID: Duplicate species (RULE #4)") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "duplicate_species";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 3.5e-9;
        
        Json::Value term1;
        term1["species"] = "H2O";
        term1["exponent"] = 1;
        term1["concentration_m3"] = 1e25;
        
        Json::Value term2;
        term2["species"] = "H2O";  // ❌ DUPLICATE
        term2["exponent"] = 1;
        term2["concentration_m3"] = 1e25;
        
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(term1);
        rxn["order"].append(term2);
        
        root["reactions"].append(rxn);
        
        REQUIRE_THROWS_WITH(
            ReactionLoader::load_from_json(root, nullptr),
            Catch::Matchers::ContainsSubstring("duplicate order term")
        );
    }
    
    SECTION("Valid: Multiple -1 fallbacks") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "two_fallbacks";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 3.5e-9;
        
        Json::Value term1;
        term1["species"] = "neutral";
        term1["exponent"] = 1;
        term1["concentration_m3"] = -1.0;
        
        Json::Value term2;
        term2["species"] = "neutral_gas";
        term2["exponent"] = 1;
        term2["concentration_m3"] = -1.0;
        
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(term1);
        rxn["order"].append(term2);
        
        root["reactions"].append(rxn);
        
        REQUIRE_NOTHROW(ReactionLoader::load_from_json(root, nullptr));
    }
}

TEST_CASE("Dimensional consistency warnings", "[reaction][warnings]") {
    SECTION("Typical SI bimolecular rate does not warn") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;

        Json::Value rxn;
        rxn["id"] = "typical_2nd_order";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 1.0e-15;

        Json::Value order_term;
        order_term["species"] = "H2O";
        order_term["exponent"] = 1;
        order_term["concentration_m3"] = 1e25;
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);

        root["reactions"].append(rxn);

        std::ostringstream captured;
        auto* old_buf = std::cout.rdbuf(captured.rdbuf());
        REQUIRE_NOTHROW(ReactionLoader::load_from_json(root, nullptr));
        std::cout.rdbuf(old_buf);
        CHECK_THAT(captured.str(), !Catch::Matchers::ContainsSubstring("outside typical"));
    }
    
    SECTION("Warning: Spontaneous decay with wrong units") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "suspicious_spontaneous";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 3.5e-9;  // Looks like [m³/s], but order=0 expects [s⁻¹]
        rxn["order"] = Json::arrayValue;  // Empty = order 0
        
        root["reactions"].append(rxn);
        
        // Should load but print warning
        REQUIRE_NOTHROW(ReactionLoader::load_from_json(root, nullptr));
    }
    
    SECTION("Warning: 2nd-order with suspicious units") {
        Json::Value root;
        root["reactions"] = Json::arrayValue;
        
        Json::Value rxn;
        rxn["id"] = "suspicious_2nd_order";
        rxn["reactant"] = "H3O+";
        rxn["product"] = "H5O2+";
        rxn["rate_constant"] = 1.0e-30;  // Looks like [m⁶/s], but order=1 expects [m³/s]
        
        Json::Value order_term;
        order_term["species"] = "H2O";
        order_term["exponent"] = 1;
        order_term["concentration_m3"] = 1e25;
        rxn["order"] = Json::arrayValue;
        rxn["order"].append(order_term);
        
        root["reactions"].append(rxn);
        
        // Should load but print warning
        REQUIRE_NOTHROW(ReactionLoader::load_from_json(root, nullptr));
    }
}

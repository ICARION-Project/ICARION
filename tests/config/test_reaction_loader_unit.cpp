// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/config/loader/ReactionLoader.h"
#include "core/config/loader/SpeciesLoader.h"
#include "core/config/types/ReactionConfig.h"
#include <filesystem>
#include <fstream>

using namespace ICARION::config;

// Helper to create temporary JSON file
struct TempReactionFile {
    std::filesystem::path path;
    
    TempReactionFile(const std::string& content) {
        path = std::filesystem::temp_directory_path() / "test_reactions.json";
        std::ofstream out(path);
        out << content;
    }
    
    ~TempReactionFile() {
        std::filesystem::remove(path);
    }
};

struct TempSpeciesFile {
    std::filesystem::path path;
    
    TempSpeciesFile(const std::string& content) {
        path = std::filesystem::temp_directory_path() / "test_species_rxn.json";
        std::ofstream out(path);
        out << content;
    }
    
    ~TempSpeciesFile() {
        std::filesystem::remove(path);
    }
};

TEST_CASE("ReactionLoader - Load simple reaction without order terms", "[reaction][loader]") {
    TempReactionFile tmp_file(R"({
        "reactions": [
            {
                "id": "rxn_test",
                "reactant": "A",
                "product": "B",
                "rate_constant_m3s": 1.5e-15
            }
        ]
    })");
    
    auto db = ReactionLoader::load(tmp_file.path);
    
    REQUIRE(db.size() == 1);
    const auto& rxn = db.reactions[0];
    
    CHECK(rxn.id == "rxn_test");
    CHECK(rxn.reactant == "A");
    CHECK(rxn.product == "B");
    CHECK_THAT(rxn.rate_constant_m3s, Catch::Matchers::WithinRel(1.5e-15, 1e-9));
    CHECK(rxn.order_terms.empty());
}

TEST_CASE("ReactionLoader - Load reaction with order terms", "[reaction][loader]") {
    TempReactionFile tmp_file(R"({
        "reactions": [
            {
                "id": "rxn_complex",
                "reactant": "H3O+",
                "product": "H5O2+",
                "rate_constant_m3s": 3.5e-15,
                "order": [
                    {
                        "species": "H2O",
                        "exponent": 1,
                        "concentration_m3": 2.5e25
                    }
                ]
            }
        ]
    })");
    
    auto db = ReactionLoader::load(tmp_file.path);
    
    REQUIRE(db.size() == 1);
    const auto& rxn = db.reactions[0];
    
    CHECK(rxn.order_terms.size() == 1);
    CHECK(rxn.order_terms[0].species == "H2O");
    CHECK(rxn.order_terms[0].exponent == 1);
    CHECK_THAT(rxn.order_terms[0].concentration_m3, Catch::Matchers::WithinRel(2.5e25, 1e-6));
}

TEST_CASE("ReactionLoader - Multiple reactions", "[reaction][loader]") {
    TempReactionFile tmp_file(R"({
        "reactions": [
            {"id": "rxn1", "reactant": "A", "product": "B", "rate_constant_m3s": 1e-15},
            {"id": "rxn2", "reactant": "B", "product": "C", "rate_constant_m3s": 2e-15},
            {"id": "rxn3", "reactant": "C", "product": "A", "rate_constant_m3s": 3e-15}
        ]
    })");
    
    auto db = ReactionLoader::load(tmp_file.path);
    
    REQUIRE(db.size() == 3);
    CHECK(db.reactions[0].id == "rxn1");
    CHECK(db.reactions[1].id == "rxn2");
    CHECK(db.reactions[2].id == "rxn3");
}

TEST_CASE("ReactionLoader - Species validation with database", "[reaction][loader][validation]") {
    // Create species database
    TempSpeciesFile species_file(R"({
        "species": {
            "H3O+": {"mass_amu": 19.02, "charge": 1},
            "H2O": {"mass_amu": 18.015, "charge": 0}
        }
    })");
    auto species_db = SpeciesLoader::load(species_file.path);
    
    // Valid reaction referencing existing species
    TempReactionFile rxn_file(R"({
        "reactions": [
            {
                "id": "valid_rxn",
                "reactant": "H3O+",
                "product": "H2O",
                "rate_constant_m3s": 1e-15,
                "order": [{"species": "H2O", "exponent": 1}]
            }
        ]
    })");
    
    // Should not throw
    CHECK_NOTHROW(ReactionLoader::load(rxn_file.path, &species_db));
}

TEST_CASE("ReactionLoader - Invalid species reference throws", "[reaction][loader][validation][error]") {
    TempSpeciesFile species_file(R"({
        "species": {
            "H3O+": {"mass_amu": 19.02, "charge": 1}
        }
    })");
    auto species_db = SpeciesLoader::load(species_file.path);
    
    TempReactionFile rxn_file(R"({
        "reactions": [
            {
                "id": "invalid_rxn",
                "reactant": "UNKNOWN",
                "product": "H3O+",
                "rate_constant_m3s": 1e-15
            }
        ]
    })");
    
    CHECK_THROWS_AS(ReactionLoader::load(rxn_file.path, &species_db), std::runtime_error);
}

TEST_CASE("ReactionLoader - Missing required fields throws", "[reaction][loader][error]") {
    TempReactionFile tmp_file(R"({
        "reactions": [
            {
                "id": "incomplete",
                "reactant": "A"
            }
        ]
    })");
    
    CHECK_THROWS_AS(ReactionLoader::load(tmp_file.path), std::runtime_error);
}

TEST_CASE("ReactionLoader - File not found throws", "[reaction][loader][error]") {
    CHECK_THROWS_AS(
        ReactionLoader::load("/nonexistent/reactions.json"),
        std::runtime_error
    );
}

TEST_CASE("Reaction - Effective rate calculation", "[reaction][calculation]") {
    Reaction rxn;
    rxn.id = "test";
    rxn.rate_constant_m3s = 1e-15;
    
    // No order terms -> k_eff = k (empty concentration map)
    std::unordered_map<std::string, double> empty_conc;
    CHECK_THAT(rxn.effective_rate_s(empty_conc), Catch::Matchers::WithinRel(1e-15, 1e-9));
    
    // Add order term: [H2O]^1 = 2e25 m⁻³
    ReactionOrderTerm term;
    term.species = "H2O";
    term.exponent = 1;
    term.concentration_m3 = 2e25;
    rxn.order_terms.push_back(term);
    
    // k_eff = k * [H2O]^1 = 1e-15 * 2e25 = 2e10
    std::unordered_map<std::string, double> conc{{"H2O", 2e25}};
    CHECK_THAT(rxn.effective_rate_s(conc), Catch::Matchers::WithinRel(2e10, 1e-3));
}

TEST_CASE("Reaction - Multiple order terms", "[reaction][calculation]") {
    Reaction rxn;
    rxn.rate_constant_m3s = 1e-15;
    
    // [A]^2 = 1e24 m⁻³, [B]^1 = 3e25 m⁻³
    ReactionOrderTerm term_a;
    term_a.species = "A";
    term_a.exponent = 2;
    term_a.concentration_m3 = 1e24;
    
    ReactionOrderTerm term_b;
    term_b.species = "B";
    term_b.exponent = 1;
    term_b.concentration_m3 = 3e25;
    
    rxn.order_terms.push_back(term_a);
    rxn.order_terms.push_back(term_b);
    
    // k_eff = k * (1e24)^2 * (3e25)^1 = 1e-15 * 1e48 * 3e25 = 3e58
    std::unordered_map<std::string, double> conc{{"A", 1e24}, {"B", 3e25}};
    CHECK_THAT(rxn.effective_rate_s(conc), Catch::Matchers::WithinRel(3e58, 1e-3));
}

TEST_CASE("ReactionDatabase - Get reactions for species", "[reaction][database]") {
    ReactionDatabase db;
    
    Reaction rxn1;
    rxn1.id = "r1";
    rxn1.reactant = "A";
    rxn1.product = "B";
    rxn1.rate_constant_m3s = 1e-15;
    
    Reaction rxn2;
    rxn2.id = "r2";
    rxn2.reactant = "A";
    rxn2.product = "C";
    rxn2.rate_constant_m3s = 2e-15;
    
    Reaction rxn3;
    rxn3.id = "r3";
    rxn3.reactant = "B";
    rxn3.product = "C";
    rxn3.rate_constant_m3s = 3e-15;
    
    db.reactions.push_back(rxn1);
    db.reactions.push_back(rxn2);
    db.reactions.push_back(rxn3);
    
    auto reactions_A = db.get_reactions_for("A");
    REQUIRE(reactions_A.size() == 2);
    CHECK(reactions_A[0]->id == "r1");
    CHECK(reactions_A[1]->id == "r2");
    
    auto reactions_B = db.get_reactions_for("B");
    REQUIRE(reactions_B.size() == 1);
    CHECK(reactions_B[0]->id == "r3");
    
    auto reactions_C = db.get_reactions_for("C");
    CHECK(reactions_C.empty());
}

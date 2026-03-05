// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/config/loader/SpeciesLoader.h"
#include "core/config/types/SpeciesConfig.h"
#include <filesystem>
#include <fstream>

using namespace ICARION::config;

// Helper to create temporary JSON file
struct TempSpeciesFile {
    std::filesystem::path path;
    
    TempSpeciesFile(const std::string& content) {
        path = std::filesystem::temp_directory_path() / "test_species.json";
        std::ofstream out(path);
        out << content;
    }
    
    ~TempSpeciesFile() {
        std::filesystem::remove(path);
    }
};

TEST_CASE("SpeciesLoader - Load valid ion species", "[species][loader]") {
    TempSpeciesFile tmp_file(R"({
        "species": {
            "H3O+": {
                "name": "Hydronium",
                "mass_amu": 19.02,
                "charge": 1,
                "mobility_cm2Vs": 2.8,
                "CCS_A2": 11.0
            }
        }
    })");
    
    auto db = SpeciesLoader::load(tmp_file.path);
    
    REQUIRE(db.size() == 1);
    REQUIRE(db.has("H3O+"));
    
    const auto& species = db.get("H3O+");
    CHECK(species.id == "H3O+");
    CHECK(species.name == "Hydronium");
    CHECK_THAT(species.mass_amu, Catch::Matchers::WithinRel(19.02, 1e-6));
    CHECK(species.charge == 1);
    CHECK(species.mobility_cm2Vs.has_value());
    CHECK_THAT(*species.mobility_cm2Vs, Catch::Matchers::WithinRel(2.8, 1e-6));
    CHECK(species.CCS_A2.has_value());
    CHECK_THAT(*species.CCS_A2, Catch::Matchers::WithinRel(11.0, 1e-6));
    
    // Check SI conversion
    CHECK_THAT(species.mass_kg, Catch::Matchers::WithinRel(19.02 * 1.66053906660e-27, 1e-6));
    CHECK_THAT(species.charge_C, Catch::Matchers::WithinRel(1.602176634e-19, 1e-6));
    CHECK_THAT(species.mobility_m2Vs, Catch::Matchers::WithinRel(2.8e-4, 1e-6));
}

TEST_CASE("SpeciesLoader - Load neutral species", "[species][loader]") {
    TempSpeciesFile tmp_file(R"({
        "species": {
            "H2O": {
                "name": "Water",
                "mass_amu": 18.015,
                "charge": 0,
                "polarizability_A3": 1.47
            }
        }
    })");
    
    auto db = SpeciesLoader::load(tmp_file.path);
    
    REQUIRE(db.size() == 1);
    const auto& species = db.get("H2O");
    
    CHECK(species.charge == 0);
    CHECK_FALSE(species.mobility_cm2Vs.has_value());
    CHECK_FALSE(species.CCS_A2.has_value());
    CHECK(species.polarizability_A3.has_value());
    CHECK_THAT(*species.polarizability_A3, Catch::Matchers::WithinRel(1.47, 1e-6));
}

TEST_CASE("SpeciesLoader - Load multiple species", "[species][loader]") {
    TempSpeciesFile tmp_file(R"({
        "species": {
            "H3O+": {"mass_amu": 19.02, "charge": 1},
            "H2O": {"mass_amu": 18.015, "charge": 0},
            "O2-": {"mass_amu": 31.998, "charge": -1}
        }
    })");
    
    auto db = SpeciesLoader::load(tmp_file.path);
    
    REQUIRE(db.size() == 3);
    CHECK(db.has("H3O+"));
    CHECK(db.has("H2O"));
    CHECK(db.has("O2-"));
}

TEST_CASE("SpeciesLoader - Missing required field throws", "[species][loader][error]") {
    TempSpeciesFile tmp_file(R"({
        "species": {
            "H3O+": {
                "mass_amu": 19.02
            }
        }
    })");
    
    CHECK_THROWS_AS(SpeciesLoader::load(tmp_file.path), std::runtime_error);
}

TEST_CASE("SpeciesLoader - Invalid mass throws", "[species][loader][error]") {
    TempSpeciesFile tmp_file(R"({
        "species": {
            "H3O+": {
                "mass_amu": -1.0,
                "charge": 1
            }
        }
    })");
    
    CHECK_THROWS_AS(SpeciesLoader::load(tmp_file.path), std::runtime_error);
}

TEST_CASE("SpeciesLoader - File not found throws", "[species][loader][error]") {
    CHECK_THROWS_AS(
        SpeciesLoader::load("/nonexistent/path.json"),
        std::runtime_error
    );
}

TEST_CASE("SpeciesLoader - Invalid JSON throws", "[species][loader][error]") {
    TempSpeciesFile tmp_file("{ invalid json }");
    CHECK_THROWS_AS(SpeciesLoader::load(tmp_file.path), std::runtime_error);
}

TEST_CASE("SpeciesConfig - Unit conversions are correct", "[species][conversion]") {
    SpeciesProperties props;
    props.id = "test";
    props.mass_amu = 20.0;
    props.charge = 2;
    props.mobility_cm2Vs = 3.0;
    props.CCS_A2 = 15.0;
    props.polarizability_A3 = 2.0;
    
    props.convert_to_SI();
    
    // Mass: amu -> kg
    CHECK_THAT(props.mass_kg, Catch::Matchers::WithinRel(20.0 * 1.66053906660e-27, 1e-9));
    
    // Charge: e -> C
    CHECK_THAT(props.charge_C, Catch::Matchers::WithinRel(2.0 * 1.602176634e-19, 1e-9));
    
    // Mobility: cm²/Vs -> m²/Vs
    CHECK_THAT(props.mobility_m2Vs, Catch::Matchers::WithinRel(3.0e-4, 1e-9));
    
    // CCS: Ų -> m²
    CHECK_THAT(props.CCS_m2, Catch::Matchers::WithinRel(15.0e-20, 1e-9));
    
    // Polarizability: ų -> m³
    CHECK_THAT(props.polarizability_m3, Catch::Matchers::WithinRel(2.0e-30, 1e-9));
}

TEST_CASE("SpeciesDatabase - Dictionary lookup", "[species][database]") {
    SpeciesDatabase db;
    
    SpeciesProperties props1;
    props1.id = "A";
    props1.mass_amu = 10.0;
    props1.charge = 1;
    props1.convert_to_SI();
    
    SpeciesProperties props2;
    props2.id = "B";
    props2.mass_amu = 20.0;
    props2.charge = -1;
    props2.convert_to_SI();
    
    db.species["A"] = props1;
    db.species["B"] = props2;
    
    CHECK(db.size() == 2);
    CHECK(db.has("A"));
    CHECK(db.has("B"));
    CHECK_FALSE(db.has("C"));
    
    CHECK(db.get("A").charge == 1);
    CHECK(db.get("B").charge == -1);
}

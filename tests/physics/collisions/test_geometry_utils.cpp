// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file test_geometry_utils.cpp
 * @brief Unit tests for geometry loading utilities
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/physics/collisions/geometryUtils.h"
#include "core/io/moleculeLoader.h"
#include <filesystem>
#include <fstream>
#include <cmath>

using namespace ICARION::physics;
using namespace ICARION::io;
namespace fs = std::filesystem;

// ============================================================================
// Test Fixtures and Helpers
// ============================================================================

namespace {
    /**
     * Create temporary test directory with molecule JSON files
     */
    class TempMoleculeDir {
    public:
        TempMoleculeDir() {
            // Create unique temp directory
            temp_dir = fs::temp_directory_path() / ("icarion_test_" + std::to_string(rand()));
            fs::create_directories(temp_dir);
        }
        
        ~TempMoleculeDir() {
            // Cleanup
            if (fs::exists(temp_dir)) {
                fs::remove_all(temp_dir);
            }
        }
        
        /**
         * Create H3O+ molecule file (4 atoms)
         * Positions in Angström, converted to meters by loader
         */
        void create_h3o_file(const std::string& filename = "H3O+.json") {
            fs::path filepath = temp_dir / filename;
            std::ofstream file(filepath);
            file << R"({
    "molecule": {
        "name": "H3O+",
        "diameter_m": 3.3e-10,
        "CCS_m2": 24.9e-19,
        "atoms": [
            {
                "element": "O",
                "pos": [0.0, 0.0, 0.0],
                "mass_u": 15.999,
                "partial_charge_e": 0.004,
                "LJ_sigma_angstrom": 2.434,
                "LJ_epsilon_eV": 0.103
            },
            {
                "element": "H",
                "pos": [0.96, 0.0, 0.0],
                "mass_u": 1.008,
                "partial_charge_e": 0.332,
                "LJ_sigma_angstrom": 2.261,
                "LJ_epsilon_eV": 0.060
            },
            {
                "element": "H",
                "pos": [-0.24, 0.93, 0.0],
                "mass_u": 1.008,
                "partial_charge_e": 0.332,
                "LJ_sigma_angstrom": 2.261,
                "LJ_epsilon_eV": 0.060
            },
            {
                "element": "H",
                "pos": [-0.24, -0.93, 0.0],
                "mass_u": 1.008,
                "partial_charge_e": 0.332,
                "LJ_sigma_angstrom": 2.261,
                "LJ_epsilon_eV": 0.060
            }
        ]
    }
})";
            file.close();
        }
        
        /**
         * Create N2 molecule file (2 atoms)
         * Positions in Angström, converted to meters by loader
         */
        void create_n2_file(const std::string& filename = "N2.json") {
            fs::path filepath = temp_dir / filename;
            std::ofstream file(filepath);
            file << R"({
    "molecule": {
        "name": "N2",
        "diameter_m": 3.8e-10,
        "CCS_m2": 27.5e-19,
        "atoms": [
            {
                "element": "N",
                "pos": [0.0, 0.0, 0.55],
                "mass_u": 14.007,
                "partial_charge_e": 0.0,
                "LJ_sigma_angstrom": 3.667,
                "LJ_epsilon_eV": 0.0086
            },
            {
                "element": "N",
                "pos": [0.0, 0.0, -0.55],
                "mass_u": 14.007,
                "partial_charge_e": 0.0,
                "LJ_sigma_angstrom": 3.667,
                "LJ_epsilon_eV": 0.0086
            }
        ]
    }
})";
            file.close();
        }
        
        /**
         * Create single-atom test molecule
         * Position in Angström, converted to meters by loader
         */
        void create_single_atom_file(const std::string& filename = "Ar.json") {
            fs::path filepath = temp_dir / filename;
            std::ofstream file(filepath);
            file << R"({
    "molecule": {
        "name": "Ar",
        "diameter_m": 3.4e-10,
        "CCS_m2": 25.0e-19,
        "atoms": [
            {
                "element": "Ar",
                "pos": [0.0, 0.0, 0.0],
                "mass_u": 39.948,
                "partial_charge_e": 0.0,
                "LJ_sigma_angstrom": 3.405,
                "LJ_epsilon_eV": 0.0103
            }
        ]
    }
})";
            file.close();
        }
        
        /**
         * Create invalid JSON file
         */
        void create_invalid_file(const std::string& filename = "invalid.json") {
            fs::path filepath = temp_dir / filename;
            std::ofstream file(filepath);
            file << "{ this is not valid JSON }";
            file.close();
        }
        
        fs::path get_path() const { return temp_dir; }
        std::string get_path_string() const { return temp_dir.string(); }
        
    private:
        fs::path temp_dir;
    };
} // anonymous namespace

// ============================================================================
// Test: convert_molecule_to_geometry()
// ============================================================================

TEST_CASE("convert_molecule_to_geometry - Basic conversion", "[geometry_utils][unit]") {
    TempMoleculeDir temp;
    temp.create_h3o_file();
    
    auto molecule = load_molecule((temp.get_path() / "H3O+.json").string());
    auto geom = convert_molecule_to_geometry(molecule);
    
    SECTION("Correct number of atoms") {
        REQUIRE(geom.first.size() == 4);  // centers
        REQUIRE(geom.second.size() == 4); // radii
    }
    
    SECTION("Position conversion (Angstrom to meters)") {
        // O atom at origin
        REQUIRE_THAT(geom.first[0].x, Catch::Matchers::WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(geom.first[0].y, Catch::Matchers::WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(geom.first[0].z, Catch::Matchers::WithinAbs(0.0, 1e-12));
        
        // H atom at [0.96, 0, 0] Å = [0.96e-10, 0, 0] m
        REQUIRE_THAT(geom.first[1].x, Catch::Matchers::WithinAbs(0.96e-10, 1e-12));
        REQUIRE_THAT(geom.first[1].y, Catch::Matchers::WithinAbs(0.0, 1e-12));
        REQUIRE_THAT(geom.first[1].z, Catch::Matchers::WithinAbs(0.0, 1e-12));
    }
    
    SECTION("Radius conversion (LJ sigma to hard-sphere)") {
        // O: σ = 2.434e-10 m, radius = σ/2 = 1.217e-10 m
        REQUIRE_THAT(geom.second[0], Catch::Matchers::WithinAbs(1.217e-10, 1e-13));
        
        // H: σ = 2.261e-10 m, radius = σ/2 = 1.1305e-10 m
        REQUIRE_THAT(geom.second[1], Catch::Matchers::WithinAbs(1.1305e-10, 1e-13));
        REQUIRE_THAT(geom.second[2], Catch::Matchers::WithinAbs(1.1305e-10, 1e-13));
        REQUIRE_THAT(geom.second[3], Catch::Matchers::WithinAbs(1.1305e-10, 1e-13));
    }
}

TEST_CASE("convert_molecule_to_geometry - Single atom", "[geometry_utils][unit]") {
    TempMoleculeDir temp;
    temp.create_single_atom_file();
    
    auto molecule = load_molecule((temp.get_path() / "Ar.json").string());
    auto geom = convert_molecule_to_geometry(molecule);
    
    REQUIRE(geom.first.size() == 1);
    REQUIRE(geom.second.size() == 1);
    
    // Ar: σ = 3.405e-10 m, radius = 1.7025e-10 m
    REQUIRE_THAT(geom.second[0], Catch::Matchers::WithinAbs(1.7025e-10, 1e-13));
}

TEST_CASE("convert_molecule_to_geometry - Empty molecule", "[geometry_utils][unit]") {
    Molecule empty_mol;
    empty_mol.name = "Empty";
    empty_mol.diameter_m = 0.0;
    empty_mol.CCS_m2 = 0.0;
    empty_mol.atoms = {};
    empty_mol.total_mass_u = 0.0;
    empty_mol.total_charge_e = 0;
    empty_mol.center_of_mass_m = Vec3{0, 0, 0};
    
    auto geom = convert_molecule_to_geometry(empty_mol);
    
    REQUIRE(geom.first.empty());
    REQUIRE(geom.second.empty());
}

// ============================================================================
// Test: load_geometry_map() - Single file mode
// ============================================================================

TEST_CASE("load_geometry_map - Single file mode", "[geometry_utils][unit]") {
    TempMoleculeDir temp;
    temp.create_h3o_file("molecules.json");
    
    std::unordered_set<std::string> species = {"H3O+"};
    auto geom_map = load_geometry_map(species, (temp.get_path() / "molecules.json").string());
    
    SECTION("Species loaded successfully") {
        REQUIRE(geom_map.count("H3O+") == 1);
        REQUIRE(geom_map["H3O+"].first.size() == 4);
        REQUIRE(geom_map["H3O+"].second.size() == 4);
    }
}

// ============================================================================
// Test: load_geometry_map() - Directory mode
// ============================================================================

TEST_CASE("load_geometry_map - Directory mode", "[geometry_utils][unit]") {
    TempMoleculeDir temp;
    temp.create_h3o_file();
    temp.create_n2_file();
    temp.create_single_atom_file();
    
    SECTION("Load multiple species from directory") {
        std::unordered_set<std::string> species = {"H3O+", "N2", "Ar"};
        auto geom_map = load_geometry_map(species, temp.get_path_string());
        
        REQUIRE(geom_map.size() == 3);
        
        // H3O+ (4 atoms)
        REQUIRE(geom_map.count("H3O+") == 1);
        REQUIRE(geom_map["H3O+"].first.size() == 4);
        
        // N2 (2 atoms)
        REQUIRE(geom_map.count("N2") == 1);
        REQUIRE(geom_map["N2"].first.size() == 2);
        
        // Ar (1 atom)
        REQUIRE(geom_map.count("Ar") == 1);
        REQUIRE(geom_map["Ar"].first.size() == 1);
    }
    
    SECTION("Load subset of available species") {
        std::unordered_set<std::string> species = {"H3O+"};
        auto geom_map = load_geometry_map(species, temp.get_path_string());
        
        REQUIRE(geom_map.size() == 1);
        REQUIRE(geom_map.count("H3O+") == 1);
        REQUIRE(geom_map["H3O+"].first.size() == 4);
    }
}

// ============================================================================
// Test: load_geometry_map() - Missing species fallback
// ============================================================================

TEST_CASE("load_geometry_map - Missing species (default: silent fallback)", "[geometry_utils][unit]") {
    TempMoleculeDir temp;
    temp.create_h3o_file();
    
    SECTION("Missing file returns empty geometry") {
        std::unordered_set<std::string> species = {"H3O+", "NonExistent"};
        auto geom_map = load_geometry_map(species, temp.get_path_string());
        
        // H3O+ loaded successfully
        REQUIRE(geom_map.count("H3O+") == 1);
        REQUIRE(geom_map["H3O+"].first.size() == 4);
        
        // NonExistent has empty geometry (fallback signal)
        REQUIRE(geom_map.count("NonExistent") == 1);
        REQUIRE(geom_map["NonExistent"].first.empty());
        REQUIRE(geom_map["NonExistent"].second.empty());
    }
    
    SECTION("Invalid JSON returns empty geometry") {
        temp.create_invalid_file("Bad.json");
        std::unordered_set<std::string> species = {"Bad"};
        auto geom_map = load_geometry_map(species, temp.get_path_string());
        
        REQUIRE(geom_map.count("Bad") == 1);
        REQUIRE(geom_map["Bad"].first.empty());
    }
}

TEST_CASE("load_geometry_map - Missing species (strict_mode=true)", "[geometry_utils][unit]") {
    TempMoleculeDir temp;
    temp.create_h3o_file();
    
    SECTION("Missing file throws in strict mode") {
        std::unordered_set<std::string> species = {"NonExistent"};
        
        REQUIRE_THROWS_AS(
            load_geometry_map(species, temp.get_path_string(), true),
            std::runtime_error
        );
    }
    
    SECTION("Invalid JSON throws in strict mode") {
        temp.create_invalid_file("Bad.json");
        std::unordered_set<std::string> species = {"Bad"};
        
        REQUIRE_THROWS_AS(
            load_geometry_map(species, temp.get_path_string(), true),
            std::runtime_error
        );
    }
    
    SECTION("Partial success: first species ok, second fails") {
        std::unordered_set<std::string> species = {"H3O+", "NonExistent"};
        
        // Should throw on NonExistent even though H3O+ loads successfully
        REQUIRE_THROWS_AS(
            load_geometry_map(species, temp.get_path_string(), true),
            std::runtime_error
        );
    }
}

// ============================================================================
// Test: Edge cases
// ============================================================================

TEST_CASE("load_geometry_map - Empty species set", "[geometry_utils][unit]") {
    TempMoleculeDir temp;
    temp.create_h3o_file();
    
    std::unordered_set<std::string> empty_species;
    auto geom_map = load_geometry_map(empty_species, temp.get_path_string());
    
    REQUIRE(geom_map.empty());
}

TEST_CASE("load_geometry_map - Nonexistent directory", "[geometry_utils][unit]") {
    std::unordered_set<std::string> species = {"H3O+"};
    
    SECTION("Silent fallback (default)") {
        auto geom_map = load_geometry_map(species, "/nonexistent/path/");
        
        // Returns empty geometry (fallback)
        REQUIRE(geom_map.count("H3O+") == 1);
        REQUIRE(geom_map["H3O+"].first.empty());
    }
    
    SECTION("Strict mode throws") {
        REQUIRE_THROWS_AS(
            load_geometry_map(species, "/nonexistent/path/", true),
            std::runtime_error
        );
    }
}

TEST_CASE("load_geometry_map - Multiple species, mixed success", "[geometry_utils][unit]") {
    TempMoleculeDir temp;
    temp.create_h3o_file();
    temp.create_n2_file();
    
    std::unordered_set<std::string> species = {"H3O+", "N2", "Missing1", "Missing2"};
    
    SECTION("Silent fallback mode") {
        auto geom_map = load_geometry_map(species, temp.get_path_string());
        
        REQUIRE(geom_map.size() == 4);
        
        // Successful loads
        REQUIRE(geom_map["H3O+"].first.size() == 4);
        REQUIRE(geom_map["N2"].first.size() == 2);
        
        // Failed loads (empty geometry)
        REQUIRE(geom_map["Missing1"].first.empty());
        REQUIRE(geom_map["Missing2"].first.empty());
    }
}

// ============================================================================
// Test: Physical correctness
// ============================================================================

TEST_CASE("convert_molecule_to_geometry - Physical constants verification", "[geometry_utils][unit]") {
    TempMoleculeDir temp;
    temp.create_n2_file();
    
    auto molecule = load_molecule((temp.get_path() / "N2.json").string());
    auto geom = convert_molecule_to_geometry(molecule);
    
    SECTION("LJ sigma to hard-sphere radius factor is 0.5") {
        // N2: σ = 3.667e-10 m
        // Expected radius = σ/2 = 1.8335e-10 m
        double expected_radius = 3.667e-10 / 2.0;
        
        REQUIRE_THAT(geom.second[0], Catch::Matchers::WithinAbs(expected_radius, 1e-13));
        REQUIRE_THAT(geom.second[1], Catch::Matchers::WithinAbs(expected_radius, 1e-13));
    }
    
    SECTION("Position scaling: Angstrom to meters") {
        // N atoms at ±0.55 Å = ±5.5e-11 m
        double expected_z = 0.55e-10;
        
        REQUIRE_THAT(geom.first[0].z, Catch::Matchers::WithinAbs(expected_z, 1e-13));
        REQUIRE_THAT(geom.first[1].z, Catch::Matchers::WithinAbs(-expected_z, 1e-13));
    }
}

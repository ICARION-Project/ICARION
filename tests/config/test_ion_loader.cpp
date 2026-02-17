// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/config/loader/IonLoader.h"
#include "core/config/types/SpeciesConfig.h"
#include "utils/constants.h"
#include <random>
#include <cmath>

using namespace ICARION::config;
using namespace ICARION::core;

TEST_CASE("IonLoader: Point position distribution", "[ion][loader]") {
    // Setup species database
    SpeciesDatabase species_db;
    SpeciesProperties h3o_plus;
    h3o_plus.mass_amu = 19.02;
    h3o_plus.mass_kg = 19.02 * 1.66054e-27;
    h3o_plus.charge = 1;
    h3o_plus.charge_C = 1.602176634e-19;
    h3o_plus.mobility_m2Vs = 2.8e-4;
    h3o_plus.CCS_m2 = 1.1e-19;
    species_db.species["H3O+"] = h3o_plus;
    
    // Setup ion config
    IonConfig config;
    IonSpeciesConfig spec;
    spec.species_id = "H3O+";
    spec.count = 10;
    spec.position.type = PositionDistribution::Point;
    spec.position.center = {0.001, 0.002, 0.003};
    spec.velocity.type = VelocityDistribution::Fixed;
    spec.velocity.value = {10, 20, 30};
    config.species.push_back(spec);
    
    // Generate ions
    std::mt19937 rng(42);
    auto result = IonLoader::generate_ions(config, species_db, rng);
    
    REQUIRE(result.success());
    REQUIRE(result.ions.size() == 10);
    
    // Check all ions at same point
    for (const auto& ion : result.ions) {
        REQUIRE(ion.species_id == "H3O+");
        REQUIRE_THAT(ion.pos.x, Catch::Matchers::WithinAbs(0.001, 1e-10));
        REQUIRE_THAT(ion.pos.y, Catch::Matchers::WithinAbs(0.002, 1e-10));
        REQUIRE_THAT(ion.pos.z, Catch::Matchers::WithinAbs(0.003, 1e-10));
        REQUIRE_THAT(ion.vel.x, Catch::Matchers::WithinAbs(10.0, 1e-10));
        REQUIRE_THAT(ion.vel.y, Catch::Matchers::WithinAbs(20.0, 1e-10));
        REQUIRE_THAT(ion.vel.z, Catch::Matchers::WithinAbs(30.0, 1e-10));
    }
}

TEST_CASE("IonLoader: Gaussian position distribution", "[ion][loader]") {
    SpeciesDatabase species_db;
    SpeciesProperties h3o_plus;
    h3o_plus.mass_amu = 19.02;
    h3o_plus.mass_kg = 19.02 * 1.66054e-27;
    h3o_plus.charge = 1;
    h3o_plus.charge_C = 1.602176634e-19;
    h3o_plus.mobility_m2Vs = 2.8e-4;
    h3o_plus.CCS_m2 = 1.1e-19;
    species_db.species["H3O+"] = h3o_plus;
    
    IonConfig config;
    IonSpeciesConfig spec;
    spec.species_id = "H3O+";
    spec.count = 1000;
    spec.position.type = PositionDistribution::Gaussian;
    spec.position.center = {0, 0, 0};
    spec.position.std_dev = {0.001, 0.002, 0.003};
    spec.velocity.type = VelocityDistribution::Fixed;
    spec.velocity.value = {0, 0, 0};
    config.species.push_back(spec);
    
    std::mt19937 rng(42);
    auto result = IonLoader::generate_ions(config, species_db, rng);
    
    REQUIRE(result.success());
    REQUIRE(result.ions.size() == 1000);
    
    // Check statistical properties
    double mean_x = 0, mean_y = 0, mean_z = 0;
    for (const auto& ion : result.ions) {
        mean_x += ion.pos.x;
        mean_y += ion.pos.y;
        mean_z += ion.pos.z;
    }
    mean_x /= result.ions.size();
    mean_y /= result.ions.size();
    mean_z /= result.ions.size();
    
    // Mean should be near center (within 3σ/√N)
    REQUIRE(std::abs(mean_x) < 0.0001);
    REQUIRE(std::abs(mean_y) < 0.0002);
    REQUIRE(std::abs(mean_z) < 0.0003);
    
    // Check std dev
    double var_x = 0, var_y = 0, var_z = 0;
    for (const auto& ion : result.ions) {
        var_x += (ion.pos.x - mean_x) * (ion.pos.x - mean_x);
        var_y += (ion.pos.y - mean_y) * (ion.pos.y - mean_y);
        var_z += (ion.pos.z - mean_z) * (ion.pos.z - mean_z);
    }
    var_x /= result.ions.size();
    var_y /= result.ions.size();
    var_z /= result.ions.size();
    
    double std_x = std::sqrt(var_x);
    double std_y = std::sqrt(var_y);
    double std_z = std::sqrt(var_z);
    
    // Std should be close to configured (within 10%)
    REQUIRE_THAT(std_x, Catch::Matchers::WithinRel(0.001, 0.1));
    REQUIRE_THAT(std_y, Catch::Matchers::WithinRel(0.002, 0.1));
    REQUIRE_THAT(std_z, Catch::Matchers::WithinRel(0.003, 0.1));
}

TEST_CASE("IonLoader: Thermal velocity always random", "[ion][loader]") {
    SpeciesDatabase species_db;
    SpeciesProperties h3o_plus;
    h3o_plus.mass_amu = 19.02;
    h3o_plus.mass_kg = 19.02 * 1.66054e-27;
    h3o_plus.charge = 1;
    h3o_plus.charge_C = 1.602176634e-19;
    h3o_plus.mobility_m2Vs = 2.8e-4;
    h3o_plus.CCS_m2 = 1.1e-19;
    species_db.species["H3O+"] = h3o_plus;
    
    IonConfig config;
    IonSpeciesConfig spec;
    spec.species_id = "H3O+";
    spec.count = 1000;
    spec.position.type = PositionDistribution::Point;
    spec.position.center = {0, 0, 0};
    spec.velocity.type = VelocityDistribution::Thermal;
    spec.velocity.temperature_K = 300;
    config.species.push_back(spec);
    
    std::mt19937 rng(42);
    auto result = IonLoader::generate_ions(config, species_db, rng);
    
    REQUIRE(result.success());
    REQUIRE(result.ions.size() == 1000);
    
    // Check that mean velocity is near zero (random directions, no drift)
    double mean_vx = 0, mean_vy = 0, mean_vz = 0;
    for (const auto& ion : result.ions) {
        mean_vx += ion.vel.x;
        mean_vy += ion.vel.y;
        mean_vz += ion.vel.z;
    }
    mean_vx /= result.ions.size();
    mean_vy /= result.ions.size();
    mean_vz /= result.ions.size();
    
    // Calculate expected thermal velocity std dev
    double sigma = std::sqrt(BOLTZMANN_CONSTANT * 300 / h3o_plus.mass_kg);
    double expected_mean_uncertainty = sigma / std::sqrt(result.ions.size());  // ~3.16
    
    // Mean velocity should be near zero (within 5σ uncertainty)
    REQUIRE(std::abs(mean_vx) < 5 * expected_mean_uncertainty);
    REQUIRE(std::abs(mean_vy) < 5 * expected_mean_uncertainty);
    REQUIRE(std::abs(mean_vz) < 5 * expected_mean_uncertainty);
}

TEST_CASE("IonLoader: Multiple species with different boundaries", "[ion][loader]") {
    SpeciesDatabase species_db;
    
    // H3O+
    SpeciesProperties h3o_plus;
    h3o_plus.mass_amu = 19.02;
    h3o_plus.mass_kg = 19.02 * 1.66054e-27;
    h3o_plus.charge = 1;
    h3o_plus.charge_C = 1.602176634e-19;
    h3o_plus.mobility_m2Vs = 2.8e-4;
    h3o_plus.CCS_m2 = 1.1e-19;
    species_db.species["H3O+"] = h3o_plus;
    
    // H5O2+
    SpeciesProperties h5o2_plus;
    h5o2_plus.mass_amu = 37.03;
    h5o2_plus.mass_kg = 37.03 * 1.66054e-27;
    h5o2_plus.charge = 1;
    h5o2_plus.charge_C = 1.602176634e-19;
    h5o2_plus.mobility_m2Vs = 2.1e-4;
    h5o2_plus.CCS_m2 = 1.45e-19;
    species_db.species["H5O2+"] = h5o2_plus;
    
    IonConfig config;
    
    // Species 1: H3O+ at x=0
    IonSpeciesConfig spec1;
    spec1.species_id = "H3O+";
    spec1.count = 100;
    spec1.position.type = PositionDistribution::Point;
    spec1.position.center = {0, 0, 0};
    spec1.velocity.type = VelocityDistribution::Fixed;
    spec1.velocity.value = {100, 0, 0};
    config.species.push_back(spec1);
    
    // Species 2: H5O2+ at x=0.01
    IonSpeciesConfig spec2;
    spec2.species_id = "H5O2+";
    spec2.count = 50;
    spec2.position.type = PositionDistribution::Point;
    spec2.position.center = {0.01, 0, 0};
    spec2.velocity.type = VelocityDistribution::Fixed;
    spec2.velocity.value = {-50, 0, 0};
    config.species.push_back(spec2);
    
    std::mt19937 rng(42);
    auto result = IonLoader::generate_ions(config, species_db, rng);
    
    REQUIRE(result.success());
    REQUIRE(result.ions.size() == 150);
    
    // Check species separation
    int h3o_count = 0, h5o2_count = 0;
    for (const auto& ion : result.ions) {
        if (ion.species_id == "H3O+") {
            h3o_count++;
            REQUIRE_THAT(ion.pos.x, Catch::Matchers::WithinAbs(0.0, 1e-10));
            REQUIRE_THAT(ion.vel.x, Catch::Matchers::WithinAbs(100.0, 1e-10));
        } else if (ion.species_id == "H5O2+") {
            h5o2_count++;
            REQUIRE_THAT(ion.pos.x, Catch::Matchers::WithinAbs(0.01, 1e-10));
            REQUIRE_THAT(ion.vel.x, Catch::Matchers::WithinAbs(-50.0, 1e-10));
        }
    }
    
    REQUIRE(h3o_count == 100);
    REQUIRE(h5o2_count == 50);
}

// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"
#include <vector>

using namespace ICARION::core;
using Catch::Approx;

// Helper to create test ions
static std::vector<IonState> create_test_ions(size_t n) {
    std::vector<IonState> ions;
    for (size_t i = 0; i < n; ++i) {
        IonState ion;
        ion.pos = {i * 0.001, i * 0.002, i * 0.003};
        ion.vel = {i * 0.1, i * 0.2, i * 0.3};
        ion.mass_kg = 3.0e-26 + i * 1e-28;
        ion.ion_charge_C = 1.6e-19;
        ion.CCS_m2 = 1.5e-18;
        ion.reduced_mobility_cm2_Vs = 2.0;
        ion.species_id = (i % 2 == 0) ? "H3O+" : "H2O";
        ion.birth_time_s = i * 1e-6;
        ion.current_domain_index = 0;
        ion.history_index = static_cast<int>(i);
        ion.t = i * 1e-9;
        ion.active = (i < 90);  // Last 10 inactive
        ion.born = true;
        ions.push_back(ion);
    }
    return ions;
}

TEST_CASE("IonEnsemble construction from legacy", "[ion_ensemble][soa]") {
    auto test_ions = create_test_ions(100);
    auto ensemble = IonEnsemble::from_legacy(test_ions);
    
    REQUIRE(ensemble.size() == 100);
    REQUIRE_FALSE(ensemble.empty());
    
    // Check first ion
    CHECK(ensemble.get_pos(0).x == Approx(0.0));
    CHECK(ensemble.get_pos(0).y == Approx(0.0));
    CHECK(ensemble.get_vel(0).x == Approx(0.0));
    
    // Check middle ion
    CHECK(ensemble.get_pos(50).x == Approx(50 * 0.001));
    CHECK(ensemble.get_vel(50).y == Approx(50 * 0.2));
    
    // Check active status
    CHECK(ensemble.is_active(0));
    CHECK(ensemble.is_active(89));
    CHECK_FALSE(ensemble.is_active(90));
    CHECK_FALSE(ensemble.is_active(99));
}

TEST_CASE("IonEnsemble conversion to legacy", "[ion_ensemble][soa]") {
    auto test_ions = create_test_ions(100);
    auto ensemble = IonEnsemble::from_legacy(test_ions);
    auto converted = ensemble.to_legacy();
    
    REQUIRE(converted.size() == test_ions.size());
    
    for (size_t i = 0; i < test_ions.size(); ++i) {
        CHECK(converted[i].pos.x == Approx(test_ions[i].pos.x));
        CHECK(converted[i].pos.y == Approx(test_ions[i].pos.y));
        CHECK(converted[i].vel.x == Approx(test_ions[i].vel.x));
        CHECK(converted[i].mass_kg == Approx(test_ions[i].mass_kg));
        CHECK(converted[i].species_id == test_ions[i].species_id);
        CHECK(converted[i].active == test_ions[i].active);
    }
}

TEST_CASE("IonEnsemble compact inactive", "[ion_ensemble][soa]") {
    auto test_ions = create_test_ions(100);
    auto ensemble = IonEnsemble::from_legacy(test_ions);
    
    REQUIRE(ensemble.size() == 100);
    
    size_t removed = ensemble.compact_inactive();
    
    CHECK(removed == 10);  // 10 inactive ions
    CHECK(ensemble.size() == 90);
    
    // All remaining should be active
    for (size_t i = 0; i < ensemble.size(); ++i) {
        CHECK(ensemble.is_active(i));
    }
}

TEST_CASE("IonEnsemble advance positions", "[ion_ensemble][soa]") {
    auto test_ions = create_test_ions(100);
    auto ensemble = IonEnsemble::from_legacy(test_ions);
    
    double dt = 1e-9;
    Vec3 pos_before = ensemble.get_pos(50);
    Vec3 vel = ensemble.get_vel(50);
    
    ensemble.advance_positions(dt);
    
    Vec3 pos_after = ensemble.get_pos(50);
    
    CHECK(pos_after.x == Approx(pos_before.x + vel.x * dt));
    CHECK(pos_after.y == Approx(pos_before.y + vel.y * dt));
    CHECK(pos_after.z == Approx(pos_before.z + vel.z * dt));
}

TEST_CASE("IonEnsemble memory footprint", "[ion_ensemble][soa]") {
    auto test_ions = create_test_ions(100);
    auto ensemble = IonEnsemble::from_legacy(test_ions);
    
    size_t memory = ensemble.memory_footprint();
    size_t per_ion = memory / ensemble.size();
    
    // Should be much less than AoS (220 bytes/ion)
    CHECK(per_ion < 150);
    
    // Should be roughly 120 bytes/ion as designed
    CHECK(per_ion > 100);
    CHECK(per_ion < 140);
    
    INFO("Memory per ion: " << per_ion << " bytes (AoS: 220 bytes)");
    INFO("Reduction: " << (1.0 - per_ion / 220.0) * 100 << "%");
}

TEST_CASE("IonEnsemble empty ensemble", "[ion_ensemble][soa]") {
    IonEnsemble ensemble;
    
    CHECK(ensemble.size() == 0);
    CHECK(ensemble.empty());
    CHECK(ensemble.memory_footprint() == 0);
}

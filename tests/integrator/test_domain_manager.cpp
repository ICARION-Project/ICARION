// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file test_domain_manager.cpp
 * @brief Unit tests for DomainManager (Catch2)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/integrator/DomainManager.h"
#include "core/param/paramUtils.h"
#include "core/types/Vec3.h"

using namespace ICARION::integrator;
using namespace ICARION;
using Catch::Approx;

InstrumentDomain create_test_domain(int index, const Vec3& origin, double length, 
                                    double radius, double end_aperture = -1.0) {
    InstrumentDomain dom;
    dom.index = index;
    dom.geom.origin_m = origin;
    dom.geom.length_m = length;
    dom.geom.radius_m = radius;
    dom.geom.end_aperture_m = end_aperture;
    
    // Identity rotation
    dom.rotation_global_to_local = Mat3::identity();
    dom.rotation_local_to_global = Mat3::identity();
    
    // Environment
    dom.env.neutral_mass_kg = 6.646e-27;  // He-4
    dom.env.temperature_K = 300.0;
    dom.env.particle_density_m_3 = 2.4e25;
    dom.env.gas_velocity_m_s = {0.0, 0.0, 0.0};
    
    return dom;
}

TEST_CASE("DomainManager: Constructor rejects empty domains") {
    std::vector<InstrumentDomain> empty;
    REQUIRE_THROWS(DomainManager(empty));
}

TEST_CASE("DomainManager: Constructor accepts valid domains") {
    std::vector<InstrumentDomain> domains = {create_test_domain(0, Vec3{0,0,0}, 0.1, 0.05)};
    REQUIRE_NOTHROW(DomainManager(domains));
}

TEST_CASE("DomainManager: Find domain by position") {
    std::vector<InstrumentDomain> domains;
    domains.push_back(create_test_domain(0, Vec3{0,0,0}, 0.1, 0.05));
    domains.push_back(create_test_domain(1, Vec3{0,0,0.1}, 0.2, 0.05, 0.02));
    DomainManager manager(domains);
    
    SECTION("Inside domain 0") {
        REQUIRE(manager.find_domain_index(Vec3{0,0,0.05}) == 0);
    }
    
    SECTION("Inside domain 1") {
        REQUIRE(manager.find_domain_index(Vec3{0,0,0.15}) == 1);
    }
    
    SECTION("Outside all domains") {
        REQUIRE(manager.find_domain_index(Vec3{0,0,-0.1}) == -1);
        REQUIRE(manager.find_domain_index(Vec3{0,0,0.5}) == -1);
        REQUIRE(manager.find_domain_index(Vec3{1.0,0,0.05}) == -1);
    }
}

TEST_CASE("DomainManager: Coordinate transforms") {
    std::vector<InstrumentDomain> domains;
    domains.push_back(create_test_domain(0, Vec3{0,0,0}, 0.1, 0.05));
    domains.push_back(create_test_domain(1, Vec3{0,0,0.1}, 0.2, 0.05));
    DomainManager manager(domains);
    
    SECTION("Global to local (identity rotation)") {
        Vec3 pos_global{0.01, 0.02, 0.03};
        Vec3 pos_local = manager.global_to_local_pos(pos_global, 0);
        REQUIRE(pos_local.x == Approx(0.01).margin(1e-9));
        REQUIRE(pos_local.y == Approx(0.02).margin(1e-9));
        REQUIRE(pos_local.z == Approx(0.03).margin(1e-9));
    }
    
    SECTION("Global to local (with origin offset)") {
        Vec3 pos_global{0.0, 0.0, 0.15};
        Vec3 pos_local = manager.global_to_local_pos(pos_global, 1);
        REQUIRE(pos_local.z == Approx(0.05).margin(1e-9));
    }
    
    SECTION("Inverse transforms") {
        Vec3 pos_orig{0.01, 0.02, 0.15};
        Vec3 pos_local = manager.global_to_local_pos(pos_orig, 1);
        Vec3 pos_back = manager.local_to_global_pos(pos_local, 1);
        REQUIRE(pos_back.x == Approx(pos_orig.x).margin(1e-9));
        REQUIRE(pos_back.y == Approx(pos_orig.y).margin(1e-9));
        REQUIRE(pos_back.z == Approx(pos_orig.z).margin(1e-9));
    }
}

TEST_CASE("DomainManager: Aperture crossing") {
    std::vector<InstrumentDomain> domains;
    domains.push_back(create_test_domain(0, Vec3{0,0,0}, 0.1, 0.05));  // No aperture
    domains.push_back(create_test_domain(1, Vec3{0,0,0.1}, 0.2, 0.05, 0.02));  // 2cm aperture
    DomainManager manager(domains);
    
    SECTION("No aperture - allows all crossings") {
        IonState ion;
        ion.active = true;
        manager.check_aperture_crossing(ion, 0, Vec3{0,0,0.09}, Vec3{0,0,0.11});
        REQUIRE(ion.active);
    }
    
    SECTION("Inside aperture - allows crossing") {
        IonState ion;
        ion.active = true;
        manager.check_aperture_crossing(ion, 1, Vec3{0,0,0.19}, Vec3{0,0,0.21});
        REQUIRE(ion.active);
    }
    
    SECTION("Outside aperture - blocks ion") {
        IonState ion;
        ion.active = true;
        manager.check_aperture_crossing(ion, 1, Vec3{0.03,0,0.19}, Vec3{0.03,0,0.21});
        REQUIRE_FALSE(ion.active);
    }
}

TEST_CASE("DomainManager: Update ion properties") {
    std::vector<InstrumentDomain> domains;
    domains.push_back(create_test_domain(0, Vec3{0,0,0}, 0.1, 0.05));
    domains.push_back(create_test_domain(1, Vec3{0,0,0.1}, 0.2, 0.05));
    DomainManager manager(domains);
    
    IonState ion;
    ion.current_domain_index = -1;
    
    manager.update_domain_properties(ion, 1);
    
    REQUIRE(ion.current_domain_index == 1);
    REQUIRE(ion.domain_temperature_K == Approx(300.0).margin(1e-6));
    REQUIRE(ion.domain_particle_density_m3 == Approx(2.4e25).margin(1e20));
}

TEST_CASE("DomainManager: Full workflow") {
    std::vector<InstrumentDomain> domains;
    domains.push_back(create_test_domain(0, Vec3{0,0,0}, 0.1, 0.05));
    domains.push_back(create_test_domain(1, Vec3{0,0,0.1}, 0.2, 0.05));
    DomainManager manager(domains);
    
    IonState ion;
    ion.pos = Vec3{0,0,0.05};
    ion.vel = Vec3{0,0,1.0};
    ion.active = true;
    
    // Find initial domain
    int idx = manager.find_domain_index(ion.pos);
    REQUIRE(idx == 0);
    
    // Update properties
    manager.update_domain_properties(ion, idx);
    REQUIRE(ion.current_domain_index == 0);
    
    // Simulate motion
    ion.pos.z = 0.15;  // Now in domain 1
    int new_idx = manager.find_domain_index(ion.pos);
    REQUIRE(new_idx == 1);
    
    manager.update_domain_properties(ion, new_idx);
    REQUIRE(ion.current_domain_index == 1);
}

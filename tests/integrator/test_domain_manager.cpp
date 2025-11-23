// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file test_domain_manager.cpp
 * @brief Unit tests for DomainManager (Catch2)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/integrator/DomainManager.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/Vec3.h"

using namespace ICARION::integrator;
using namespace ICARION;
using Catch::Approx;

config::DomainConfig create_test_domain(int index, const Vec3& origin, double length, 
                                         double radius, double end_aperture = -1.0) {
    config::DomainConfig dom;
    dom.domain_index = index;
    dom.name = "test_domain_" + std::to_string(index);
    dom.instrument = config::Instrument::IMS;  // Default test instrument
    
    // Geometry (SSOT accessor)
    dom.geometry.origin_m = origin;
    dom.geometry.length_m = length;
    dom.geometry.radius_m = radius;
    dom.geometry.end_aperture_m = end_aperture;
    
    // Identity rotation
    dom.rotation_global_to_local = Mat3::identity();
    dom.rotation_local_to_global = Mat3::identity();
    
    // Environment (SSOT accessor)
    dom.environment.gas_species = "He";
    dom.environment.temperature_K = 300.0;
    dom.environment.pressure_Pa = 101325.0;
    dom.environment.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    dom.environment.compute_derived_properties();  // Compute gas_mass_kg, particle_density_m_3
    
    return dom;
}

TEST_CASE("DomainManager: Constructor rejects empty domains") {
    std::vector<config::DomainConfig> empty;
    REQUIRE_THROWS(DomainManager(empty));
}

TEST_CASE("DomainManager: Constructor accepts valid domains") {
    std::vector<config::DomainConfig> domains = {create_test_domain(0, Vec3{0,0,0}, 0.1, 0.05)};
    REQUIRE_NOTHROW(DomainManager(domains));
}

TEST_CASE("DomainManager: Find domain by position") {
    std::vector<config::DomainConfig> domains;
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
    std::vector<config::DomainConfig> domains;
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
    std::vector<config::DomainConfig> domains;
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
    std::vector<config::DomainConfig> domains;
    domains.push_back(create_test_domain(0, Vec3{0,0,0}, 0.1, 0.05));
    domains.push_back(create_test_domain(1, Vec3{0,0,0.1}, 0.2, 0.05));
    DomainManager manager(domains);
    
    IonState ion;
    ion.current_domain_index = -1;
    
    manager.update_domain_properties(ion, 1);
    
    REQUIRE(ion.current_domain_index == 1);
    REQUIRE(ion.domain_temperature_K == Approx(300.0).margin(1e-6));
    // particle_density is computed from P/(kB*T), approximately 2.44e25 for 1 atm at 300K
    REQUIRE(ion.domain_particle_density_m3 == Approx(2.44e25).margin(1e23));
    REQUIRE(ion.domain_neutral_mass_kg == Approx(6.646e-27).margin(1e-30));  // He-4 mass
}

TEST_CASE("DomainManager: Full workflow") {
    std::vector<config::DomainConfig> domains;
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

TEST_CASE("DomainManager: Orbitrap hyperbolic geometry") {
    config::DomainConfig orbitrap;
    orbitrap.domain_index = 0;
    orbitrap.instrument = config::Instrument::Orbitrap;
    
    // Typical Orbitrap dimensions (Thermo Q Exactive)
    // Note: r_char > r_out for realistic Orbitrap geometry!
    orbitrap.geometry.radius_in_m = 0.012;   // 12 mm inner electrode
    orbitrap.geometry.radius_out_m = 0.020;  // 20 mm outer electrode  
    orbitrap.geometry.radius_char_m = 0.025; // 25 mm characteristic radius (> r_out!)
    orbitrap.geometry.origin_m = Vec3{0,0,0};
    
    orbitrap.rotation_global_to_local = Mat3::identity();
    orbitrap.rotation_local_to_global = Mat3::identity();
    
    orbitrap.environment.gas_species = "N2";
    orbitrap.environment.temperature_K = 300.0;
    orbitrap.environment.pressure_Pa = 1e-8;  // Ultra-high vacuum
    orbitrap.environment.compute_derived_properties();
    
    std::vector<config::DomainConfig> domains = {orbitrap};
    DomainManager manager(domains);
    
    SECTION("Ion inside Orbitrap (between hyperbolas)") {
        // At z=0, allowed radii: r_in = 12 mm, r_out = 20 mm
        REQUIRE(manager.find_domain_index(Vec3{0.015, 0, 0}) == 0);  // r=15mm, inside
        REQUIRE(manager.find_domain_index(Vec3{0.017, 0, 0}) == 0);  // r=17mm, inside
    }
    
    SECTION("Ion too close to inner electrode") {
        REQUIRE(manager.find_domain_index(Vec3{0.010, 0, 0}) == -1);  // r=10mm < r_in(0)
    }
    
    SECTION("Ion outside outer electrode") {
        REQUIRE(manager.find_domain_index(Vec3{0.025, 0, 0}) == -1);  // r=25mm > r_out(0)
    }
    
    SECTION("Hyperbolic constraint at z != 0") {
        // The electrode surfaces r_in(z) and r_out(z) depend on z via the
        // logarithmic-hyperbolic equation: z² = 0.5·(r² - R²) + R_m² · ln(R/r)
        // For r_char=25mm > r_out=20mm, this z-dependence is very weak
        // (gap changes by <0.01% for z=0..20mm), but it's NOT zero!
        double z = 0.007;
        REQUIRE(manager.find_domain_index(Vec3{0.015, 0, z}) == 0);  // Still inside
        REQUIRE(manager.find_domain_index(Vec3{0.011, 0, z}) == -1); // Too close to inner
    }
}

TEST_CASE("DomainManager: Orbitrap boundary termination uses midpoint approximation") {
    config::DomainConfig orbitrap;
    orbitrap.domain_index = 0;
    orbitrap.instrument = config::Instrument::Orbitrap;
    orbitrap.geometry.radius_in_m = 0.012;
    orbitrap.geometry.radius_out_m = 0.020;
    orbitrap.geometry.radius_char_m = 0.025;  // r_char > r_out (realistic)
    orbitrap.geometry.origin_m = Vec3{0,0,0};
    orbitrap.rotation_global_to_local = Mat3::identity();
    orbitrap.rotation_local_to_global = Mat3::identity();
    orbitrap.environment.temperature_K = 300.0;
    orbitrap.environment.pressure_Pa = 1e-8;
    orbitrap.environment.compute_derived_properties();
    
    std::vector<config::DomainConfig> domains = {orbitrap};
    DomainManager manager(domains);
    
    IonState ion;
    ion.active = true;
    ion.vel = Vec3{0, 0, 10.0};  // Will be set to zero
    
    Vec3 pos_before{0.015, 0, 0};  // Inside
    Vec3 pos_after{0.025, 0, 0};   // Outside (hit outer electrode)
    
    manager.terminate_ion_at_boundary(ion, 0, pos_before, pos_after);
    
    // Should use midpoint approximation for Orbitrap
    Vec3 expected_midpoint = (pos_before + pos_after) * 0.5;
    REQUIRE(ion.pos.x == Approx(expected_midpoint.x).margin(1e-9));
    REQUIRE(ion.pos.y == Approx(expected_midpoint.y).margin(1e-9));
    REQUIRE(ion.pos.z == Approx(expected_midpoint.z).margin(1e-9));
    
    // Velocity should be zero
    REQUIRE(ion.vel.x == 0.0);
    REQUIRE(ion.vel.y == 0.0);
    REQUIRE(ion.vel.z == 0.0);
    
    // Ion should be deactivated
    REQUIRE_FALSE(ion.active);
}

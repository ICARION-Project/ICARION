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
    
    // Orbitrap hyperbolic boundary configuration
    // Hyperboloid equation: z² - r²/2 = C, where C = -0.5 * R²
    orbitrap.geometry.radius_in_m = 0.010;   // 10 mm inner electrode
    orbitrap.geometry.radius_out_m = 0.015;  // 15 mm outer electrode  
    orbitrap.geometry.radius_char_m = 0.020; // Characteristic radius (used for fields)
    orbitrap.geometry.origin_m = Vec3{0,0,0};
    
    // Compute hyperbolic boundary constants
    orbitrap.geometry.orbitrap_C_in  = -0.5 * 0.010 * 0.010;  // -5e-5
    orbitrap.geometry.orbitrap_C_out = -0.5 * 0.015 * 0.015;  // -1.125e-4
    
    orbitrap.rotation_global_to_local = Mat3::identity();
    orbitrap.rotation_local_to_global = Mat3::identity();
    
    orbitrap.environment.gas_species = "N2";
    orbitrap.environment.temperature_K = 300.0;
    orbitrap.environment.pressure_Pa = 1e-8;  // Ultra-high vacuum
    orbitrap.environment.compute_derived_properties();
    
    std::vector<config::DomainConfig> domains = {orbitrap};
    DomainManager manager(domains);
    
    SECTION("Ion inside Orbitrap (between hyperbolas)") {
        // At z=0: C = -r²/2, so inside if C_in <= -r²/2 <= C_out
        // -5e-5 <= -r²/2 <= -1.125e-4  =>  sqrt(2*1.125e-4) <= r <= sqrt(2*5e-5)
        // 0.015 <= r <= 0.010  (WAIT, that's backwards! Let me recalculate)
        // Actually: -r²/2 must be between C_in and C_out
        // C_in < C_out (both negative, C_in is MORE negative)
        // -5e-5 < -1.125e-4 is FALSE! So C_in = -1.125e-4, C_out = -5e-5
        // Wait, let me recalculate: C = -0.5*R² means SMALLER R gives LESS negative C
        // So C_in (inner, smaller R) = -0.5*0.01² = -5e-5 (less negative)
        // C_out (outer, larger R) = -0.5*0.015² = -1.125e-4 (more negative)
        // Inside means: C_out <= C <= C_in (from more negative to less negative)
        REQUIRE(manager.find_domain_index(Vec3{0.012, 0, 0}) == 0);  // r=12mm, inside
        REQUIRE(manager.find_domain_index(Vec3{0.013, 0, 0}) == 0);  // r=13mm, inside
    }
    
    SECTION("Ion too close to inner electrode") {
        REQUIRE(manager.find_domain_index(Vec3{0.009, 0, 0}) == -1);  // r=9mm < r_in
    }
    
    SECTION("Ion outside outer electrode") {
        REQUIRE(manager.find_domain_index(Vec3{0.016, 0, 0}) == -1);  // r=16mm > r_out
    }
    
    SECTION("Hyperbolic constraint at z != 0") {
        // At z != 0, hyperboloid allows larger radii for larger |z|
        // Example: At z=5mm, compute allowed r range:
        // C_in: z² - r²/2 = -5e-5  →  r² = 2*(z² + 5e-5) = 2*(2.5e-5 + 5e-5) = 1.5e-4  →  r ≈ 12.25mm
        // C_out: z² - r²/2 = -1.125e-4  →  r² = 2*(z² + 1.125e-4) = 2*(2.5e-5 + 1.125e-4) = 2.75e-4  →  r ≈ 16.58mm
        // So at z=5mm, allowed range is approximately 12.25mm <= r <= 16.58mm
        double z = 0.005;  // 5mm
        
        // r=13mm should be inside (13 > 12.25)
        REQUIRE(manager.find_domain_index(Vec3{0.013, 0, z}) == 0);
        
        // r=12mm should be outside (12 < 12.25)  
        REQUIRE(manager.find_domain_index(Vec3{0.012, 0, z}) == -1);
        
        // r=17mm should be outside (17 > 16.58)
        REQUIRE(manager.find_domain_index(Vec3{0.017, 0, z}) == -1);
    }
}

TEST_CASE("DomainManager: Orbitrap boundary termination uses analytical ray-hyperboloid intersection") {
    config::DomainConfig orbitrap;
    orbitrap.domain_index = 0;
    orbitrap.instrument = config::Instrument::Orbitrap;
    orbitrap.geometry.radius_in_m = 0.010;
    orbitrap.geometry.radius_out_m = 0.015;
    orbitrap.geometry.radius_char_m = 0.020;
    orbitrap.geometry.origin_m = Vec3{0,0,0};
    
    // Compute hyperbolic boundary constants
    orbitrap.geometry.orbitrap_C_in  = -0.5 * 0.010 * 0.010;
    orbitrap.geometry.orbitrap_C_out = -0.5 * 0.015 * 0.015;
    
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
    
    Vec3 pos_before{0.012, 0, 0};  // Inside (r=12mm)
    Vec3 pos_after{0.020, 0, 0};   // Outside (r=20mm, hit outer electrode at r=15mm)
    
    manager.terminate_ion_at_boundary(ion, 0, pos_before, pos_after);
    
    // Should use analytical ray-hyperboloid intersection
    // Ray travels from r=12mm to r=20mm along x-axis (z=0)
    // At z=0: C = -r²/2, boundary is at r=15mm where C = C_out = -1.125e-4
    // Intersection should be at approximately (0.015, 0, 0)
    REQUIRE(ion.pos.x == Approx(0.015).margin(1e-6));
    REQUIRE(ion.pos.y == Approx(0.0).margin(1e-9));
    REQUIRE(ion.pos.z == Approx(0.0).margin(1e-9));
    
    // Velocity should be zero
    REQUIRE(ion.vel.x == 0.0);
    REQUIRE(ion.vel.y == 0.0);
    REQUIRE(ion.vel.z == 0.0);
    
    // Ion should be deactivated
    REQUIRE_FALSE(ion.active);
}

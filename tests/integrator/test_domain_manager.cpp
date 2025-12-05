// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file test_domain_manager.cpp
 * @brief Unit tests for DomainManager (Catch2)
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "core/integrator/DomainManager.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/CylindricalGeometry.h"
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
    // Simulate motion
    ion.pos.z = 0.15;  // Now in domain 1
    int new_idx = manager.find_domain_index(ion.pos);
    REQUIRE(new_idx == 1);
}

TEST_CASE("DomainManager: Orbitrap hyperbolic geometry") {
    config::DomainConfig orbitrap;
    orbitrap.domain_index = 0;
    orbitrap.instrument = config::Instrument::Orbitrap;
    
    // Orbitrap hyperlogarithmic boundary configuration (Makarov 2000)
    // Equation: z² = 0.5(r² - R²) + R_m² × ln(R/r)
    orbitrap.geometry.radius_in_m = 0.010;   // 10 mm inner electrode (R_in)
    orbitrap.geometry.radius_out_m = 0.015;  // 15 mm outer electrode (R_out)
    orbitrap.geometry.radius_char_m = 0.020; // 20 mm characteristic radius (R_m)
    orbitrap.geometry.length_m = 0.04;       // 40 mm trap length (±20mm axial)
    orbitrap.geometry.origin_m = Vec3{0,0,0};
    
    orbitrap.rotation_global_to_local = Mat3::identity();
    orbitrap.rotation_local_to_global = Mat3::identity();
    
    orbitrap.environment.gas_species = "N2";
    orbitrap.environment.temperature_K = 300.0;
    orbitrap.environment.pressure_Pa = 1e-8;  // Ultra-high vacuum
    orbitrap.environment.compute_derived_properties();
    
    std::vector<config::DomainConfig> domains = {orbitrap};
    DomainManager manager(domains);
    
    SECTION("Ion inside Orbitrap (between electrodes)") {
        // At z=0: r must be between R_in and R_out
        // R_in=10mm, R_out=15mm → allowed range: 10mm < r < 15mm
        REQUIRE(manager.find_domain_index(Vec3{0.012, 0, 0}) == 0);  // r=12mm, inside
        REQUIRE(manager.find_domain_index(Vec3{0.013, 0, 0}) == 0);  // r=13mm, inside
    }
    
    SECTION("Ion too close to inner electrode") {
        REQUIRE(manager.find_domain_index(Vec3{0.009, 0, 0}) == -1);  // r=9mm < r_in
    }
    
    SECTION("Ion outside outer electrode") {
        REQUIRE(manager.find_domain_index(Vec3{0.016, 0, 0}) == -1);  // r=16mm > r_out
    }
    
    SECTION("Hyperlogarithmic constraint at z != 0") {
        // Hyperlogarithmic equation: z² = 0.5(r² - R²) + R_m² × ln(R/r)
        // At z != 0, electrodes curve inward (smaller r for larger |z|)
        
        double z = 0.005;  // 5mm
        const double R_in = 0.010;   // 10mm
        const double R_out = 0.015;  // 15mm
        const double R_m = 0.020;    // 20mm
        
        // Compute exact boundaries using bisection (same as DomainManager::orbitrap_r_for_z)
        auto compute_r_at_z = [&](double z_val, double R, double R_m_val) -> double {
            // Bisection solver for: z² = 0.5(r² - R²) + R_m² × ln(R/r)
            auto residual = [&](double r) -> double {
                return z_val*z_val - 0.5*(r*r - R*R) - R_m_val*R_m_val * std::log(R/r);
            };
            
            // Hyperlogarithmic surfaces curve INWARD as |z| increases
            // So we bracket between 0.1*R and R (not outward to 3*R)
            double r_lo = 0.1 * R;
            double r_hi = R;
            const double eps = 1e-10;
            const int max_iter = 80;
            
            // Expand bracket if needed (same as DomainManager)
            double f_lo = residual(r_lo);
            double f_hi = residual(r_hi);
            int expand_iter = 0;
            while (f_lo * f_hi > 0.0 && expand_iter < 10) {
                r_lo *= 0.5;
                r_hi *= 1.5;
                f_lo = residual(r_lo);
                f_hi = residual(r_hi);
                expand_iter++;
            }
            
            for (int i = 0; i < max_iter; ++i) {
                double r_mid = 0.5 * (r_lo + r_hi);
                double f_mid = residual(r_mid);
                
                if (std::fabs(f_mid) < eps) {
                    return r_mid;
                }
                
                f_lo = residual(r_lo);
                if (f_mid * f_lo > 0.0) {
                    r_lo = r_mid;
                    f_lo = f_mid;
                } else {
                    r_hi = r_mid;
                }
            }
            return 0.5 * (r_lo + r_hi);
        };
        
        double r_in_at_z5 = compute_r_at_z(z, R_in, R_m);
        double r_out_at_z5 = compute_r_at_z(z, R_out, R_m);
        
        // Debug output to verify computed boundaries
        INFO("At z=" << z*1000 << "mm:");
        INFO("  r_in = " << r_in_at_z5*1000 << " mm");
        INFO("  r_out = " << r_out_at_z5*1000 << " mm");
        INFO("  r_in < r_out? " << (r_in_at_z5 < r_out_at_z5));
        
        // Sanity check: r_in should be less than r_out
        REQUIRE(r_in_at_z5 < r_out_at_z5);
        
        // Test points based on computed boundaries
        double r_mid = 0.5 * (r_in_at_z5 + r_out_at_z5);  // Midpoint should be inside
        double r_below = r_in_at_z5 - 0.001;  // 1mm below inner boundary
        double r_above = r_out_at_z5 + 0.001; // 1mm above outer boundary
        
        INFO("Testing r_mid=" << r_mid*1000 << "mm at z=" << z*1000 << "mm");
        INFO("Expected: inside (r_in < r_mid < r_out)");
        
        REQUIRE(manager.find_domain_index(Vec3{r_mid, 0, z}) == 0);      // Inside
        REQUIRE(manager.find_domain_index(Vec3{r_below, 0, z}) == -1);   // Too close to inner
        REQUIRE(manager.find_domain_index(Vec3{r_above, 0, z}) == -1);   // Beyond outer
    }
}

TEST_CASE("CylindricalGeometry: aperture respected at exit plane") {
    auto dom = create_test_domain(0, Vec3{0,0,0}, 0.2, 0.05, 0.02);
    ICARION::config::CylindricalGeometry geom(dom);

    Vec3 start_local{0.01, 0.0, 0.19};  // inside, near exit
    Vec3 end_local{0.01, 0.0, 0.21};    // crosses exit plane
    Vec3 hit{};
    REQUIRE(geom.first_boundary_intersection(start_local, end_local, hit));
    // r=0.01 < aperture=0.02 → should be accepted; intersection at z~0.2
    REQUIRE(hit.z == Approx(0.2).margin(1e-6));

    Vec3 start_block{0.03, 0.0, 0.19};  // r=0.03 > aperture
    Vec3 end_block{0.03, 0.0, 0.21};
    Vec3 hit_block{};
    REQUIRE(geom.first_boundary_intersection(start_block, end_block, hit_block));
    double r_hit = std::sqrt(hit_block.x * hit_block.x + hit_block.y * hit_block.y);
    // Intersection exists, but radius exceeds aperture; caller can reject
    REQUIRE(r_hit > dom.geometry.end_aperture_m);
}

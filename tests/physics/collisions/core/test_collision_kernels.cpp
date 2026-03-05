// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file test_collision_kernels.cpp
 * @brief Unit tests for CollisionKernels module
 * 
 * Tests critical collision physics:
 * - EHSS collision (geometry-resolved)
 * - HSS collision (isotropic)
 * - OU velocity updates (Langevin dynamics)
 * - Momentum conservation
 * - Energy conservation
 * - Regression against old implementation
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/physics/collisions/core/CollisionKernels.h"
#include "core/utils/mathUtils.h"
#include "utils/constants.h"
#include <cmath>
#include <vector>

using namespace ICARION::physics::collision_core;
using namespace ICARION::core;
using ICARION::physics::PhysicsRng;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

constexpr double TOL = 1e-10;

// Helper: Compute kinetic energy
double kinetic_energy(const Vec3& v, double mass) {
    return 0.5 * mass * dot(v, v);
}

// Helper: Compute total momentum
Vec3 total_momentum(const Vec3& v1, double m1, const Vec3& v2, double m2) {
    return v1 * m1 + v2 * m2;
}

TEST_CASE("CollisionKernels: HSS collision basics", "[collision][kernels][hss]") {
    SECTION("Momentum conservation") {
        PhysicsRng rng(42);
        
        Vec3 v_ion{100.0, 50.0, -20.0};
        Vec3 v_neutral{10.0, -5.0, 15.0};
        double m_ion = 28.0 * AMU_TO_KG;
        double m_neutral = 28.0 * AMU_TO_KG;
        
        Vec3 p_before = total_momentum(v_ion, m_ion, v_neutral, m_neutral);
        
        Vec3 v_ion_after = CollisionKernels::hss_collision(v_ion, v_neutral, m_ion, m_neutral, rng);
        
        // Momentum conservation (assuming neutral unchanged in collision frame)
        Vec3 p_after = total_momentum(v_ion_after, m_ion, v_neutral, m_neutral);
        
        REQUIRE_THAT(p_after.x, WithinAbs(p_before.x, 1e-12));
        REQUIRE_THAT(p_after.y, WithinAbs(p_before.y, 1e-12));
        REQUIRE_THAT(p_after.z, WithinAbs(p_before.z, 1e-12));
    }
    
    SECTION("Energy conservation") {
        PhysicsRng rng(42);
        
        Vec3 v_ion{100.0, 50.0, -20.0};
        Vec3 v_neutral{10.0, -5.0, 15.0};
        double m_ion = 28.0 * AMU_TO_KG;
        double m_neutral = 28.0 * AMU_TO_KG;
        
        // Energy before collision
        double E_before = 0.5 * m_ion * dot(v_ion, v_ion) + 0.5 * m_neutral * dot(v_neutral, v_neutral);
        
        // Momentum before collision
        Vec3 p_before = v_ion * m_ion + v_neutral * m_neutral;
        
        Vec3 v_ion_after = CollisionKernels::hss_collision(v_ion, v_neutral, m_ion, m_neutral, rng);
        
        // Compute neutral velocity after from momentum conservation
        Vec3 v_neutral_after = (p_before - v_ion_after * m_ion) * (1.0 / m_neutral);
        
        // Energy after collision
        double E_after = 0.5 * m_ion * dot(v_ion_after, v_ion_after) + 0.5 * m_neutral * dot(v_neutral_after, v_neutral_after);
        
        // Energy should be conserved (elastic collision)
        REQUIRE_THAT(E_after, WithinRel(E_before, 1e-10));
    }
    
    SECTION("Isotropic scattering (statistical test)") {
        PhysicsRng rng(42);
        
        Vec3 v_ion{100.0, 0.0, 0.0};  // Ion moving in +x direction
        Vec3 v_neutral{0.0, 0.0, 0.0};  // Neutral at rest
        double m_ion = 28.0 * AMU_TO_KG;
        double m_neutral = 28.0 * AMU_TO_KG;
        
        // Sample many collisions
        const int N = 5000;
        Vec3 mean_v{0, 0, 0};
        
        for (int i = 0; i < N; ++i) {
            Vec3 v_after = CollisionKernels::hss_collision(v_ion, v_neutral, m_ion, m_neutral, rng);
            mean_v = mean_v + v_after;
        }
        mean_v = mean_v / N;
        
        // Mean velocity should be near COM velocity (50, 0, 0) for equal masses
        Vec3 v_com = (v_ion * m_ion + v_neutral * m_neutral) / (m_ion + m_neutral);
        
        REQUIRE_THAT(mean_v.x, WithinAbs(v_com.x, 2.0));  // Statistical tolerance
        REQUIRE_THAT(mean_v.y, WithinAbs(v_com.y, 2.0));
        REQUIRE_THAT(mean_v.z, WithinAbs(v_com.z, 2.0));
    }
    
    SECTION("Determinism with same seed") {
        PhysicsRng rng1(123);
        PhysicsRng rng2(123);
        
        Vec3 v_ion{100.0, 50.0, -20.0};
        Vec3 v_neutral{10.0, -5.0, 15.0};
        double m_ion = 28.0 * AMU_TO_KG;
        double m_neutral = 28.0 * AMU_TO_KG;
        
        Vec3 v1 = CollisionKernels::hss_collision(v_ion, v_neutral, m_ion, m_neutral, rng1);
        Vec3 v2 = CollisionKernels::hss_collision(v_ion, v_neutral, m_ion, m_neutral, rng2);
        
        REQUIRE_THAT(v1.x, WithinAbs(v2.x, TOL));
        REQUIRE_THAT(v1.y, WithinAbs(v2.y, TOL));
        REQUIRE_THAT(v1.z, WithinAbs(v2.z, TOL));
    }
}

TEST_CASE("CollisionKernels: EHSS collision basics", "[collision][kernels][ehss]") {
    SECTION("Simple single-atom geometry") {
        PhysicsRng rng(42);
        
        // Simple geometry: one atom at origin
        std::vector<Vec3> centers = {{0.0, 0.0, 0.0}};
        std::vector<double> radii = {1.5e-10};  // 1.5 Å
        
        Vec3 v_ion{100.0, 0.0, 0.0};
        Vec3 v_neutral{0.0, 0.0, 0.0};
        double m_ion = 28.0 * AMU_TO_KG;
        double m_neutral = 28.0 * AMU_TO_KG;
        double r_ion = 1.5e-10;
        
        Vec3 v_after = CollisionKernels::ehss_collision(
            v_ion, v_neutral, m_ion, m_neutral, r_ion,
            centers, radii, rng
        );
        
        // Should have collided (velocity changed)
        REQUIRE(norm(v_after - v_ion) > 1e-6);
    }
    
    SECTION("Momentum conservation with geometry") {
        PhysicsRng rng(42);
        
        // H2O-like geometry (simplified)
        std::vector<Vec3> centers = {
            {0.0, 0.0, 0.0},        // O
            {0.96e-10, 0.0, 0.0},   // H1
            {-0.24e-10, 0.93e-10, 0.0}  // H2
        };
        std::vector<double> radii = {1.4e-10, 1.2e-10, 1.2e-10};
        
        Vec3 v_ion{200.0, 50.0, -30.0};
        Vec3 v_neutral{20.0, -10.0, 5.0};
        double m_ion = 18.0 * AMU_TO_KG;  // H3O+
        double m_neutral = 18.0 * AMU_TO_KG;  // H2O
        double r_ion = 1.5e-10;
        
        Vec3 p_before = total_momentum(v_ion, m_ion, v_neutral, m_neutral);
        
        Vec3 v_after = CollisionKernels::ehss_collision(
            v_ion, v_neutral, m_ion, m_neutral, r_ion,
            centers, radii, rng
        );
        
        Vec3 p_after = total_momentum(v_after, m_ion, v_neutral, m_neutral);
        
        REQUIRE_THAT(p_after.x, WithinAbs(p_before.x, 1e-12));
        REQUIRE_THAT(p_after.y, WithinAbs(p_before.y, 1e-12));
        REQUIRE_THAT(p_after.z, WithinAbs(p_before.z, 1e-12));
    }
    
    SECTION("Energy conservation with geometry") {
        PhysicsRng rng(42);
        
        std::vector<Vec3> centers = {{0.0, 0.0, 0.0}};
        std::vector<double> radii = {1.5e-10};
        
        Vec3 v_ion{200.0, 50.0, -30.0};
        Vec3 v_neutral{20.0, -10.0, 5.0};
        double m_ion = 28.0 * AMU_TO_KG;
        double m_neutral = 28.0 * AMU_TO_KG;
        double r_ion = 1.5e-10;
        
        // Energy before collision
        double E_before = 0.5 * m_ion * dot(v_ion, v_ion) + 0.5 * m_neutral * dot(v_neutral, v_neutral);
        
        // Momentum before collision
        Vec3 p_before = v_ion * m_ion + v_neutral * m_neutral;
        
        Vec3 v_after = CollisionKernels::ehss_collision(
            v_ion, v_neutral, m_ion, m_neutral, r_ion,
            centers, radii, rng
        );
        
        // Compute neutral velocity after from momentum conservation
        Vec3 v_neutral_after = (p_before - v_after * m_ion) * (1.0 / m_neutral);
        
        // Energy after collision
        double E_after = 0.5 * m_ion * dot(v_after, v_after) + 0.5 * m_neutral * dot(v_neutral_after, v_neutral_after);
        
        // Energy should be conserved (elastic collision)
        REQUIRE_THAT(E_after, WithinRel(E_before, 1e-10));
    }
    
    SECTION("No collision with empty geometry") {
        PhysicsRng rng(42);
        
        std::vector<Vec3> centers;
        std::vector<double> radii;
        
        Vec3 v_ion{100.0, 0.0, 0.0};
        Vec3 v_neutral{0.0, 0.0, 0.0};
        double m_ion = 28.0 * AMU_TO_KG;
        double m_neutral = 28.0 * AMU_TO_KG;
        double r_ion = 1.5e-10;
        
        Vec3 v_after = CollisionKernels::ehss_collision(
            v_ion, v_neutral, m_ion, m_neutral, r_ion,
            centers, radii, rng
        );
        
        // Velocity unchanged (no collision)
        REQUIRE_THAT(v_after.x, WithinAbs(v_ion.x, TOL));
        REQUIRE_THAT(v_after.y, WithinAbs(v_ion.y, TOL));
        REQUIRE_THAT(v_after.z, WithinAbs(v_ion.z, TOL));
    }
}

TEST_CASE("CollisionKernels: OU velocity updates", "[collision][kernels][ou]") {
    SECTION("Thermalization to gas velocity (damping mode)") {
        PhysicsRng rng(42);
        
        IonState ion;
        ion.mass_kg = 28.0 * AMU_TO_KG;
        ion.vel = {500.0, 0.0, 0.0};  // Hot ion
        
        double dt = 1e-8;  // 10 ns
        double gamma = 1e6;  // 1 MHz damping
        double T = 300.0;  // K
        Vec3 u_gas{0.0, 0.0, 0.0};
        
        // Apply many OU steps
        for (int i = 0; i < 1000; ++i) {
            CollisionKernels::ou_velocity_update(ion, rng, dt, gamma, T, u_gas, true);
        }
        
        // Should have thermalized (velocity near gas velocity)
        double v_thermal = std::sqrt(BOLTZMANN_CONSTANT * T / ion.mass_kg);
        
        REQUIRE(norm(ion.vel) < 5.0 * v_thermal);  // Within 5σ
    }
    
    SECTION("Thermal kicks only (no damping)") {
        PhysicsRng rng(42);
        
        IonState ion;
        ion.mass_kg = 28.0 * AMU_TO_KG;
        ion.vel = {0.0, 0.0, 0.0};
        
        double dt = 1e-8;
        double gamma = 1e6;
        double T = 300.0;
        Vec3 u_gas{0.0, 0.0, 0.0};
        
        // Apply thermal kicks only
        CollisionKernels::ou_velocity_update(ion, rng, dt, gamma, T, u_gas, false);
        
        // Velocity should have changed (thermal noise)
        REQUIRE(norm(ion.vel) > 0.0);
    }
    
    SECTION("Determinism with same seed") {
        PhysicsRng rng1(123);
        PhysicsRng rng2(123);
        
        IonState ion1, ion2;
        ion1.mass_kg = ion2.mass_kg = 28.0 * AMU_TO_KG;
        ion1.vel = ion2.vel = {100.0, 50.0, -20.0};
        
        double dt = 1e-8;
        double gamma = 1e6;
        double T = 300.0;
        Vec3 u_gas{0.0, 0.0, 0.0};
        
        CollisionKernels::ou_velocity_update(ion1, rng1, dt, gamma, T, u_gas, true);
        CollisionKernels::ou_velocity_update(ion2, rng2, dt, gamma, T, u_gas, true);
        
        REQUIRE_THAT(ion1.vel.x, WithinAbs(ion2.vel.x, TOL));
        REQUIRE_THAT(ion1.vel.y, WithinAbs(ion2.vel.y, TOL));
        REQUIRE_THAT(ion1.vel.z, WithinAbs(ion2.vel.z, TOL));
    }
    
    SECTION("Velocity variance approaches thermal equilibrium") {
        PhysicsRng rng(42);
        
        double T = 300.0;
        double m = 28.0 * AMU_TO_KG;
        double expected_var = BOLTZMANN_CONSTANT * T / m;
        
        const int N = 5000;
        std::vector<double> vx_samples;
        vx_samples.reserve(N);
        
        for (int i = 0; i < N; ++i) {
            IonState ion;
            ion.mass_kg = m;
            ion.vel = {0.0, 0.0, 0.0};
            
            // Thermalize
            double dt = 1e-8;
            double gamma = 1e6;
            Vec3 u_gas{0.0, 0.0, 0.0};
            
            for (int j = 0; j < 100; ++j) {
                CollisionKernels::ou_velocity_update(ion, rng, dt, gamma, T, u_gas, true);
            }
            
            vx_samples.push_back(ion.vel.x);
        }
        
        // Compute variance
        double mean_vx = 0.0;
        for (double vx : vx_samples) {
            mean_vx += vx;
        }
        mean_vx /= N;
        
        double var_vx = 0.0;
        for (double vx : vx_samples) {
            var_vx += (vx - mean_vx) * (vx - mean_vx);
        }
        var_vx /= N;
        
        // Should match thermal variance (within ~20% for N=5000, 100 steps)
        // Note: More steps or longer dt would improve convergence
        REQUIRE_THAT(var_vx, WithinRel(expected_var, 0.2));
    }
}

TEST_CASE("CollisionKernels: Regression tests", "[collision][kernels][regression]") {
    SECTION("HSS collision matches old implementation") {
        // Use same seed for both
        PhysicsRng rng_old(42);
        PhysicsRng rng_new(42);
        
        Vec3 v_ion{100.0, 50.0, -20.0};
        Vec3 v_neutral{10.0, -5.0, 15.0};
        double m_ion = 28.0 * AMU_TO_KG;
        double m_neutral = 28.0 * AMU_TO_KG;
        
        // Old implementation (manual)
        Vec3 vrel_old = v_ion - v_neutral;
        double vrel_mag_old = norm(vrel_old);
        
        double u1 = rng_old.uniform01();
        double u2 = rng_old.uniform01();
        double cosT = 2.0*u1 - 1.0;
        double sinT = std::sqrt(std::max(0.0, 1.0 - cosT*cosT));
        double phi = 2.0*M_PI*u2;
        Vec3 new_dir{sinT*std::cos(phi), sinT*std::sin(phi), cosT};
        Vec3 vrel_scattered_old = new_dir * vrel_mag_old;
        
        double mt = m_ion + m_neutral;
        Vec3 Vcom_old = (v_ion * m_ion + v_neutral * m_neutral) / mt;
        Vec3 v_old = Vcom_old + vrel_scattered_old * (m_neutral / mt);
        
        // New implementation
        Vec3 v_new = CollisionKernels::hss_collision(v_ion, v_neutral, m_ion, m_neutral, rng_new);
        
        // Should be bit-for-bit identical
        REQUIRE_THAT(v_new.x, WithinAbs(v_old.x, TOL));
        REQUIRE_THAT(v_new.y, WithinAbs(v_old.y, TOL));
        REQUIRE_THAT(v_new.z, WithinAbs(v_old.z, TOL));
    }
}

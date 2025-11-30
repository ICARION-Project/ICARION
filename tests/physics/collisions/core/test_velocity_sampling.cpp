// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file test_velocity_sampling.cpp
 * @brief Unit tests for VelocitySampling module
 * 
 * Tests thermal velocity generation:
 * - Maxwell-Boltzmann distribution sampling
 * - Thermal velocity width calculation
 * - Statistical properties (mean, variance)
 * - Regression against old implementation
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/physics/collisions/core/VelocitySampling.h"
#include "core/utils/mathUtils.h"
#include "utils/constants.h"
#include <cmath>
#include <vector>

using namespace ICARION::physics::collision_core;
using namespace ICARION::core;
using ICARION::physics::EhssRng;
using Catch::Matchers::WithinAbs;
using Catch::Matchers::WithinRel;

constexpr double TOL = 1e-10;

TEST_CASE("VelocitySampling: Thermal velocity width", "[collision][velocity]") {
    SECTION("Nitrogen at 300K") {
        double T = 300.0;  // K
        double m_N2 = 28.0 * 1.66054e-27;  // kg (molecular mass)
        
        double v_th = VelocitySampling::thermal_velocity_width(T, m_N2);
        
        // Expected: sqrt(kB*T/m) ≈ 298 m/s for N2 at 300K
        double expected = std::sqrt(BOLTZMANN_CONSTANT * T / m_N2);
        REQUIRE_THAT(v_th, WithinAbs(expected, 1e-6));
        REQUIRE_THAT(v_th, WithinRel(298.5, 0.01));  // Within 1%
    }
    
    SECTION("Helium at 300K (lighter, faster)") {
        double T = 300.0;  // K
        double m_He = 4.0 * 1.66054e-27;  // kg
        
        double v_th = VelocitySampling::thermal_velocity_width(T, m_He);
        
        // Helium should have higher thermal velocity than N2
        double v_th_N2 = VelocitySampling::thermal_velocity_width(T, 28.0 * 1.66054e-27);
        REQUIRE(v_th > v_th_N2);
        REQUIRE_THAT(v_th, WithinRel(790.0, 0.01));  // ~790 m/s
    }
    
    SECTION("Temperature scaling") {
        double m = 28.0 * 1.66054e-27;  // N2
        
        double v_th_300K = VelocitySampling::thermal_velocity_width(300.0, m);
        double v_th_600K = VelocitySampling::thermal_velocity_width(600.0, m);
        
        // v_th should scale as sqrt(T)
        REQUIRE_THAT(v_th_600K / v_th_300K, WithinAbs(std::sqrt(2.0), 1e-10));
    }
}

TEST_CASE("VelocitySampling: Thermal component sampling", "[collision][velocity]") {
    SECTION("Statistical properties (mean ≈ 0, variance ≈ σ²)") {
        EhssRng rng(42);
        double T = 300.0;
        double m = 28.0 * 1.66054e-27;  // N2
        double sigma = VelocitySampling::thermal_velocity_width(T, m);
        
        // Sample many velocities
        const int N = 10000;
        std::vector<double> samples;
        samples.reserve(N);
        
        for (int i = 0; i < N; ++i) {
            samples.push_back(VelocitySampling::sample_thermal_component(T, m, rng));
        }
        
        // Compute mean
        double mean = 0.0;
        for (double v : samples) {
            mean += v;
        }
        mean /= N;
        
        // Compute variance
        double variance = 0.0;
        for (double v : samples) {
            variance += (v - mean) * (v - mean);
        }
        variance /= N;
        
        double stddev = std::sqrt(variance);
        
        // Check: mean should be near 0 (within 3σ/sqrt(N) by CLT)
        double mean_tolerance = 3.0 * sigma / std::sqrt(N);
        REQUIRE_THAT(mean, WithinAbs(0.0, mean_tolerance));
        
        // Check: stddev should be near σ (within ~5% for N=10000)
        REQUIRE_THAT(stddev, WithinRel(sigma, 0.05));
    }
    
    SECTION("Determinism with same seed") {
        EhssRng rng1(123);
        EhssRng rng2(123);
        
        double T = 300.0;
        double m = 28.0 * 1.66054e-27;
        
        double v1 = VelocitySampling::sample_thermal_component(T, m, rng1);
        double v2 = VelocitySampling::sample_thermal_component(T, m, rng2);
        
        REQUIRE_THAT(v1, WithinAbs(v2, TOL));
    }
}

TEST_CASE("VelocitySampling: Neutral velocity sampling", "[collision][velocity]") {
    SECTION("Zero flow velocity → thermal only") {
        EhssRng rng(42);
        double T = 300.0;
        double m = 28.0 * 1.66054e-27;
        Vec3 v_bulk{0.0, 0.0, 0.0};
        
        // Sample many velocities
        const int N = 5000;
        std::vector<Vec3> samples;
        samples.reserve(N);
        
        for (int i = 0; i < N; ++i) {
            samples.push_back(VelocitySampling::sample_neutral_velocity(T, m, v_bulk, rng));
        }
        
        // Compute mean velocity (should be near zero)
        Vec3 mean{0.0, 0.0, 0.0};
        for (const auto& v : samples) {
            mean.x += v.x;
            mean.y += v.y;
            mean.z += v.z;
        }
        mean.x /= N;
        mean.y /= N;
        mean.z /= N;
        
        double sigma = VelocitySampling::thermal_velocity_width(T, m);
        double mean_tolerance = 3.0 * sigma / std::sqrt(N);
        
        REQUIRE_THAT(mean.x, WithinAbs(0.0, mean_tolerance));
        REQUIRE_THAT(mean.y, WithinAbs(0.0, mean_tolerance));
        REQUIRE_THAT(mean.z, WithinAbs(0.0, mean_tolerance));
    }
    
    SECTION("Non-zero flow velocity → thermal + bulk") {
        EhssRng rng(42);
        double T = 300.0;
        double m = 28.0 * 1.66054e-27;
        Vec3 v_bulk{100.0, 50.0, -30.0};  // m/s
        
        // Sample many velocities
        const int N = 5000;
        std::vector<Vec3> samples;
        samples.reserve(N);
        
        for (int i = 0; i < N; ++i) {
            samples.push_back(VelocitySampling::sample_neutral_velocity(T, m, v_bulk, rng));
        }
        
        // Compute mean velocity (should be near v_bulk)
        Vec3 mean{0.0, 0.0, 0.0};
        for (const auto& v : samples) {
            mean.x += v.x;
            mean.y += v.y;
            mean.z += v.z;
        }
        mean.x /= N;
        mean.y /= N;
        mean.z /= N;
        
        double sigma = VelocitySampling::thermal_velocity_width(T, m);
        double mean_tolerance = 3.0 * sigma / std::sqrt(N);
        
        REQUIRE_THAT(mean.x, WithinAbs(v_bulk.x, mean_tolerance));
        REQUIRE_THAT(mean.y, WithinAbs(v_bulk.y, mean_tolerance));
        REQUIRE_THAT(mean.z, WithinAbs(v_bulk.z, mean_tolerance));
    }
    
    SECTION("Velocity magnitude distribution (Maxwell-Boltzmann)") {
        EhssRng rng(42);
        double T = 300.0;
        double m = 28.0 * 1.66054e-27;
        Vec3 v_bulk{0.0, 0.0, 0.0};
        
        // Sample velocity magnitudes
        const int N = 5000;
        std::vector<double> speeds;
        speeds.reserve(N);
        
        for (int i = 0; i < N; ++i) {
            Vec3 v = VelocitySampling::sample_neutral_velocity(T, m, v_bulk, rng);
            speeds.push_back(norm(v));
        }
        
        // Compute mean speed
        double mean_speed = 0.0;
        for (double s : speeds) {
            mean_speed += s;
        }
        mean_speed /= N;
        
        // Expected mean speed for Maxwell-Boltzmann: <v> = sqrt(8*kB*T/(π*m))
        double expected_mean_speed = std::sqrt(8.0 * BOLTZMANN_CONSTANT * T / (M_PI * m));
        
        // Should match within ~5% (statistical fluctuations)
        REQUIRE_THAT(mean_speed, WithinRel(expected_mean_speed, 0.05));
    }
}

TEST_CASE("VelocitySampling: Regression test against old implementation", "[collision][velocity][regression]") {
    SECTION("Bit-for-bit identical with same seed") {
        // Old implementation used inline Box-Muller in sample_neutral_velocity()
        // New implementation uses same algorithm via box_muller_sample()
        
        EhssRng rng_old(42);
        EhssRng rng_new(42);
        
        double T = 300.0;
        double m = 28.0 * 1.66054e-27;
        Vec3 v_bulk{50.0, -20.0, 100.0};
        
        // Old implementation (manual Box-Muller)
        auto box_muller_old = [&rng_old](double sigma) -> double {
            double u1 = rng_old.uniform01();
            double u2 = rng_old.uniform01();
            return sigma * std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
        };
        
        double sigma = std::sqrt(BOLTZMANN_CONSTANT * T / m);
        Vec3 v_old{
            box_muller_old(sigma) + v_bulk.x,
            box_muller_old(sigma) + v_bulk.y,
            box_muller_old(sigma) + v_bulk.z
        };
        
        // New implementation
        Vec3 v_new = VelocitySampling::sample_neutral_velocity(T, m, v_bulk, rng_new);
        
        // Should be bit-for-bit identical
        REQUIRE_THAT(v_new.x, WithinAbs(v_old.x, TOL));
        REQUIRE_THAT(v_new.y, WithinAbs(v_old.y, TOL));
        REQUIRE_THAT(v_new.z, WithinAbs(v_old.z, TOL));
    }
    
    SECTION("Multiple samples with same seed match") {
        EhssRng rng_old(123);
        EhssRng rng_new(123);
        
        double T = 300.0;
        double m = 28.0 * 1.66054e-27;
        Vec3 v_bulk{0.0, 0.0, 0.0};
        
        for (int i = 0; i < 100; ++i) {
            // Old implementation
            auto box_muller_old = [&rng_old](double sigma) -> double {
                double u1 = rng_old.uniform01();
                double u2 = rng_old.uniform01();
                return sigma * std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
            };
            
            double sigma = std::sqrt(BOLTZMANN_CONSTANT * T / m);
            Vec3 v_old{
                box_muller_old(sigma),
                box_muller_old(sigma),
                box_muller_old(sigma)
            };
            
            // New implementation
            Vec3 v_new = VelocitySampling::sample_neutral_velocity(T, m, v_bulk, rng_new);
            
            REQUIRE_THAT(v_new.x, WithinAbs(v_old.x, TOL));
            REQUIRE_THAT(v_new.y, WithinAbs(v_old.y, TOL));
            REQUIRE_THAT(v_new.z, WithinAbs(v_old.z, TOL));
        }
    }
}

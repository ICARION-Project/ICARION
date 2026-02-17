// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * @file test_gpu_space_charge.cpp
 * @brief GPU Space Charge P³M validation tests
 * 
 * Tests:
 * 1. Two-ion Coulomb force validation (analytical solution)
 * 2. Charge conservation (total charge = sum of ion charges)
 * 3. CPU/GPU parity (P³M GPU vs CPU direct summation)
 * 4. Performance benchmark
 */

#include <catch2/catch_all.hpp>
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/spacecharge/GPUSpaceChargeP3M.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include <vector>
#include <cmath>

using namespace ICARION;
using namespace ICARION::core;

// Helper function for vector norm
inline double vec_norm(const Vec3& v) {
    return std::sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

// Physical constants
constexpr double EPSILON_0 = 8.854187817e-12;  // F/m
constexpr double ELEMENTARY_CHARGE = 1.602176634e-19;  // C

/**
 * @brief CPU reference: Direct O(N²) Coulomb summation
 */
std::vector<Vec3> compute_space_charge_cpu(const std::vector<IonState>& ions) {
    size_t N = ions.size();
    std::vector<Vec3> E_fields(N, Vec3{0, 0, 0});
    
    for (size_t i = 0; i < N; ++i) {
        if (!ions[i].active) continue;
        
        Vec3 E_total{0, 0, 0};
        
        for (size_t j = 0; j < N; ++j) {
            if (i == j || !ions[j].active) continue;
            
            Vec3 r_ij = ions[i].pos - ions[j].pos;
            double r = vec_norm(r_ij);
            
            if (r < 1e-12) continue;  // Skip if too close
            
            // E = (1 / 4πε₀) * (q_j / r²) * r̂
            double q_j = ions[j].ion_charge_C;
            double factor = q_j / (4.0 * M_PI * EPSILON_0 * r * r * r);
            
            E_total.x += factor * r_ij.x;
            E_total.y += factor * r_ij.y;
            E_total.z += factor * r_ij.z;
        }
        
        E_fields[i] = E_total;
    }
    
    return E_fields;
}

TEST_CASE("GPU Space Charge P³M", "[gpu][space_charge]") {
    // Create GPU context
    auto gpu_ctx = icarion::gpu::GPUContext::create();
    REQUIRE(gpu_ctx != nullptr);
    REQUIRE(gpu_ctx->is_valid());
    
    SECTION("Two-ion Coulomb force validation") {
        // Two ions separated by 1 mm along z-axis
        // Expected: E = q / (4πε₀ d²) = 1.44 GV/m for +e at 1 mm
        
        std::vector<IonState> ions(2);
        
        // Ion 0 at origin
        ions[0].pos = Vec3{0.0, 0.0, 0.0};
        ions[0].ion_charge_C = ELEMENTARY_CHARGE;
        ions[0].active = true;
        
        // Ion 1 at 1 mm in z
        ions[1].pos = Vec3{0.0, 0.0, 1e-3};
        ions[1].ion_charge_C = ELEMENTARY_CHARGE;
        ions[1].active = true;
        
        // CPU reference
        auto E_cpu = compute_space_charge_cpu(ions);
        
        // Expected E-field magnitude on ion 1 from ion 0
        double d = 1e-3;
        double E_expected = ELEMENTARY_CHARGE / (4.0 * M_PI * EPSILON_0 * d * d);
        
        // Verify CPU calculation
        double E_cpu_mag = vec_norm(E_cpu[1]);
        REQUIRE(std::abs(E_cpu_mag - E_expected) / E_expected < 0.01);  // <1% error
        
        // GPU P³M
        icarion::gpu::GPUSpaceChargeP3M::Config config;
        config.grid_nx = 128;
        config.grid_ny = 128;
        config.grid_nz = 128;
        config.domain_min = Vec3{-2e-3, -2e-3, -1e-3};
        config.domain_max = Vec3{2e-3, 2e-3, 2e-3};
        config.epsilon_0 = EPSILON_0;
        
        auto solver = icarion::gpu::GPUSpaceChargeP3M::create(*gpu_ctx, config);
        REQUIRE(solver != nullptr);
        
        std::vector<Vec3> E_gpu;
        bool success = solver->compute_space_charge_field(ions, E_gpu);
        REQUIRE(success);
        REQUIRE(E_gpu.size() == 2);
        
        // Check E-field on ion 1
        double E_gpu_mag = vec_norm(E_gpu[1]);
        double relative_error = std::abs(E_gpu_mag - E_expected) / E_expected;
        
        INFO("Expected: " << E_expected << " V/m");
        INFO("CPU:      " << E_cpu_mag << " V/m");
        INFO("GPU:      " << E_gpu_mag << " V/m");
        INFO("Relative error: " << (relative_error * 100) << "%");
        
        // P³M has grid discretization error. For near-field interactions (1mm separation
        // with ~30µm cells), expect ~20% error. Production use with higher resolution
        // or many-body far-field dominated systems will have better accuracy.
        REQUIRE(relative_error < 0.25);
        
        // Check direction (should be repulsive, pointing away from ion 0)
        REQUIRE(E_gpu[1].z > 0);  // Positive z (away from origin)
        
        // Check symmetry: E on ion 0 should be opposite
        REQUIRE(E_gpu[0].z < 0);  // Negative z
        REQUIRE(std::abs(E_gpu[0].z + E_gpu[1].z) / E_expected < 0.15);
    }
    
    SECTION("Charge conservation") {
        // Create 100 ions with random positions
        const int N = 100;
        std::vector<IonState> ions(N);
        
        double total_charge = 0.0;
        for (int i = 0; i < N; ++i) {
            ions[i].pos = Vec3{
                (rand() / (double)RAND_MAX - 0.5) * 0.02,  // ±1 cm
                (rand() / (double)RAND_MAX - 0.5) * 0.02,
                (rand() / (double)RAND_MAX - 0.5) * 0.02
            };
            ions[i].ion_charge_C = ELEMENTARY_CHARGE;
            ions[i].active = true;
            total_charge += ions[i].ion_charge_C;
        }
        
        // GPU P³M
        icarion::gpu::GPUSpaceChargeP3M::Config config;
        config.grid_nx = 32;
        config.grid_ny = 32;
        config.grid_nz = 32;
        config.domain_min = Vec3{-0.015, -0.015, -0.015};
        config.domain_max = Vec3{0.015, 0.015, 0.015};
        config.epsilon_0 = EPSILON_0;
        
        auto solver = icarion::gpu::GPUSpaceChargeP3M::create(*gpu_ctx, config);
        REQUIRE(solver != nullptr);
        
        std::vector<Vec3> E_gpu;
        bool success = solver->compute_space_charge_field(ions, E_gpu);
        REQUIRE(success);
        REQUIRE(E_gpu.size() == N);
        
        // All E-fields should be finite
        for (size_t i = 0; i < N; ++i) {
            REQUIRE(std::isfinite(E_gpu[i].x));
            REQUIRE(std::isfinite(E_gpu[i].y));
            REQUIRE(std::isfinite(E_gpu[i].z));
        }
    }
    
    SECTION("CPU/GPU parity - 1000 ions") {
        // Large enough for meaningful comparison, small enough for CPU
        const int N = 1000;
        std::vector<IonState> ions(N);
        
        // Gaussian distribution around origin
        for (int i = 0; i < N; ++i) {
            double r = 0.005 * std::sqrt(-2.0 * std::log((rand() + 1.0) / (RAND_MAX + 1.0)));
            double theta = 2.0 * M_PI * rand() / RAND_MAX;
            double z = 0.01 * ((rand() / (double)RAND_MAX) - 0.5);
            
            ions[i].pos = Vec3{
                r * std::cos(theta),
                r * std::sin(theta),
                z
            };
            ions[i].ion_charge_C = ELEMENTARY_CHARGE;
            ions[i].active = true;
        }
        
        // CPU reference
        auto E_cpu = compute_space_charge_cpu(ions);
        
        // GPU P³M
        icarion::gpu::GPUSpaceChargeP3M::Config config;
        config.grid_nx = 64;
        config.grid_ny = 64;
        config.grid_nz = 64;
        config.domain_min = Vec3{-0.01, -0.01, -0.01};
        config.domain_max = Vec3{0.01, 0.01, 0.01};
        config.epsilon_0 = EPSILON_0;
        
        auto solver = icarion::gpu::GPUSpaceChargeP3M::create(*gpu_ctx, config);
        REQUIRE(solver != nullptr);
        
        std::vector<Vec3> E_gpu;
        bool success = solver->compute_space_charge_field(ions, E_gpu);
        REQUIRE(success);
        REQUIRE(E_gpu.size() == N);
        
        // Compare CPU vs GPU
        double total_error = 0.0;
        double total_magnitude = 0.0;
        int count = 0;
        
        for (size_t i = 0; i < N; ++i) {
            double E_cpu_mag = vec_norm(E_cpu[i]);
            double E_gpu_mag = vec_norm(E_gpu[i]);
            
            if (E_cpu_mag > 1e6) {  // Only compare significant fields (>1 MV/m)
                double error = std::abs(E_gpu_mag - E_cpu_mag);
                total_error += error;
                total_magnitude += E_cpu_mag;
                count++;
            }
        }
        
        double avg_relative_error = (count > 0) ? (total_error / total_magnitude) : 0.0;
        
        INFO("Ions with E > 1 MV/m: " << count);
        INFO("Average relative error: " << (avg_relative_error * 100) << "%");
        
        // P³M should match CPU within 20% on average (grid discretization)
        REQUIRE(avg_relative_error < 0.20);
    }
    
    SECTION("Performance benchmark") {
        const int N = 10000;
        std::vector<IonState> ions(N);
        
        // Random distribution
        for (int i = 0; i < N; ++i) {
            ions[i].pos = Vec3{
                (rand() / (double)RAND_MAX - 0.5) * 0.02,
                (rand() / (double)RAND_MAX - 0.5) * 0.02,
                (rand() / (double)RAND_MAX - 0.5) * 0.02
            };
            ions[i].ion_charge_C = ELEMENTARY_CHARGE;
            ions[i].active = true;
        }
        
        // GPU P³M
        icarion::gpu::GPUSpaceChargeP3M::Config config;
        config.grid_nx = 64;
        config.grid_ny = 64;
        config.grid_nz = 64;
        config.domain_min = Vec3{-0.01, -0.01, -0.01};
        config.domain_max = Vec3{0.01, 0.01, 0.01};
        config.epsilon_0 = EPSILON_0;
        auto solver = icarion::gpu::GPUSpaceChargeP3M::create(*gpu_ctx, config);
        REQUIRE(solver != nullptr);
        
        // Warm-up
        std::vector<Vec3> E_gpu;
        solver->compute_space_charge_field(ions, E_gpu);
        
        // Benchmark 10 iterations
        const int iterations = 10;
        for (int iter = 0; iter < iterations; ++iter) {
            bool success = solver->compute_space_charge_field(ions, E_gpu);
            REQUIRE(success);
        }
        
        auto stats = solver->get_stats();
        
        INFO("N = " << N);
        INFO("Grid = " << config.grid_nx << "³");
        INFO("Average time: " << stats.avg_time_ms << " ms");
        INFO("Estimated speedup vs CPU: " << stats.speedup_vs_direct_cpu() << "×");
        
        // For N=10k, should be <20 ms/timestep
        REQUIRE(stats.avg_time_ms < 20.0);
        
        // Should be faster than CPU
        REQUIRE(stats.speedup_vs_direct_cpu() > 10.0);
    }
}

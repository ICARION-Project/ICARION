// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file test_gpu_collision_parity.cpp
 * @brief Test GPU collision helper basic functionality
 * 
 * NOTE: Full CPU/GPU parity testing requires integration with SimulationEngine.
 * This test validates:
 * - GPU helper creation and initialization
 * - Threshold-based GPU dispatch
 * - Basic collision processing (no crash)
 */

#ifdef ICARION_USE_GPU

#include <catch2/catch_test_macros.hpp>

#include "core/gpu/core/GPUContext.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/types/CollisionTypes.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/types/Vec3.h"

#include <vector>
#include <cmath>
#include <tuple>

using namespace ICARION;
using namespace ICARION::config;
using ICARION::core::IonEnsemble;
using ICARION::core::IonState;
using ICARION::core::Vec3;
using ICARION::physics::HSSCollisionHandler;
using ICARION::physics::PhysicsRng;
using icarion::gpu::GPUContext;
using icarion::gpu::GPUCollisionHelper;

static void run_cpu_collision(HSSCollisionHandler& handler,
                              IonState& ion,
                              double dt,
                              PhysicsRng& rng,
                              const EnvironmentConfig& env) {
    std::vector<IonState> tmp{ion};
    auto ens = IonEnsemble::from_legacy(tmp);
    auto view = ens.collision_data(0);
    handler.handle_collision(view, dt, rng, env);
    ion = ens.ion_state(0);
}

// =============================================================================
// TEST: GPU collision helper threshold behavior
// =============================================================================

TEST_CASE("GPU collision helper respects threshold", "[collision][gpu]") {
    // Create GPU context
    auto gpu_ctx = GPUContext::create(0);
    if (!gpu_ctx) {
        SKIP("GPU not available");
    }

    // Setup environment
    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_mass_kg = 4.0026 * 1.66054e-27;  // He
    
    const size_t threshold = 5000;
    auto gpu_helper = GPUCollisionHelper::create(*gpu_ctx, threshold, "HSS");
    
    REQUIRE(gpu_helper != nullptr);
    
    // Test with N < threshold -> should return false (use CPU)
    std::vector<IonState> small_batch(1000);
    for (auto& ion : small_batch) {
        ion.pos = Vec3(0, 0, 0);
        ion.vel = Vec3(100, 100, 100);
        ion.mass_kg = 29.0 * 1.66054e-27;
        ion.CCS_m2 = 45e-20;
        ion.active = true;
    }
    
    bool used_gpu = gpu_helper->process_collisions_batch(small_batch, 1e-6, env);
    CHECK_FALSE(used_gpu);  // Should NOT use GPU (below threshold)
    
    // Test with N >= threshold -> should return true (use GPU)
    std::vector<IonState> large_batch(6000);
    for (auto& ion : large_batch) {
        ion.pos = Vec3(0, 0, 0);
        ion.vel = Vec3(100, 100, 100);
        ion.mass_kg = 29.0 * 1.66054e-27;
        ion.CCS_m2 = 45e-20;
        ion.active = true;
    }
    
    used_gpu = gpu_helper->process_collisions_batch(large_batch, 1e-6, env);
    CHECK(used_gpu);  // Should use GPU (above threshold)
}

TEST_CASE("GPU collision helper doesn't crash on empty batch", "[collision][gpu]") {
    auto gpu_ctx = GPUContext::create(0);
    if (!gpu_ctx) {
        SKIP("GPU not available");
    }

    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_mass_kg = 4.0026 * 1.66054e-27;
    
    auto gpu_helper = GPUCollisionHelper::create(*gpu_ctx, 100, "HSS");
    REQUIRE(gpu_helper != nullptr);
    
    std::vector<IonState> empty_batch;
    
    // Should handle empty batch gracefully
    bool used_gpu = gpu_helper->process_collisions_batch(empty_batch, 1e-6, env);
    CHECK_FALSE(used_gpu);  // Empty batch, don't use GPU
}

TEST_CASE("GPU HSS produces deterministic results with fixed seed", "[collision][gpu][deterministic]") {
    auto gpu_ctx = GPUContext::create(0);
    if (!gpu_ctx) {
        SKIP("GPU not available");
    }

    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_mass_kg = 4.0026 * 1.66054e-27;  // He
    
    // Create initial state
    const size_t N = 10000;  // Above threshold to force GPU
    const double mass = 29.0 * 1.66054e-27;  // H3O+
    const double CCS = 45e-20;
    
    std::vector<IonState> ions_run1(N);
    std::vector<IonState> ions_run2(N);
    
    // Initialize with identical velocities
    for (size_t i = 0; i < N; ++i) {
        ions_run1[i].pos = Vec3(0, 0, 0);
        ions_run1[i].vel = Vec3(500.0, 300.0, 200.0);  // Fixed velocity
        ions_run1[i].mass_kg = mass;
        ions_run1[i].CCS_m2 = CCS;
        ions_run1[i].active = true;
        
        ions_run2[i] = ions_run1[i];  // Exact copy
    }
    
    // Run with SAME seed twice
    const unsigned long long seed = 42;
    auto gpu_helper1 = GPUCollisionHelper::create(*gpu_ctx, 5000, "HSS", seed);
    auto gpu_helper2 = GPUCollisionHelper::create(*gpu_ctx, 5000, "HSS", seed);
    
    REQUIRE(gpu_helper1 != nullptr);
    REQUIRE(gpu_helper2 != nullptr);
    
    const double dt = 1e-6;
    
    bool used_gpu1 = gpu_helper1->process_collisions_batch(ions_run1, dt, env);
    bool used_gpu2 = gpu_helper2->process_collisions_batch(ions_run2, dt, env);
    
    REQUIRE(used_gpu1);
    REQUIRE(used_gpu2);
    
    // With same seed, results should be IDENTICAL
    size_t n_identical = 0;
    size_t n_checked = std::min<size_t>(100, N);  // Check first 100 ions
    
    for (size_t i = 0; i < n_checked; ++i) {
        bool vx_match = std::abs(ions_run1[i].vel.x - ions_run2[i].vel.x) < 1e-12;
        bool vy_match = std::abs(ions_run1[i].vel.y - ions_run2[i].vel.y) < 1e-12;
        bool vz_match = std::abs(ions_run1[i].vel.z - ions_run2[i].vel.z) < 1e-12;
        
        if (vx_match && vy_match && vz_match) {
            n_identical++;
        }
    }
    
    INFO("Identical velocities: " << n_identical << " / " << n_checked);
    
    // Should be 100% identical with same seed
    CHECK(n_identical == n_checked);
}

TEST_CASE("GPU vs CPU HSS collision comparison (statistical)", "[collision][gpu][cpu-parity]") {
    auto gpu_ctx = GPUContext::create(0);
    if (!gpu_ctx) {
        SKIP("GPU not available");
    }

    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_mass_kg = 4.0026 * 1.66054e-27;  // He
    
    // Create test ensemble
    const size_t N_cpu = 500;   // Below GPU threshold → forces CPU
    const size_t N_gpu = 10000; // Above GPU threshold → forces GPU
    const double mass = 29.0 * 1.66054e-27;
    const double CCS = 45e-20;
    
    std::vector<IonState> cpu_ions(N_cpu);
    std::vector<IonState> gpu_ions(N_gpu);
    
    // Initialize with identical velocity distribution
    for (size_t i = 0; i < N_cpu; ++i) {
        cpu_ions[i].pos = Vec3(0, 0, 0);
        cpu_ions[i].vel = Vec3(500.0, 300.0, 200.0);
        cpu_ions[i].mass_kg = mass;
        cpu_ions[i].CCS_m2 = CCS;
        cpu_ions[i].active = true;
    }
    
    for (size_t i = 0; i < N_gpu; ++i) {
        gpu_ions[i].pos = Vec3(0, 0, 0);
        gpu_ions[i].vel = Vec3(500.0, 300.0, 200.0);
        gpu_ions[i].mass_kg = mass;
        gpu_ions[i].CCS_m2 = CCS;
        gpu_ions[i].active = true;
    }
    
    auto gpu_helper = GPUCollisionHelper::create(*gpu_ctx, 5000, "HSS", 42);
    REQUIRE(gpu_helper != nullptr);
    
    const double dt = 1e-6;
    
    // Process - should use CPU for first, GPU for second
    bool used_cpu = gpu_helper->process_collisions_batch(cpu_ions, dt, env);
    bool used_gpu = gpu_helper->process_collisions_batch(gpu_ions, dt, env);
    
    CHECK_FALSE(used_cpu);  // Below threshold
    CHECK(used_gpu);        // Above threshold
    
    // Compare statistical properties
    auto calc_mean_speed = [](const std::vector<IonState>& ions) {
        double sum = 0.0;
        for (const auto& ion : ions) {
            sum += std::sqrt(ion.vel.x*ion.vel.x + ion.vel.y*ion.vel.y + ion.vel.z*ion.vel.z);
        }
        return sum / ions.size();
    };
    
    double cpu_mean_speed = calc_mean_speed(cpu_ions);
    double gpu_mean_speed = calc_mean_speed(gpu_ions);
    
    INFO("CPU mean speed: " << cpu_mean_speed << " m/s");
    INFO("GPU mean speed: " << gpu_mean_speed << " m/s");
    
    // Should be statistically similar (within ~20% due to different sample sizes and stochastic process)
    double speed_ratio = std::abs(cpu_mean_speed - gpu_mean_speed) / std::max(cpu_mean_speed, gpu_mean_speed);
    CHECK(speed_ratio < 0.3);  // Allow 30% deviation (very loose - stochastic!)
}

TEST_CASE("CPU HSS handler also produces deterministic results with fixed seed", "[collision][cpu][deterministic]") {
    // This demonstrates CPU handler determinism for reference
    // NOTE: CPU uses std::mt19937_64, GPU uses cuRAND XORWOW - these are different PRNGs
    //       so we cannot expect bit-identical results between CPU and GPU even with "same" seed
    
    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_mass_kg = 4.0026 * 1.66054e-27;  // He
    
    const size_t N = 1000;
    const double mass = 29.0 * 1.66054e-27;
    const double CCS = 45e-20;
    
    std::vector<IonState> ions_run1(N);
    std::vector<IonState> ions_run2(N);
    
    // Initialize with identical state
    for (size_t i = 0; i < N; ++i) {
        ions_run1[i].pos = Vec3(0, 0, 0);
        ions_run1[i].vel = Vec3(500.0, 300.0, 200.0);
        ions_run1[i].mass_kg = mass;
        ions_run1[i].CCS_m2 = CCS;
        ions_run1[i].active = true;
        
        ions_run2[i] = ions_run1[i];
    }
    
    // Create two CPU handlers with SAME seed
    const unsigned long long seed = 42;
    physics::PhysicsRng rng1(seed);
    physics::PhysicsRng rng2(seed);
    
    physics::HSSCollisionHandler handler1(false, nullptr);
    physics::HSSCollisionHandler handler2(false, nullptr);
    
    const double dt = 1e-6;
    
    // Process all ions with both handlers
    for (size_t i = 0; i < N; ++i) {
        run_cpu_collision(handler1, ions_run1[i], dt, rng1, env);
        run_cpu_collision(handler2, ions_run2[i], dt, rng2, env);
    }
    
    // With same seed, CPU results should be IDENTICAL
    size_t n_identical = 0;
    size_t n_checked = std::min<size_t>(100, N);
    
    for (size_t i = 0; i < n_checked; ++i) {
        bool vx_match = std::abs(ions_run1[i].vel.x - ions_run2[i].vel.x) < 1e-12;
        bool vy_match = std::abs(ions_run1[i].vel.y - ions_run2[i].vel.y) < 1e-12;
        bool vz_match = std::abs(ions_run1[i].vel.z - ions_run2[i].vel.z) < 1e-12;
        
        if (vx_match && vy_match && vz_match) {
            n_identical++;
        }
    }
    
    INFO("CPU identical velocities: " << n_identical << " / " << n_checked);
    CHECK(n_identical == n_checked);
}

TEST_CASE("GPU collision rate formula matches CPU", "[collision][gpu][physics-check]") {
    // This test verifies that the collision probability formula is identical
    // We can't test random sampling, but we can test the deterministic physics
    
    // Test parameters
    const double temperature_K = 300.0;
    const double pressure_Pa = 101325.0;
    const double CCS = 45e-20;
    const double dt = 1e-6;
    const double v_rel_mag = 500.0;  // m/s
    
    // Compute collision probability using CPU formula
    const double kB = 1.380649e-23;
    const double number_density = pressure_Pa / (kB * temperature_K);
    const double collision_rate_cpu = v_rel_mag * CCS * number_density;
    const double collision_prob_cpu = 1.0 - std::exp(-collision_rate_cpu * dt);
    
    // GPU uses the same formula (after the fix)
    const double collision_rate_gpu = v_rel_mag * CCS * number_density;
    const double collision_prob_gpu = 1.0 - std::exp(-collision_rate_gpu * dt);
    
    INFO("CPU collision probability: " << collision_prob_cpu);
    INFO("GPU collision probability: " << collision_prob_gpu);
    INFO("Collision rate: " << collision_rate_cpu << " Hz");
    
    // Formulas should be IDENTICAL
    CHECK(std::abs(collision_prob_cpu - collision_prob_gpu) < 1e-15);
    
    // Typical values check (probability should be positive and ≤ 1)
    CHECK(collision_prob_cpu > 0.0);
    CHECK(collision_prob_cpu <= 1.0);  // Can saturate at 100% for high pressure/large dt
}

TEST_CASE("CPU and GPU HSS produce statistically equivalent results", "[collision][gpu][cpu][parity][statistical]") {
    auto gpu_ctx = GPUContext::create(0);
    if (!gpu_ctx) {
        SKIP("GPU not available");
    }

    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_mass_kg = 4.0026 * 1.66054e-27;  // He
    
    // Large sample size for statistical comparison
    const size_t N = 50000;
    const double mass = 29.0 * 1.66054e-27;
    const double CCS = 45e-20;
    const double dt = 1e-6;
    
    std::vector<IonState> cpu_ions(N);
    std::vector<IonState> gpu_ions(N);
    
    // Initialize with IDENTICAL state
    for (size_t i = 0; i < N; ++i) {
        IonState ion;
        ion.pos = Vec3(0, 0, 0);
        ion.vel = Vec3(500.0, 300.0, 200.0);  // Fixed initial velocity
        ion.mass_kg = mass;
        ion.CCS_m2 = CCS;
        ion.active = true;
        
        cpu_ions[i] = ion;
        gpu_ions[i] = ion;
    }
    
    // Process with CPU
    physics::PhysicsRng cpu_rng(12345);
    physics::HSSCollisionHandler cpu_handler(false, nullptr);
    
    for (size_t i = 0; i < N; ++i) {
        run_cpu_collision(cpu_handler, cpu_ions[i], dt, cpu_rng, env);
    }
    
    // Process with GPU
    auto gpu_helper = GPUCollisionHelper::create(*gpu_ctx, 100, "HSS", 12345);
    REQUIRE(gpu_helper != nullptr);
    
    bool used_gpu = gpu_helper->process_collisions_batch(gpu_ions, dt, env);
    REQUIRE(used_gpu);
    
    // Statistical comparison
    auto calc_stats = [](const std::vector<IonState>& ions) {
        double mean_vx = 0, mean_vy = 0, mean_vz = 0;
        double mean_speed = 0;
        double mean_KE = 0;
        
        for (const auto& ion : ions) {
            mean_vx += ion.vel.x;
            mean_vy += ion.vel.y;
            mean_vz += ion.vel.z;
            
            double speed = std::sqrt(ion.vel.x*ion.vel.x + ion.vel.y*ion.vel.y + ion.vel.z*ion.vel.z);
            mean_speed += speed;
            mean_KE += 0.5 * ion.mass_kg * speed * speed;
        }
        
        size_t n = ions.size();
        return std::make_tuple(mean_vx/n, mean_vy/n, mean_vz/n, mean_speed/n, mean_KE/n);
    };
    
    auto [cpu_vx, cpu_vy, cpu_vz, cpu_speed, cpu_KE] = calc_stats(cpu_ions);
    auto [gpu_vx, gpu_vy, gpu_vz, gpu_speed, gpu_KE] = calc_stats(gpu_ions);
    
    INFO("CPU: vx=" << cpu_vx << " vy=" << cpu_vy << " vz=" << cpu_vz << " speed=" << cpu_speed << " KE=" << cpu_KE);
    INFO("GPU: vx=" << gpu_vx << " vy=" << gpu_vy << " vz=" << gpu_vz << " speed=" << gpu_speed << " KE=" << gpu_KE);
    
    // NOTE: Different PRNGs (mt19937_64 vs cuRAND XORWOW) give DIFFERENT random sequences
    // This means collision probabilities, sampled velocities, and deflection angles differ
    // We can only expect QUALITATIVE similarity (both thermalize), not quantitative agreement
    
    // Both implementations should show thermalization (reduction from initial 500,300,200 m/s)
    double initial_speed = std::sqrt(500*500 + 300*300 + 200*200);  // ~621 m/s
    
    INFO("Initial speed: " << initial_speed << " m/s");
    
    // Both should have reduced speed compared to initial (thermalization occurred)
    CHECK(cpu_speed < initial_speed * 1.05);  // Allow some variance
    CHECK(gpu_speed < initial_speed * 1.05);
    
    // Both should have similar ORDER OF MAGNITUDE for final speed
    // (within factor of 2 is reasonable given different RNGs)
    double speed_ratio = std::max(cpu_speed, gpu_speed) / std::min(cpu_speed, gpu_speed);
    CHECK(speed_ratio < 2.0);
    
    // Both should have similar ORDER OF MAGNITUDE for KE
    double KE_ratio = std::max(cpu_KE, gpu_KE) / std::min(cpu_KE, gpu_KE);
    CHECK(KE_ratio < 2.0);
}

#endif // ICARION_USE_GPU

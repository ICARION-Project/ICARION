// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file test_rk45_boris_parity.cpp
 * @brief CPU/GPU parity tests for RK45 and Boris integrators
 * 
 * Validates that GPU implementations produce identical results to CPU
 * for the same inputs and parameters.
 */

#ifdef ICARION_USE_GPU

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/gpu/core/GPUContext.h"
#include "core/gpu/core/GPUIntegrationHelper.h"
#include "core/integrator/strategies/RK45Strategy.h"
#include "core/integrator/strategies/BorisStrategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/types/Vec3.h"

#include <vector>
#include <cmath>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;
using Catch::Matchers::WithinAbs;
using ICARION::core::IonEnsemble;
using ICARION::core::IonState;
using ICARION::core::Vec3;

// =============================================================================
// HELPER: Create identical ion sets
// =============================================================================
struct IonBatch {
    std::vector<IonState> cpu_ions;
    std::vector<IonState> gpu_ions;
};

IonBatch create_identical_ions(size_t count, Vec3 pos0, Vec3 vel0, double mass = 1.0e-26) {
    IonBatch batch;
    batch.cpu_ions.resize(count);
    batch.gpu_ions.resize(count);
    
    for (size_t i = 0; i < count; ++i) {
        IonState ion;
        ion.pos = pos0;
        ion.vel = vel0;
        ion.mass_kg = mass;
        ion.ion_charge_C = 1.602e-19;
        ion.CCS_m2 = 45e-20;
        ion.active = true;
        ion.born = true;
        ion.t = 0.0;
        
        batch.cpu_ions[i] = ion;
        batch.gpu_ions[i] = ion;
    }
    
    return batch;
}

// Helper to run CPU RK45 for entire ensemble
static void rk45_cpu_step_batch(RK45Strategy& strategy,
                                IonEnsemble& ensemble,
                                double t,
                                double dt,
                                const ForceRegistry& forces) {
    for (size_t idx = 0; idx < ensemble.size(); ++idx) {
        strategy.step(ensemble, idx, t, dt, forces);
    }
}

static void boris_cpu_step_batch(BorisStrategy& strategy,
                                 IonEnsemble& ensemble,
                                 double t,
                                 double dt,
                                 const ForceRegistry& forces) {
    for (size_t idx = 0; idx < ensemble.size(); ++idx) {
        strategy.step(ensemble, idx, t, dt, forces);
    }
}

// =============================================================================
// HELPER: Compare ion states
// =============================================================================
void compare_ions(const IonState& cpu, const IonState& gpu, 
                  double pos_tol = 1e-10, double vel_tol = 1e-10) {
    REQUIRE_THAT(gpu.pos.x, WithinAbs(cpu.pos.x, pos_tol));
    REQUIRE_THAT(gpu.pos.y, WithinAbs(cpu.pos.y, pos_tol));
    REQUIRE_THAT(gpu.pos.z, WithinAbs(cpu.pos.z, pos_tol));
    
    REQUIRE_THAT(gpu.vel.x, WithinAbs(cpu.vel.x, vel_tol));
    REQUIRE_THAT(gpu.vel.y, WithinAbs(cpu.vel.y, vel_tol));
    REQUIRE_THAT(gpu.vel.z, WithinAbs(cpu.vel.z, vel_tol));
}

// =============================================================================
// RK45 PARITY TESTS
// =============================================================================

TEST_CASE("RK45: CPU/GPU parity - Free particle", "[rk45][parity][gpu]") {
    // Setup
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    REQUIRE(gpu_ctx->is_valid());
    
    auto gpu_helper = icarion::gpu::GPUIntegrationHelper::create(*gpu_ctx, 100);
    REQUIRE(gpu_helper != nullptr);
    
    auto rk45_cpu = std::make_shared<RK45Strategy>();
    
    // Create dummy domain and force registry (no forces)
    DomainConfig domain;
    domain.solver = SolverType::RK45;
    ForceRegistry forces(domain);  // Pass domain in constructor
    
    // Test parameters
    const size_t N = 1000;
    const double dt = 1e-6;  // 1 μs
    const Vec3 pos0{0.001, 0.002, 0.003};  // mm scale
    const Vec3 vel0{1000.0, 500.0, 200.0};  // m/s
    
    auto batch = create_identical_ions(N, pos0, vel0);
    
    auto cpu_ensemble = IonEnsemble::from_legacy(batch.cpu_ions);
    rk45_cpu_step_batch(*rk45_cpu, cpu_ensemble, 0.0, dt, forces);
    batch.cpu_ions = cpu_ensemble.to_legacy();
    
    // GPU integration (batch)
    bool success = gpu_helper->integrate_batch_rk45(batch.gpu_ions, dt, 0.0, nullptr);
    REQUIRE(success);
    
    // Compare results (should be exact for free particle)
    INFO("Comparing " << N << " ions after 1 timestep");
    for (size_t i = 0; i < N; ++i) {
        INFO("Ion " << i);
        compare_ions(batch.cpu_ions[i], batch.gpu_ions[i], 1e-12, 1e-12);
    }
}

TEST_CASE("RK45: CPU/GPU parity - Multi-step trajectory", "[rk45][parity][gpu]") {
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    REQUIRE(gpu_ctx->is_valid());
    
    auto gpu_helper = icarion::gpu::GPUIntegrationHelper::create(*gpu_ctx, 100);
    auto rk45_cpu = std::make_shared<RK45Strategy>();
    
    DomainConfig domain;
    domain.solver = SolverType::RK45;
    ForceRegistry forces(domain);
    
    const size_t N = 500;
    const double dt = 1e-6;
    const int steps = 10;
    const Vec3 pos0{0.0, 0.0, 0.0};
    const Vec3 vel0{500.0, 300.0, 100.0};
    
    auto batch = create_identical_ions(N, pos0, vel0);
    
    auto cpu_ensemble = IonEnsemble::from_legacy(batch.cpu_ions);
    
    // Multi-step integration
    for (int step = 0; step < steps; ++step) {
        double t = step * dt;
        
        // CPU
        rk45_cpu_step_batch(*rk45_cpu, cpu_ensemble, t, dt, forces);
        
        // GPU
        bool success = gpu_helper->integrate_batch_rk45(batch.gpu_ions, dt, t, nullptr);
        REQUIRE(success);
    }
    
    batch.cpu_ions = cpu_ensemble.to_legacy();
    
    // Compare final states
    INFO("Comparing after " << steps << " timesteps");
    for (size_t i = 0; i < N; ++i) {
        INFO("Ion " << i);
        // Accumulated error over 10 steps → slightly relaxed tolerance
        compare_ions(batch.cpu_ions[i], batch.gpu_ions[i], 1e-11, 1e-11);
    }
}

TEST_CASE("RK45: CPU/GPU parity - Various initial conditions", "[rk45][parity][gpu]") {
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    REQUIRE(gpu_ctx->is_valid());
    
    auto gpu_helper = icarion::gpu::GPUIntegrationHelper::create(*gpu_ctx, 100);
    auto rk45_cpu = std::make_shared<RK45Strategy>();
    
    DomainConfig domain;
    domain.solver = SolverType::RK45;
    ForceRegistry forces(domain);
    
    const double dt = 1e-6;
    
    SECTION("Zero velocity") {
        auto batch = create_identical_ions(500, Vec3{0.001, 0.0, 0.0}, Vec3{0.0, 0.0, 0.0});
        auto cpu_ensemble = IonEnsemble::from_legacy(batch.cpu_ions);
        rk45_cpu_step_batch(*rk45_cpu, cpu_ensemble, 0.0, dt, forces);
        batch.cpu_ions = cpu_ensemble.to_legacy();
        
        bool success = gpu_helper->integrate_batch_rk45(batch.gpu_ions, dt, 0.0, nullptr);
        REQUIRE(success);
        
        for (size_t i = 0; i < batch.cpu_ions.size(); ++i) {
            compare_ions(batch.cpu_ions[i], batch.gpu_ions[i]);
        }
    }
    
    SECTION("High velocity") {
        auto batch = create_identical_ions(500, Vec3{0.0, 0.0, 0.0}, Vec3{1e6, 5e5, 2e5});
        auto cpu_ensemble = IonEnsemble::from_legacy(batch.cpu_ions);
        rk45_cpu_step_batch(*rk45_cpu, cpu_ensemble, 0.0, dt, forces);
        batch.cpu_ions = cpu_ensemble.to_legacy();
        
        bool success = gpu_helper->integrate_batch_rk45(batch.gpu_ions, dt, 0.0, nullptr);
        REQUIRE(success);
        
        for (size_t i = 0; i < batch.cpu_ions.size(); ++i) {
            compare_ions(batch.cpu_ions[i], batch.gpu_ions[i], 1e-9, 1e-9);  // Relaxed for high v
        }
    }
    
    SECTION("Various masses") {
        std::vector<double> masses = {1e-27, 1e-26, 1e-25, 5e-26};
        
        for (double mass : masses) {
            INFO("Testing mass: " << mass << " kg");
            auto batch = create_identical_ions(200, Vec3{0.0, 0.0, 0.0}, Vec3{1000.0, 0.0, 0.0}, mass);
            auto cpu_ensemble = IonEnsemble::from_legacy(batch.cpu_ions);
            rk45_cpu_step_batch(*rk45_cpu, cpu_ensemble, 0.0, dt, forces);
            batch.cpu_ions = cpu_ensemble.to_legacy();
            
            bool success = gpu_helper->integrate_batch_rk45(batch.gpu_ions, dt, 0.0, nullptr);
            REQUIRE(success);
            
            for (size_t i = 0; i < batch.cpu_ions.size(); ++i) {
                compare_ions(batch.cpu_ions[i], batch.gpu_ions[i]);
            }
        }
    }
}

// =============================================================================
// BORIS PARITY TESTS
// =============================================================================

TEST_CASE("Boris: CPU/GPU parity - Free particle", "[boris][parity][gpu]") {
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    REQUIRE(gpu_ctx->is_valid());
    
    auto gpu_helper = icarion::gpu::GPUIntegrationHelper::create(*gpu_ctx, 100);
    auto boris_cpu = std::make_shared<BorisStrategy>();
    
    DomainConfig domain;
    domain.solver = SolverType::Boris;
    ForceRegistry forces(domain);
    
    const size_t N = 1000;
    const double dt = 1e-7;  // 100 ns (smaller for Boris)
    const Vec3 pos0{0.0, 0.0, 0.0};
    const Vec3 vel0{1000.0, 500.0, 200.0};
    
    auto batch = create_identical_ions(N, pos0, vel0);
    
    auto cpu_ensemble = IonEnsemble::from_legacy(batch.cpu_ions);
    boris_cpu_step_batch(*boris_cpu, cpu_ensemble, 0.0, dt, forces);
    batch.cpu_ions = cpu_ensemble.to_legacy();
    
    // GPU integration
    bool success = gpu_helper->integrate_batch_boris(batch.gpu_ions, dt, 0.0, nullptr);
    REQUIRE(success);
    
    // Compare results
    for (size_t i = 0; i < N; ++i) {
        INFO("Ion " << i);
        compare_ions(batch.cpu_ions[i], batch.gpu_ions[i], 1e-12, 1e-12);
    }
}

TEST_CASE("Boris: CPU/GPU parity - Multi-step trajectory", "[boris][parity][gpu]") {
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    REQUIRE(gpu_ctx->is_valid());
    
    auto gpu_helper = icarion::gpu::GPUIntegrationHelper::create(*gpu_ctx, 100);
    auto boris_cpu = std::make_shared<BorisStrategy>();
    
    DomainConfig domain;
    domain.solver = SolverType::Boris;
    ForceRegistry forces(domain);
    
    const size_t N = 500;
    const double dt = 1e-7;
    const int steps = 20;
    const Vec3 pos0{0.0, 0.0, 0.0};
    const Vec3 vel0{500.0, 300.0, 100.0};
    
    auto batch = create_identical_ions(N, pos0, vel0);
    
    auto cpu_ensemble = IonEnsemble::from_legacy(batch.cpu_ions);
    
    // Multi-step integration
    for (int step = 0; step < steps; ++step) {
        double t = step * dt;
        boris_cpu_step_batch(*boris_cpu, cpu_ensemble, t, dt, forces);
        
        // GPU
        bool success = gpu_helper->integrate_batch_boris(batch.gpu_ions, dt, t, nullptr);
        REQUIRE(success);
    }
    
    batch.cpu_ions = cpu_ensemble.to_legacy();
    
    // Compare final states
    INFO("Comparing after " << steps << " timesteps");
    for (size_t i = 0; i < N; ++i) {
        INFO("Ion " << i);
        compare_ions(batch.cpu_ions[i], batch.gpu_ions[i], 1e-11, 1e-11);
    }
}

TEST_CASE("Boris: CPU/GPU parity - Energy conservation", "[boris][parity][gpu]") {
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    REQUIRE(gpu_ctx->is_valid());
    
    auto gpu_helper = icarion::gpu::GPUIntegrationHelper::create(*gpu_ctx, 100);
    auto boris_cpu = std::make_shared<BorisStrategy>();
    
    DomainConfig domain;
    domain.solver = SolverType::Boris;
    ForceRegistry forces(domain);
    
    const size_t N = 500;
    const double dt = 1e-7;
    const int steps = 100;
    const Vec3 pos0{0.0, 0.0, 0.0};
    const Vec3 vel0{1000.0, 1000.0, 1000.0};
    const double mass = 29.0 * 1.66054e-27;
    
    auto batch = create_identical_ions(N, pos0, vel0, mass);
    
    // Initial energy
    double E0_cpu = 0.5 * mass * (vel0.x*vel0.x + vel0.y*vel0.y + vel0.z*vel0.z);
    
    auto cpu_ensemble = IonEnsemble::from_legacy(batch.cpu_ions);
    
    // Multi-step integration
    for (int step = 0; step < steps; ++step) {
        double t = step * dt;
        boris_cpu_step_batch(*boris_cpu, cpu_ensemble, t, dt, forces);
        
        bool success = gpu_helper->integrate_batch_boris(batch.gpu_ions, dt, t, nullptr);
        REQUIRE(success);
    }
    
    batch.cpu_ions = cpu_ensemble.to_legacy();
    
    // Check energy conservation for both CPU and GPU
    auto& ion_cpu = batch.cpu_ions[0];
    auto& ion_gpu = batch.gpu_ions[0];
    
    double E_cpu = 0.5 * mass * (ion_cpu.vel.x*ion_cpu.vel.x + 
                                  ion_cpu.vel.y*ion_cpu.vel.y + 
                                  ion_cpu.vel.z*ion_cpu.vel.z);
    double E_gpu = 0.5 * mass * (ion_gpu.vel.x*ion_gpu.vel.x + 
                                  ion_gpu.vel.y*ion_gpu.vel.y + 
                                  ion_gpu.vel.z*ion_gpu.vel.z);
    
    INFO("Initial energy: " << E0_cpu << " J");
    INFO("Final CPU energy: " << E_cpu << " J");
    INFO("Final GPU energy: " << E_gpu << " J");
    
    // Both should conserve energy (free particle)
    REQUIRE_THAT(E_cpu, WithinAbs(E0_cpu, 1e-15 * E0_cpu));
    REQUIRE_THAT(E_gpu, WithinAbs(E0_cpu, 1e-15 * E0_cpu));
    
    // CPU and GPU should match
    REQUIRE_THAT(E_gpu, WithinAbs(E_cpu, 1e-15 * E0_cpu));
}

TEST_CASE("Boris: CPU/GPU parity - Timestep variation", "[boris][parity][gpu]") {
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    REQUIRE(gpu_ctx->is_valid());
    
    auto gpu_helper = icarion::gpu::GPUIntegrationHelper::create(*gpu_ctx, 100);
    auto boris_cpu = std::make_shared<BorisStrategy>();
    
    DomainConfig domain;
    domain.solver = SolverType::Boris;
    ForceRegistry forces(domain);
    
    std::vector<double> timesteps = {1e-8, 5e-8, 1e-7, 5e-7, 1e-6};
    
    for (double dt : timesteps) {
        INFO("Testing timestep: " << dt << " s");
        
        auto batch = create_identical_ions(300, Vec3{0.0, 0.0, 0.0}, Vec3{1000.0, 0.0, 0.0});
        auto cpu_ensemble = IonEnsemble::from_legacy(batch.cpu_ions);
        boris_cpu_step_batch(*boris_cpu, cpu_ensemble, 0.0, dt, forces);
        batch.cpu_ions = cpu_ensemble.to_legacy();
        
        bool success = gpu_helper->integrate_batch_boris(batch.gpu_ions, dt, 0.0, nullptr);
        REQUIRE(success);
        
        for (size_t i = 0; i < batch.cpu_ions.size(); ++i) {
            compare_ions(batch.cpu_ions[i], batch.gpu_ions[i], 1e-11, 1e-11);
        }
    }
}

#endif // ICARION_USE_GPU

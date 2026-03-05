// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifdef ICARION_USE_GPU

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/gpu/core/GPUContext.h"
#include "core/gpu/core/GPUIntegrationHelper.h"
#include "core/integrator/strategies/RK4Strategy.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/config/types/DomainConfig.h"
#include "core/types/IonEnsemble.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

#include <vector>
#include <cmath>
#include <algorithm>

using namespace ICARION;
using namespace ICARION::integrator;
using namespace ICARION::physics;
using namespace ICARION::config;
using ICARION::core::IonEnsemble;
using Catch::Matchers::WithinAbs;
using icarion::gpu::GPUContext;
using icarion::gpu::GPUIntegrationHelper;

// =============================================================================
// HELPER: Create identical ion sets for CPU/GPU comparison
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
        ion.active = true;
        ion.born = true;
        ion.t = 0.0;
        
        batch.cpu_ions[i] = ion;
        batch.gpu_ions[i] = ion;
    }
    
    return batch;
}

void advance_cpu_batch(RK4Strategy& integrator,
                       ForceRegistry& forces,
                       std::vector<IonState>& ions,
                       double t,
                       double dt) {
    IonEnsemble ensemble = IonEnsemble::from_legacy(ions);
    for (size_t i = 0; i < ions.size(); ++i) {
        if (!ensemble.is_active(i)) {
            continue;
        }
        integrator.step(ensemble, i, t, dt, forces);
    }
    ions = ensemble.to_legacy();
}

// =============================================================================
// HELPER: Compare ion states
// =============================================================================
void compare_ion_states(const IonState& cpu, const IonState& gpu, 
                       double pos_tol = 1e-6, double vel_tol = 1e-6) {
    // Position parity
    REQUIRE_THAT(gpu.pos.x, WithinAbs(cpu.pos.x, pos_tol));
    REQUIRE_THAT(gpu.pos.y, WithinAbs(cpu.pos.y, pos_tol));
    REQUIRE_THAT(gpu.pos.z, WithinAbs(cpu.pos.z, pos_tol));
    
    // Velocity parity
    REQUIRE_THAT(gpu.vel.x, WithinAbs(cpu.vel.x, vel_tol));
    REQUIRE_THAT(gpu.vel.y, WithinAbs(cpu.vel.y, vel_tol));
    REQUIRE_THAT(gpu.vel.z, WithinAbs(cpu.vel.z, vel_tol));
    
    // Time parity (exact)
    REQUIRE(gpu.t == cpu.t);
    
    // State flags
    REQUIRE(gpu.active == cpu.active);
    REQUIRE(gpu.born == cpu.born);
}

// =============================================================================
// TEST: Single Step Parity - Uniform Motion
// =============================================================================
TEST_CASE("CPU/GPU parity - single RK4 step, uniform motion", "[parity][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    // Setup GPU
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    auto gpu_helper = GPUIntegrationHelper::create(*ctx, 10);
    REQUIRE(gpu_helper != nullptr);
    
    // Setup CPU (RK4 with empty force registry)
    DomainConfig domain_cfg;
    ForceRegistry cpu_forces(domain_cfg);
    RK4Strategy cpu_integrator;
    
    SECTION("Zero forces, 100 ions (at threshold)") {
        // Test with 100 ions (exactly at GPU threshold)
        auto batch = create_identical_ions(100, Vec3{1.0, 2.0, 3.0}, Vec3{0.5, 0.0, -0.2});
        const double dt = 0.001;
        const double t0 = 0.0;
        
        // CPU integration
        advance_cpu_batch(cpu_integrator, cpu_forces, batch.cpu_ions, t0, dt);
        
        // GPU integration
        bool success = gpu_helper->integrate_batch_rk4(batch.gpu_ions, dt, t0);
        REQUIRE(success);
        
        // Compare results
        for (size_t i = 0; i < 100; ++i) {
            compare_ion_states(batch.cpu_ions[i], batch.gpu_ions[i]);
        }
    }
    
    SECTION("Zero forces, 100 ions") {
        auto batch = create_identical_ions(100, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0});
        const double dt = 0.01;
        const double t0 = 0.0;
        
        // CPU integration (per-ion loop)
        advance_cpu_batch(cpu_integrator, cpu_forces, batch.cpu_ions, t0, dt);
        
        // GPU integration (batch)
        bool success = gpu_helper->integrate_batch_rk4(batch.gpu_ions, dt, t0);
        REQUIRE(success);
        
        // Compare all ions
        for (size_t i = 0; i < 100; ++i) {
            compare_ion_states(batch.cpu_ions[i], batch.gpu_ions[i]);
        }
    }
}

// =============================================================================
// TEST: Multi-Step Parity
// =============================================================================
TEST_CASE("CPU/GPU parity - multi-step trajectory", "[parity][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    auto gpu_helper = GPUIntegrationHelper::create(*ctx, 10);
    REQUIRE(gpu_helper != nullptr);
    
    DomainConfig domain_cfg;
    ForceRegistry cpu_forces(domain_cfg);
    RK4Strategy cpu_integrator;
    
    SECTION("100 steps, 50 ions") {
        auto batch = create_identical_ions(50, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.5, 0.0});
        const double dt = 0.001;
        const int num_steps = 100;
        
        for (int step = 0; step < num_steps; ++step) {
            double t = step * dt;
            
            // CPU per-ion
            advance_cpu_batch(cpu_integrator, cpu_forces, batch.cpu_ions, t, dt);
            
            // GPU batch
            bool success = gpu_helper->integrate_batch_rk4(batch.gpu_ions, dt, t);
            REQUIRE(success);
        }
        
        // Compare final states
        for (size_t i = 0; i < 50; ++i) {
            // After 100 steps, accumulated error might be larger
            compare_ion_states(batch.cpu_ions[i], batch.gpu_ions[i], 1e-5, 1e-5);
        }
    }
}

// =============================================================================
// TEST: Varied Timesteps
// =============================================================================
TEST_CASE("CPU/GPU parity - different timestep sizes", "[parity][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    auto gpu_helper = GPUIntegrationHelper::create(*ctx, 10);
    REQUIRE(gpu_helper != nullptr);
    
    DomainConfig domain_cfg;
    ForceRegistry cpu_forces(domain_cfg);
    RK4Strategy cpu_integrator;
    
    std::vector<double> timesteps = {1e-3, 1e-6, 1e-9, 1e-12};
    
    for (double dt : timesteps) {
        SECTION("dt = " + std::to_string(dt)) {
            auto batch = create_identical_ions(20, Vec3{1.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0});
            
            // CPU
            advance_cpu_batch(cpu_integrator, cpu_forces, batch.cpu_ions, 0.0, dt);
            
            // GPU
            bool success = gpu_helper->integrate_batch_rk4(batch.gpu_ions, dt, 0.0);
            REQUIRE(success);
            
            // Compare
            for (size_t i = 0; i < 20; ++i) {
                compare_ion_states(batch.cpu_ions[i], batch.gpu_ions[i]);
            }
        }
    }
}

// =============================================================================
// TEST: Varied Initial Conditions
// =============================================================================
TEST_CASE("CPU/GPU parity - diverse initial conditions", "[parity][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    auto gpu_helper = GPUIntegrationHelper::create(*ctx, 10);
    REQUIRE(gpu_helper != nullptr);
    
    DomainConfig domain_cfg;
    ForceRegistry cpu_forces(domain_cfg);
    RK4Strategy cpu_integrator;
    
    SECTION("Mixed positions and velocities") {
        const size_t n_ions = 30;
        auto batch = create_identical_ions(n_ions, Vec3{0, 0, 0}, Vec3{0, 0, 0});
        
        // Set diverse initial conditions
        for (size_t i = 0; i < n_ions; ++i) {
            double x = 0.01 * i;
            double vx = 1.0 + 0.1 * i;
            
            batch.cpu_ions[i].pos = Vec3{x, 0, 0};
            batch.cpu_ions[i].vel = Vec3{vx, 0, 0};
            
            batch.gpu_ions[i].pos = Vec3{x, 0, 0};
            batch.gpu_ions[i].vel = Vec3{vx, 0, 0};
        }
        
        const double dt = 0.01;
        
        // CPU
        advance_cpu_batch(cpu_integrator, cpu_forces, batch.cpu_ions, 0.0, dt);
        
        // GPU
        bool success = gpu_helper->integrate_batch_rk4(batch.gpu_ions, dt, 0.0);
        REQUIRE(success);
        
        // Compare all
        for (size_t i = 0; i < n_ions; ++i) {
            compare_ion_states(batch.cpu_ions[i], batch.gpu_ions[i]);
        }
    }
}

// =============================================================================
// TEST: Inactive Ions Handling
// =============================================================================
TEST_CASE("CPU/GPU parity - inactive ions", "[parity][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    auto gpu_helper = GPUIntegrationHelper::create(*ctx, 10);
    REQUIRE(gpu_helper != nullptr);
    
    DomainConfig domain_cfg;
    ForceRegistry cpu_forces(domain_cfg);
    RK4Strategy cpu_integrator;
    
    SECTION("Mix of active and inactive ions") {
        auto batch = create_identical_ions(50, Vec3{1.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        
        // Mark every 5th ion as inactive
        for (size_t i = 0; i < 50; i += 5) {
            batch.cpu_ions[i].active = false;
            batch.gpu_ions[i].active = false;
        }
        
        const double dt = 0.01;
        
        // CPU (RK4Strategy skips inactive ions)
        advance_cpu_batch(cpu_integrator, cpu_forces, batch.cpu_ions, 0.0, dt);
        
        // GPU (kernel skips inactive ions)
        bool success = gpu_helper->integrate_batch_rk4(batch.gpu_ions, dt, 0.0);
        REQUIRE(success);
        
        // Compare all (inactive should remain unchanged)
        for (size_t i = 0; i < 50; ++i) {
            compare_ion_states(batch.cpu_ions[i], batch.gpu_ions[i]);
        }
    }
}

// =============================================================================
// TEST: Numerical Stability
// =============================================================================
TEST_CASE("CPU/GPU parity - numerical stability", "[parity][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    auto gpu_helper = GPUIntegrationHelper::create(*ctx, 10);
    REQUIRE(gpu_helper != nullptr);
    
    DomainConfig domain_cfg;
    ForceRegistry cpu_forces(domain_cfg);
    RK4Strategy cpu_integrator;
    
    SECTION("Very small values (near zero)") {
        auto batch = create_identical_ions(20, 
                                          Vec3{1e-12, 1e-12, 1e-12}, 
                                          Vec3{1e-12, 0, 0});
        const double dt = 1e-9;
        
        // CPU
        advance_cpu_batch(cpu_integrator, cpu_forces, batch.cpu_ions, 0.0, dt);
        
        // GPU
        bool success = gpu_helper->integrate_batch_rk4(batch.gpu_ions, dt, 0.0);
        REQUIRE(success);
        
        // Compare with tight tolerance
        for (size_t i = 0; i < 20; ++i) {
            compare_ion_states(batch.cpu_ions[i], batch.gpu_ions[i], 1e-15, 1e-15);
        }
    }
    
    SECTION("Large values") {
        auto batch = create_identical_ions(20, 
                                          Vec3{1e6, 1e6, 1e6}, 
                                          Vec3{1e3, 0, 0});
        const double dt = 1e-6;
        
        // CPU
        advance_cpu_batch(cpu_integrator, cpu_forces, batch.cpu_ions, 0.0, dt);
        
        // GPU
        bool success = gpu_helper->integrate_batch_rk4(batch.gpu_ions, dt, 0.0);
        REQUIRE(success);
        
        // Compare (relative error due to large magnitudes)
        for (size_t i = 0; i < 20; ++i) {
            double rel_tol = 1e-9;  // Relative tolerance
            REQUIRE_THAT(batch.gpu_ions[i].pos.x, 
                        WithinAbs(batch.cpu_ions[i].pos.x, 
                                 std::abs(batch.cpu_ions[i].pos.x) * rel_tol + 1e-6));
        }
    }
}

// =============================================================================
// TEST: Performance Comparison (Informational)
// =============================================================================
TEST_CASE("CPU/GPU performance comparison", "[parity][integration][.performance]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    auto gpu_helper = GPUIntegrationHelper::create(*ctx, 1000);
    REQUIRE(gpu_helper != nullptr);
    
    DomainConfig domain_cfg;
    ForceRegistry cpu_forces(domain_cfg);
    RK4Strategy cpu_integrator;
    
    SECTION("Large batch (10k ions, 100 steps)") {
        const size_t n_ions = 10000;
        const int n_steps = 100;
        const double dt = 1e-9;
        
        auto batch = create_identical_ions(n_ions, Vec3{0, 0, 0}, Vec3{1, 0, 0});
        
        // CPU timing
        auto cpu_start = std::chrono::high_resolution_clock::now();
        for (int step = 0; step < n_steps; ++step) {
            advance_cpu_batch(cpu_integrator, cpu_forces, batch.cpu_ions, step * dt, dt);
        }
        auto cpu_end = std::chrono::high_resolution_clock::now();
        double cpu_ms = std::chrono::duration<double, std::milli>(cpu_end - cpu_start).count();
        
        // GPU timing
        auto gpu_start = std::chrono::high_resolution_clock::now();
        for (int step = 0; step < n_steps; ++step) {
            gpu_helper->integrate_batch_rk4(batch.gpu_ions, dt, step * dt);
        }
        auto gpu_end = std::chrono::high_resolution_clock::now();
        double gpu_ms = std::chrono::duration<double, std::milli>(gpu_end - gpu_start).count();
        
        // Verify parity despite timing
        for (size_t i = 0; i < std::min(n_ions, size_t(100)); ++i) {
            compare_ion_states(batch.cpu_ions[i], batch.gpu_ions[i], 1e-4, 1e-4);
        }
        
        // Report timing (informational)
        INFO("CPU time: " << cpu_ms << " ms");
        INFO("GPU time: " << gpu_ms << " ms");
        INFO("Speedup: " << cpu_ms / gpu_ms << "x");
        
        // GPU should be faster for large batches
        REQUIRE(gpu_ms < cpu_ms);
    }
}

#endif  // ICARION_USE_GPU

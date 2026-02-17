// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#ifdef ICARION_USE_GPU

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "core/gpu/core/GPUContext.h"
#include "core/gpu/core/GPUIntegrationHelper.h"
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

#include <vector>
#include <cmath>

using namespace ICARION;
using Catch::Matchers::WithinAbs;
using icarion::gpu::GPUContext;
using icarion::gpu::GPUIntegrationHelper;

// =============================================================================
// HELPER: Create test ions
// =============================================================================
std::vector<IonState> create_test_ions(size_t count, Vec3 pos0, Vec3 vel0, double mass = 1.0e-26) {
    std::vector<IonState> ions(count);
    for (size_t i = 0; i < count; ++i) {
        ions[i].pos = pos0;
        ions[i].vel = vel0;
        ions[i].mass_kg = mass;
        ions[i].ion_charge_C = 1.602e-19;  // Single proton charge
        ions[i].active = true;
        ions[i].born = true;
        ions[i].t = 0.0;
    }
    return ions;
}

// =============================================================================
// TEST: GPU Context Creation
// =============================================================================
TEST_CASE("GPUContext creation and properties", "[gpu][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    SECTION("Create context for device 0") {
        auto ctx = GPUContext::create(0);
        REQUIRE(ctx != nullptr);
        
        const auto& props = ctx->get_properties();
        REQUIRE(props.device_id == 0);
        REQUIRE(props.compute_capability_major > 0);
        REQUIRE(props.total_memory > 0);
        REQUIRE(!props.name.empty());
    }
}

// =============================================================================
// TEST: GPUIntegrationHelper Creation
// =============================================================================
TEST_CASE("GPUIntegrationHelper creation", "[gpu][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    
    SECTION("Create helper with default threshold") {
        auto helper = GPUIntegrationHelper::create(*ctx, 5000);
        REQUIRE(helper != nullptr);
    }
    
    SECTION("Create helper with custom threshold") {
        auto helper = GPUIntegrationHelper::create(*ctx, 1000);
        REQUIRE(helper != nullptr);
    }
}

// =============================================================================
// TEST: Batch Integration - Uniform Motion (Zero Forces)
// =============================================================================
TEST_CASE("GPU batch integration - uniform motion", "[gpu][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    
    auto helper = GPUIntegrationHelper::create(*ctx, 100);
    REQUIRE(helper != nullptr);
    
    SECTION("Single timestep, zero forces") {
        // Initial conditions: x=1m, v=2m/s, a=0
        auto ions = create_test_ions(100, Vec3{1.0, 0.0, 0.0}, Vec3{2.0, 0.0, 0.0});
        const double dt = 0.01;  // 10ms
        const double t0 = 0.0;
        
        // Integrate (with zero fields, acceleration = 0)
        bool success = helper->integrate_batch_rk4(ions, dt, t0);
        REQUIRE(success);
        
        // Update ion times (GPUIntegrationHelper doesn't do this)
        for (auto& ion : ions) {
            ion.t = t0 + dt;
        }
        
        // Expected: x = x0 + v*dt = 1.0 + 2.0*0.01 = 1.02m
        for (const auto& ion : ions) {
            REQUIRE_THAT(ion.pos.x, WithinAbs(1.02, 1e-6));
            REQUIRE_THAT(ion.pos.y, WithinAbs(0.0, 1e-9));
            REQUIRE_THAT(ion.pos.z, WithinAbs(0.0, 1e-9));
            
            // Velocity unchanged (no forces)
            REQUIRE_THAT(ion.vel.x, WithinAbs(2.0, 1e-6));
            REQUIRE_THAT(ion.vel.y, WithinAbs(0.0, 1e-9));
            REQUIRE_THAT(ion.vel.z, WithinAbs(0.0, 1e-9));
            
            // Time updated
            REQUIRE_THAT(ion.t, WithinAbs(t0 + dt, 1e-9));
        }
    }
    
    SECTION("Multiple timesteps") {
        auto ions = create_test_ions(100, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        const double dt = 0.01;
        
        // Integrate 10 steps
        for (int step = 0; step < 10; ++step) {
            double t = step * dt;
            bool success = helper->integrate_batch_rk4(ions, dt, t);
            REQUIRE(success);
            
            // Update times (caller's responsibility)
            for (auto& ion : ions) {
                ion.t = t + dt;
            }
            
            // Check position after each step
            INFO("Step " << step << ": pos.x = " << ions[0].pos.x);
        }
        
        // Expected: x = v*t = 1.0 * 0.1 = 0.1m
        for (const auto& ion : ions) {
            REQUIRE_THAT(ion.pos.x, WithinAbs(0.1, 1e-5));
            REQUIRE_THAT(ion.t, WithinAbs(0.1, 1e-9));
        }
    }
}

// =============================================================================
// TEST: Edge Cases
// =============================================================================
TEST_CASE("GPU batch integration - edge cases", "[gpu][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    
    auto helper = GPUIntegrationHelper::create(*ctx, 100);
    REQUIRE(helper != nullptr);
    
    SECTION("Empty batch") {
        std::vector<IonState> empty_ions;
        bool success = helper->integrate_batch_rk4(empty_ions, 0.01, 0.0);
        // Empty batch returns false (expected - below threshold)
        REQUIRE(!success);
    }
    
    SECTION("Single ion - below threshold") {
        // Single ion is below threshold (100), so GPU returns false
        auto ions = create_test_ions(1, Vec3{1.0, 2.0, 3.0}, Vec3{0.5, 0.0, 0.0});
        bool success = helper->integrate_batch_rk4(ions, 0.01, 0.0);
        // Returns false because N < threshold
        REQUIRE(!success);
    }
    
    SECTION("At threshold - 100 ions") {
        auto ions = create_test_ions(100, Vec3{1.0, 0.0, 0.0}, Vec3{0.5, 0.0, 0.0});
        bool success = helper->integrate_batch_rk4(ions, 0.01, 0.0);
        REQUIRE(success);
        
        REQUIRE_THAT(ions[0].pos.x, WithinAbs(1.005, 1e-6));
        REQUIRE_THAT(ions[0].pos.y, WithinAbs(0.0, 1e-9));
        REQUIRE_THAT(ions[0].pos.z, WithinAbs(0.0, 1e-9));
    }
    
    SECTION("Large batch (above threshold)") {
        auto ions = create_test_ions(10000, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 1.0, 1.0});
        bool success = helper->integrate_batch_rk4(ions, 0.01, 0.0);
        REQUIRE(success);
        
        // Verify all ions integrated
        for (const auto& ion : ions) {
            REQUIRE_THAT(ion.pos.x, WithinAbs(0.01, 1e-6));
            REQUIRE_THAT(ion.pos.y, WithinAbs(0.01, 1e-6));
            REQUIRE_THAT(ion.pos.z, WithinAbs(0.01, 1e-6));
        }
    }
    
    SECTION("Inactive ions preserved") {
        auto ions = create_test_ions(100, Vec3{1.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        ions[50].active = false;
        Vec3 inactive_pos = ions[50].pos;
        
        bool success = helper->integrate_batch_rk4(ions, 0.01, 0.0);
        REQUIRE(success);
        
        // Inactive ion should not move
        REQUIRE(ions[50].pos.x == inactive_pos.x);
        REQUIRE(ions[50].pos.y == inactive_pos.y);
        REQUIRE(ions[50].pos.z == inactive_pos.z);
    }
}

// =============================================================================
// TEST: Statistics
// =============================================================================
TEST_CASE("GPU integration statistics", "[gpu][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    
    auto helper = GPUIntegrationHelper::create(*ctx, 100);
    REQUIRE(helper != nullptr);
    
    SECTION("Stats after integration") {
        auto ions = create_test_ions(1000, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        
        // Integrate multiple times
        for (int i = 0; i < 5; ++i) {
            helper->integrate_batch_rk4(ions, 0.01, i * 0.01);
        }
        
        auto stats = helper->get_stats();
        REQUIRE(stats.gpu_integrations == 5);
        REQUIRE(stats.total_ions_gpu == 5000);
        REQUIRE(stats.total_time_ms > 0.0);
    }
}

// =============================================================================
// TEST: Memory Reuse
// =============================================================================
TEST_CASE("GPU memory buffer reuse", "[gpu][integration]") {
    if (!GPUContext::is_cuda_available()) {
        SKIP("CUDA not available");
    }
    
    auto ctx = GPUContext::create(0);
    REQUIRE(ctx != nullptr);
    
    auto helper = GPUIntegrationHelper::create(*ctx, 100);
    REQUIRE(helper != nullptr);
    
    SECTION("Multiple batches with same size") {
        auto ions = create_test_ions(500, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        
        // First batch allocates
        bool success1 = helper->integrate_batch_rk4(ions, 0.01, 0.0);
        REQUIRE(success1);
        
        // Second batch reuses
        bool success2 = helper->integrate_batch_rk4(ions, 0.01, 0.01);
        REQUIRE(success2);
        
        auto stats = helper->get_stats();
        REQUIRE(stats.gpu_integrations == 2);
    }
    
    SECTION("Growing batch sizes") {
        // Start small
        auto ions = create_test_ions(100, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        bool success1 = helper->integrate_batch_rk4(ions, 0.01, 0.0);
        REQUIRE(success1);
        
        // Grow to 500 (should reallocate)
        ions = create_test_ions(500, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        bool success2 = helper->integrate_batch_rk4(ions, 0.01, 0.0);
        REQUIRE(success2);
        
        // Back to 300 (should reuse 500-sized buffer)
        ions = create_test_ions(300, Vec3{0.0, 0.0, 0.0}, Vec3{1.0, 0.0, 0.0});
        bool success3 = helper->integrate_batch_rk4(ions, 0.01, 0.0);
        REQUIRE(success3);
        
        auto stats = helper->get_stats();
        REQUIRE(stats.gpu_integrations == 3);
    }
}

#endif  // ICARION_USE_GPU

// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
/**
 * @file test_gpu_boundaries.cpp
 * @brief GPU boundary checking CPU/GPU parity test
 * 
 * Validates that GPU boundary checks match CPU results exactly.
 */

#include <catch2/catch_all.hpp>
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/core/GPUIntegrationHelper.h"
#include "core/types/IonState.h"
#include <vector>
#include <cmath>

using namespace ICARION;
using namespace ICARION::core;

// CPU reference implementation
void check_boundaries_cpu(
    std::vector<IonState>& ions,
    double length_m,
    double radius_m,
    bool is_last_domain
) {
    constexpr double EPSILON = 1e-12;
    
    for (auto& ion : ions) {
        if (!ion.active) continue;
        
        // Axial lower bound
        if (ion.pos.z < -EPSILON) {
            ion.active = false;
            continue;
        }
        
        // Axial upper bound
        bool z_valid;
        if (is_last_domain) {
            z_valid = (ion.pos.z < length_m);
        } else {
            z_valid = (ion.pos.z <= length_m + EPSILON);
        }
        
        if (!z_valid) {
            ion.active = false;
            continue;
        }
        
        // Radial bound
        double r = std::sqrt(ion.pos.x*ion.pos.x + ion.pos.y*ion.pos.y);
        if (r > radius_m + EPSILON) {
            ion.active = false;
            continue;
        }
    }
}

TEST_CASE("GPU Boundary Checking", "[gpu][boundary]") {
    // Create GPU context
    auto gpu_ctx = icarion::gpu::GPUContext::create();
    REQUIRE(gpu_ctx != nullptr);
    REQUIRE(gpu_ctx->is_valid());
    
    // Create helper
    auto gpu_helper = icarion::gpu::GPUIntegrationHelper::create(*gpu_ctx, 100);
    REQUIRE(gpu_helper != nullptr);
    REQUIRE(gpu_helper->is_enabled());
    
    // Domain geometry
    const double LENGTH = 0.1;   // 10 cm
    const double RADIUS = 0.01;  // 1 cm
    
    SECTION("All ions inside boundaries") {
        const int N = 1000;
        
        std::vector<IonState> ions_cpu;
        std::vector<IonState> ions_gpu;
        
        // Create ions inside domain (well within boundaries)
        for (int i = 0; i < N; ++i) {
            IonState ion;
            ion.pos = {
                0.003 * std::sin(2.0 * M_PI * i / N),  // r < 0.003 << 0.01 (radius)
                0.003 * std::cos(2.0 * M_PI * i / N),
                0.02 + 0.06 * i / N  // z ∈ [0.02, 0.08] ⊂ [0, 0.1]
            };
            ion.vel = {0, 0, 0};
            ion.active = true;
            ion.mass_kg = 19.0 * 1.66054e-27;
            ion.ion_charge_C = 1.60218e-19;
            
            ions_cpu.push_back(ion);
            ions_gpu.push_back(ion);
        }
        
        // CPU check
        check_boundaries_cpu(ions_cpu, LENGTH, RADIUS, true);
        
        // GPU check
        bool gpu_success = gpu_helper->check_boundaries_batch(ions_gpu, LENGTH, RADIUS, true);
        REQUIRE(gpu_success);
        
        // Compare results
        int cpu_active = 0, gpu_active = 0;
        for (int i = 0; i < N; ++i) {
            REQUIRE(ions_cpu[i].active == ions_gpu[i].active);
            if (ions_cpu[i].active) ++cpu_active;
            if (ions_gpu[i].active) ++gpu_active;
        }
        
        REQUIRE(cpu_active == N);  // All should be active
        REQUIRE(gpu_active == N);
    }
    
    SECTION("Ions outside axial bounds") {
        const int N = 100;
        
        std::vector<IonState> ions_cpu;
        std::vector<IonState> ions_gpu;
        
        // Create ions at various z positions
        for (int i = 0; i < N; ++i) {
            IonState ion;
            ion.pos = {0.0, 0.0, -0.01 + 0.12 * i / N};  // From z=-0.01 to z=0.11
            ion.vel = {0, 0, 0};
            ion.active = true;
            ion.mass_kg = 19.0 * 1.66054e-27;
            ion.ion_charge_C = 1.60218e-19;
            
            ions_cpu.push_back(ion);
            ions_gpu.push_back(ion);
        }
        
        // CPU check
        check_boundaries_cpu(ions_cpu, LENGTH, RADIUS, true);
        
        // GPU check
        bool gpu_success = gpu_helper->check_boundaries_batch(ions_gpu, LENGTH, RADIUS, true);
        REQUIRE(gpu_success);
        
        // Compare results
        for (int i = 0; i < N; ++i) {
            REQUIRE(ions_cpu[i].active == ions_gpu[i].active);
        }
        
        // Count active ions
        int cpu_active = 0, gpu_active = 0;
        for (int i = 0; i < N; ++i) {
            if (ions_cpu[i].active) ++cpu_active;
            if (ions_gpu[i].active) ++gpu_active;
        }
        
        REQUIRE(cpu_active > 0);   // Some inside
        REQUIRE(cpu_active < N);   // Some outside
        REQUIRE(gpu_active == cpu_active);
    }
    
    SECTION("Ions outside radial bounds") {
        const int N = 100;
        
        std::vector<IonState> ions_cpu;
        std::vector<IonState> ions_gpu;
        
        // Create ions at various radial positions
        for (int i = 0; i < N; ++i) {
            double r = 0.015 * i / N;  // From r=0 to r=0.015 (radius=0.01)
            IonState ion;
            ion.pos = {r, 0.0, 0.05};
            ion.vel = {0, 0, 0};
            ion.active = true;
            ion.mass_kg = 19.0 * 1.66054e-27;
            ion.ion_charge_C = 1.60218e-19;
            
            ions_cpu.push_back(ion);
            ions_gpu.push_back(ion);
        }
        
        // CPU check
        check_boundaries_cpu(ions_cpu, LENGTH, RADIUS, true);
        
        // GPU check
        bool gpu_success = gpu_helper->check_boundaries_batch(ions_gpu, LENGTH, RADIUS, true);
        REQUIRE(gpu_success);
        
        // Compare results
        for (int i = 0; i < N; ++i) {
            REQUIRE(ions_cpu[i].active == ions_gpu[i].active);
        }
        
        // Count active ions
        int cpu_active = 0, gpu_active = 0;
        for (int i = 0; i < N; ++i) {
            if (ions_cpu[i].active) ++cpu_active;
            if (ions_gpu[i].active) ++gpu_active;
        }
        
        REQUIRE(cpu_active > 0);   // Some inside
        REQUIRE(cpu_active < N);   // Some outside
        REQUIRE(gpu_active == cpu_active);
    }
    
    SECTION("Last domain vs transition domain") {
        const int N = 100;
        
        // Create ions exactly at z=length (boundary case)
        std::vector<IonState> ions_last;
        std::vector<IonState> ions_transition;
        
        for (int i = 0; i < N; ++i) {
            IonState ion;
            ion.pos = {0.0, 0.0, LENGTH + 1e-14 * i};  // Slightly beyond length
            ion.vel = {0, 0, 0};
            ion.active = true;
            ion.mass_kg = 19.0 * 1.66054e-27;
            ion.ion_charge_C = 1.60218e-19;
            
            ions_last.push_back(ion);
            ions_transition.push_back(ion);
        }
        
        // Last domain (strict check)
        check_boundaries_cpu(ions_last, LENGTH, RADIUS, true);
        
        // Transition domain (lenient check)
        check_boundaries_cpu(ions_transition, LENGTH, RADIUS, false);
        
        // Count active
        int last_active = 0, transition_active = 0;
        for (int i = 0; i < N; ++i) {
            if (ions_last[i].active) ++last_active;
            if (ions_transition[i].active) ++transition_active;
        }
        
        // Transition domain should allow more ions (within epsilon)
        REQUIRE(transition_active >= last_active);
    }
}

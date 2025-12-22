// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file test_gpu_field_interpolation.cpp
 * @brief Test GPU field interpolation with real field arrays
 * 
 * Simplified test focusing on field provider extraction and basic integration.
 */

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "core/gpu/core/GPUIntegrationHelper.h"
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/fields/FieldArrayGPU_conversion.h"
#include "fieldsolver/utils/GridFieldProvider.h"
#include "core/io/fieldArrayLoader.h"
#include "core/types/IonState.h"
#include "utils/constants.h"
#include <cmath>
#include <vector>
#include <memory>

using namespace ICARION::core;
using namespace icarion::gpu;

TEST_CASE("GPU Field Interpolation: Uniform Field Analytical", "[gpu][fields][integration]") {
    // Create GPU context
    if (!GPUContext::is_cuda_available()) {
        SKIP("GPU not available");
    }
    
    auto gpu_context = GPUContext::create(0);
    REQUIRE(gpu_context != nullptr);
    
    auto gpu_helper = GPUIntegrationHelper::create(*gpu_context, 10);  // Low threshold for testing
    REQUIRE(gpu_helper != nullptr);
    
    // Load uniform field array (Ez = 1000 V/m)
    const std::string field_file = "examples/field_arrays/uniform_field.h5";
    
    FieldArray field_array = load_field_array(field_file);
    REQUIRE(field_array.nx > 0);  // Check if loaded
    
    // Verify field is valid for GPU
    REQUIRE(is_field_valid_for_gpu(field_array));
    
    // Create field provider
    auto field_provider = std::make_unique<GridFieldProvider>(&field_array);
    
    // Create test ions: 100 H3O+ (all at same position for uniform field)
    const int N = 100;
    std::vector<IonState> ions(N);
    for (auto& ion : ions) {
        ion.pos = {0.0, 0.0, 0.01};  // Center of field (z=1cm)
        ion.vel = {0.0, 0.0, 0.0};   // Start at rest
        ion.mass_kg = 19.0 * AMU_TO_KG;
        ion.ion_charge_C = ELEM_CHARGE_C;
        ion.active = true;
        ion.t = 0.0;
        ion.species_id = "H3O+";
        ion.current_domain_index = 0;
    }
    
    // Analytical expectation: F = qE, a = qE/m, z(t) = z0 + 0.5*a*t^2
    const double E_z = 1000.0;  // V/m
    const double q = ELEM_CHARGE_C;
    const double m = 19.0 * AMU_TO_KG;
    const double a_z = q * E_z / m;  // Acceleration [m/s^2]
    
    // Integrate one step (dt = 1 us)
    const double dt = 1e-6;  // 1 microsecond
    const double t = 0.0;
    
    // GPU integration with field
    bool success = gpu_helper->integrate_batch_rk4(ions, dt, t, field_provider.get());
    REQUIRE(success);
    
    // Expected position: z = z0 + 0.5*a*dt^2
    const double z_expected = 0.01 + 0.5 * a_z * dt * dt;
    
    // Expected velocity: v = a*dt
    const double vz_expected = a_z * dt;
    
    // Check results (relaxed tolerances for texture interpolation + float precision)
    REQUIRE_THAT(ions[0].pos.z, Catch::Matchers::WithinAbs(z_expected, 1e-8));
    REQUIRE_THAT(ions[0].vel.z, Catch::Matchers::WithinAbs(vz_expected, 1e-6));
    
    // x and y should remain unchanged (Ex = Ey = 0)
    REQUIRE_THAT(ions[0].pos.x, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(ions[0].pos.y, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(ions[0].vel.x, Catch::Matchers::WithinAbs(0.0, 1e-12));
    REQUIRE_THAT(ions[0].vel.y, Catch::Matchers::WithinAbs(0.0, 1e-12));
}

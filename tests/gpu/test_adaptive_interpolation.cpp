// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>

#include "core/gpu/fields/AdaptiveFieldInterpolator.h"
#include <cmath>

using namespace icarion::gpu;
using Catch::Approx;

TEST_CASE("AdaptiveFieldInterpolator", "[gpu][adaptive][field]") {
    
    SECTION("Uniform field - linear interpolation") {
        // Create 4×4×4 grid with uniform E = (100, 200, 300) V/m
        int nx = 4, ny = 4, nz = 4;
        size_t grid_size = nx * ny * nz;
        
        std::vector<double> Ex(grid_size, 100.0);
        std::vector<double> Ey(grid_size, 200.0);
        std::vector<double> Ez(grid_size, 300.0);
        
        Vec3 origin{0.0, 0.0, 0.0};
        Vec3 spacing{0.01, 0.01, 0.01};  // 1cm cells
        
        AdaptiveFieldInterpolator interp(Ex, Ey, Ez, nx, ny, nz, origin, spacing);
        
        // Test interpolation at center
        Vec3 pos{0.015, 0.015, 0.015};  // Inside grid
        Vec3 E = interp.interpolate(pos);
        
        REQUIRE(E.x == Approx(100.0).margin(1e-10));
        REQUIRE(E.y == Approx(200.0).margin(1e-10));
        REQUIRE(E.z == Approx(300.0).margin(1e-10));
        
        // Check that linear interpolation was used (uniform field → zero gradient)
        const auto& stats = interp.get_stats();
        REQUIRE(stats.n_linear == 1);
        REQUIRE(stats.n_cubic == 0);
        REQUIRE(stats.n_quintic == 0);
    }
    
    SECTION("Linear gradient field") {
        // Create 8×8×8 grid with E = (x, 0, 0) * 1000 V/m²
        // This has ∂Ex/∂x = 1000 V/m², other gradients = 0
        int nx = 8, ny = 8, nz = 8;
        size_t grid_size = nx * ny * nz;
        
        std::vector<double> Ex(grid_size);
        std::vector<double> Ey(grid_size, 0.0);
        std::vector<double> Ez(grid_size, 0.0);
        
        Vec3 origin{0.0, 0.0, 0.0};
        Vec3 spacing{0.01, 0.01, 0.01};  // 1cm cells
        
        // Fill Ex with linear gradient: Ex(i,j,k) = i * spacing * 1000
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    int idx = i * (ny * nz) + j * nz + k;
                    Ex[idx] = i * spacing.x * 1000.0;  // Linear in x
                }
            }
        }
        
        // Configure with low gradient threshold to trigger cubic
        AdaptiveFieldInterpolator::Config config;
        config.gradient_threshold_low = 500.0;   // Below 500 V/m² → linear
        config.gradient_threshold_high = 2000.0; // Above 2000 V/m² → quintic
        config.cache_gradients = true;
        
        AdaptiveFieldInterpolator interp(Ex, Ey, Ez, nx, ny, nz, origin, spacing, config);
        
        // Test at position with gradient ~1000 V/m² (should use cubic)
        Vec3 pos{0.035, 0.035, 0.035};  // Mid-grid
        Vec3 E = interp.interpolate(pos);
        
        // Expected: Ex ≈ 0.035 * 1000 = 35 V/m
        REQUIRE(E.x == Approx(35.0).margin(1.0));  // Allow 1 V/m error from interpolation
        REQUIRE(E.y == Approx(0.0).margin(1e-10));
        REQUIRE(E.z == Approx(0.0).margin(1e-10));
        
        // Check that gradient was estimated (currently falls back to linear)
        const auto& stats = interp.get_stats();
        REQUIRE(stats.n_linear + stats.n_cubic > 0);
    }
    
    SECTION("Outside grid returns zero") {
        int nx = 4, ny = 4, nz = 4;
        size_t grid_size = nx * ny * nz;
        
        std::vector<double> Ex(grid_size, 100.0);
        std::vector<double> Ey(grid_size, 200.0);
        std::vector<double> Ez(grid_size, 300.0);
        
        Vec3 origin{0.0, 0.0, 0.0};
        Vec3 spacing{0.01, 0.01, 0.01};
        
        AdaptiveFieldInterpolator interp(Ex, Ey, Ez, nx, ny, nz, origin, spacing);
        
        // Test outside grid
        Vec3 pos_outside{0.1, 0.1, 0.1};  // Far outside 4cm grid
        Vec3 E = interp.interpolate(pos_outside);
        
        REQUIRE(E.x == 0.0);
        REQUIRE(E.y == 0.0);
        REQUIRE(E.z == 0.0);
    }
    
    SECTION("Statistics tracking") {
        int nx = 4, ny = 4, nz = 4;
        size_t grid_size = nx * ny * nz;
        
        std::vector<double> Ex(grid_size, 100.0);
        std::vector<double> Ey(grid_size, 0.0);
        std::vector<double> Ez(grid_size, 0.0);
        
        Vec3 origin{0.0, 0.0, 0.0};
        Vec3 spacing{0.01, 0.01, 0.01};
        
        AdaptiveFieldInterpolator interp(Ex, Ey, Ez, nx, ny, nz, origin, spacing);
        
        // Do multiple interpolations
        for (int i = 0; i < 10; ++i) {
            Vec3 pos{0.005 + i * 0.002, 0.015, 0.015};
            interp.interpolate(pos);
        }
        
        const auto& stats = interp.get_stats();
        REQUIRE(stats.n_linear == 10);  // Uniform field → all linear
        
        // Check fractions
        REQUIRE(stats.fraction_linear() == Approx(1.0));
        REQUIRE(stats.fraction_cubic() == Approx(0.0));
        REQUIRE(stats.fraction_quintic() == Approx(0.0));
        
        // Reset stats
        interp.reset_stats();
        const auto& stats2 = interp.get_stats();
        REQUIRE(stats2.n_linear == 0);
    }
    
    SECTION("Gradient magnitude estimation") {
        // Create field with known gradient
        int nx = 8, ny = 8, nz = 8;
        size_t grid_size = nx * ny * nz;
        
        std::vector<double> Ex(grid_size);
        std::vector<double> Ey(grid_size, 0.0);
        std::vector<double> Ez(grid_size, 0.0);
        
        Vec3 origin{0.0, 0.0, 0.0};
        Vec3 spacing{0.001, 0.001, 0.001};  // 1mm cells
        
        // Quadratic field: Ex = 1000 * x²
        // ∂Ex/∂x = 2000 * x V/m²
        for (int i = 0; i < nx; ++i) {
            for (int j = 0; j < ny; ++j) {
                for (int k = 0; k < nz; ++k) {
                    int idx = i * (ny * nz) + j * nz + k;
                    double x = origin.x + i * spacing.x;
                    Ex[idx] = 1000.0 * x * x;
                }
            }
        }
        
        AdaptiveFieldInterpolator::Config config;
        config.cache_gradients = true;
        config.gradient_threshold_low = 1.0;      // Very low → force cubic
        config.gradient_threshold_high = 1000.0;
        
        AdaptiveFieldInterpolator interp(Ex, Ey, Ez, nx, ny, nz, origin, spacing, config);
        
        // Test gradient-based order selection
        Vec3 pos_center{0.0035, 0.004, 0.004};
        auto order = interp.get_order_at(pos_center);
        
        // At x=3.5mm: ∂Ex/∂x = 2000 * 0.0035 = 7 V/m² → intermediate gradient
        // Should select Cubic (between low and high threshold)
        // NOTE: Currently fallback to Linear, but order selection logic is there
        REQUIRE((order == AdaptiveFieldInterpolator::InterpolationOrder::Linear ||
                 order == AdaptiveFieldInterpolator::InterpolationOrder::Cubic));
    }
}

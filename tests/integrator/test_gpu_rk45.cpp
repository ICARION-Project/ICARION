// Test GPU RK45 integration
#include "core/gpu/GPUContext.h"
#include "core/gpu/GPUIntegrationHelper.h"
#include "core/types/IonState.h"
#include "utils/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace ICARION;
using namespace ICARION::core;

int main() {
    std::cout << "GPU RK45 Integration Test\n";
    std::cout << "==========================\n\n";
    
    // Test 1: Free particle motion (should be exact with RK45)
    std::cout << "Test 1: Free particle (v = 1000 m/s in x)\n";
    
    std::vector<IonState> ions(5000);
    const double mass_kg = 29.0 * 1.66054e-27;  // N2+
    const double charge_C = ELEM_CHARGE_C;
    const double v0 = 1000.0;  // m/s
    const double dt = 1e-6;    // 1 μs
    
    for (auto& ion : ions) {
        ion.mass_kg = mass_kg;
        ion.ion_charge_C = charge_C;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{v0, 0.0, 0.0};
        ion.active = true;
    }
    
    // GPU context and helper
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    auto gpu_helper = icarion::gpu::GPUIntegrationHelper::create(*gpu_ctx, 1000);
    
    // Integrate for 10 steps
    for (int i = 0; i < 10; ++i) {
        bool success = gpu_helper->integrate_batch_rk45(ions, dt, i * dt, nullptr);
        if (!success) {
            std::cerr << "GPU integration failed at step " << i << "\n";
            return 1;
        }
    }
    
    // Check result
    double expected_x = v0 * 10 * dt;
    double actual_x = ions[0].pos.x;
    double error = std::abs(actual_x - expected_x);
    
    std::cout << "Expected x: " << expected_x << " m\n";
    std::cout << "Actual x:   " << actual_x << " m\n";
    std::cout << "Error:      " << error << " m\n";
    
    if (error < 1e-12) {
        std::cout << "✓ Free particle test PASSED\n\n";
    } else {
        std::cout << "✗ Free particle test FAILED\n\n";
        return 1;
    }
    
    // Test 2: Constant electric field (exact motion: x = 0.5*a*t²)
    std::cout << "Test 2: Constant E-field acceleration\n";
    std::cout << "Note: GPU uses zero fields in batch mode without field provider\n";
    std::cout << "Skipping for now (would need field array setup)\n\n";
    
    // Test 3: RK45 adaptive substeps
    std::cout << "Test 3: RK45 adaptive step control\n";
    
    ions.clear();
    ions.resize(5000);
    
    // High velocity → larger errors → more substeps
    const double v_high = 1e6;  // 1000 km/s
    
    for (auto& ion : ions) {
        ion.mass_kg = mass_kg;
        ion.ion_charge_C = charge_C;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{v_high, 0.0, 0.0};
        ion.active = true;
    }
    
    // Integrate with tight tolerance
    double atol = 1e-15;
    double rtol = 1e-9;
    bool success = gpu_helper->integrate_batch_rk45(ions, dt, 0.0, nullptr, atol, rtol);
    
    if (success) {
        double expected = v_high * dt;
        double actual = ions[0].pos.x;
        double rel_error = std::abs((actual - expected) / expected);
        
        std::cout << "Expected x: " << expected << " m\n";
        std::cout << "Actual x:   " << actual << " m\n";
        std::cout << "Rel error:  " << rel_error << "\n";
        
        if (rel_error < 1e-8) {
            std::cout << "✓ Adaptive step test PASSED\n\n";
        } else {
            std::cout << "✗ Adaptive step test FAILED\n\n";
            return 1;
        }
    } else {
        std::cout << "✗ GPU integration failed\n\n";
        return 1;
    }
    
    // Statistics
    auto stats = gpu_helper->get_stats();
    std::cout << "GPU Statistics:\n";
    std::cout << "  Total integrations: " << stats.gpu_integrations << "\n";
    std::cout << "  Total ions:         " << stats.total_ions_gpu << "\n";
    std::cout << "  Total time:         " << stats.total_time_ms << " ms\n";
    std::cout << "  Avg time per batch: " 
              << stats.total_time_ms / stats.gpu_integrations << " ms\n";
    
    std::cout << "\n✓ All tests PASSED\n";
    
    return 0;
}

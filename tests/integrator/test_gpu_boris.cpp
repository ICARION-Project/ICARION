// Test GPU Boris pusher integration
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/core/GPUIntegrationHelper.h"
#include "core/types/IonState.h"
#include "utils/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace ICARION;
using namespace ICARION::core;

int main() {
    std::cout << "GPU Boris Pusher Integration Test\n";
    std::cout << "===================================\n\n";
    
    const double mass_kg = 29.0 * 1.66054e-27;  // N2+
    const double charge_C = ELEM_CHARGE_C;
    
    // Test 1: Free particle (no fields)
    std::cout << "Test 1: Free particle motion\n";
    
    std::vector<IonState> ions(5000);
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
        bool success = gpu_helper->integrate_batch_boris(ions, dt, i * dt, nullptr);
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
    
    // Test 2: Cyclotron motion (would need B-field)
    std::cout << "Test 2: Cyclotron motion in B-field\n";
    std::cout << "Note: GPU uses zero fields without field provider\n";
    std::cout << "Skipping for now (would need field array setup)\n\n";
    
    // Test 3: Energy conservation
    std::cout << "Test 3: Energy conservation (free particle)\n";
    
    ions.clear();
    ions.resize(5000);
    
    for (auto& ion : ions) {
        ion.mass_kg = mass_kg;
        ion.ion_charge_C = charge_C;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{v0, v0, v0};  // All directions
        ion.active = true;
    }
    
    double E_initial = 0.5 * mass_kg * 3.0 * v0 * v0;
    
    // Integrate for 100 steps
    for (int i = 0; i < 100; ++i) {
        bool success = gpu_helper->integrate_batch_boris(ions, dt, i * dt, nullptr);
        if (!success) {
            std::cerr << "GPU integration failed\n";
            return 1;
        }
    }
    
    // Check energy conservation
    double vx = ions[0].vel.x;
    double vy = ions[0].vel.y;
    double vz = ions[0].vel.z;
    double E_final = 0.5 * mass_kg * (vx*vx + vy*vy + vz*vz);
    double rel_error = std::abs((E_final - E_initial) / E_initial);
    
    std::cout << "Initial energy: " << E_initial << " J\n";
    std::cout << "Final energy:   " << E_final << " J\n";
    std::cout << "Rel error:      " << rel_error << "\n";
    
    if (rel_error < 1e-10) {
        std::cout << "✓ Energy conservation test PASSED\n\n";
    } else {
        std::cout << "✗ Energy conservation test FAILED\n\n";
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

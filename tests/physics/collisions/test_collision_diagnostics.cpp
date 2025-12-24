// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

// Detailed diagnostics for GPU vs CPU collision physics
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/types/IonState.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>
#include <random>
#include <vector>

using namespace ICARION::physics;
using namespace ICARION;
using namespace ICARION::config;

constexpr double EV_TO_JOULES = 1.60218e-19;

// Test 1: Compare neutral velocity sampling
void test_neutral_sampling() {
    std::cout << "\n=== TEST 1: NEUTRAL VELOCITY SAMPLING ===" << std::endl;
    
    double T = 300.0;
    double m_neutral = 4.0 * 1.66054e-27;  // He
    double kB = 1.380649e-23;
    
    // Expected sigma for Maxwell-Boltzmann
    double sigma_expected = std::sqrt(kB * T / m_neutral);
    std::cout << "Expected sigma: " << sigma_expected << " m/s" << std::endl;
    
    // Sample many neutrals from CPU
    std::mt19937_64 rng(42);
    std::normal_distribution<double> dist(0.0, sigma_expected);
    
    std::vector<double> cpu_samples(10000);
    for (auto& v : cpu_samples) {
        v = dist(rng);
    }
    
    // Compute statistics
    double mean_cpu = 0.0, var_cpu = 0.0;
    for (double v : cpu_samples) mean_cpu += v;
    mean_cpu /= cpu_samples.size();
    for (double v : cpu_samples) var_cpu += (v - mean_cpu) * (v - mean_cpu);
    var_cpu /= cpu_samples.size();
    double sigma_cpu = std::sqrt(var_cpu);
    
    std::cout << "CPU samples:" << std::endl;
    std::cout << "  Mean: " << mean_cpu << " m/s (expected 0)" << std::endl;
    std::cout << "  Sigma: " << sigma_cpu << " m/s (expected " << sigma_expected << ")" << std::endl;
    std::cout << "  Ratio: " << sigma_cpu / sigma_expected << std::endl;
    
    // Check thermal energy
    double mean_KE_cpu = 0.0;
    for (double v : cpu_samples) {
        mean_KE_cpu += 0.5 * m_neutral * (v*v + v*v + v*v);  // 3 components
    }
    mean_KE_cpu /= cpu_samples.size();
    double expected_thermal = 1.5 * kB * T;
    
    std::cout << "\nThermal energy check:" << std::endl;
    std::cout << "  Mean KE: " << mean_KE_cpu / EV_TO_JOULES << " eV" << std::endl;
    std::cout << "  Expected: " << expected_thermal / EV_TO_JOULES << " eV" << std::endl;
    std::cout << "  Ratio: " << mean_KE_cpu / expected_thermal << std::endl;
}

// Test 2: Single collision energy conservation
void test_single_collision_energy() {
    std::cout << "\n=== TEST 2: SINGLE COLLISION ENERGY CONSERVATION ===" << std::endl;
    
    double m_ion = 29.0 * 1.66054e-27;
    double m_neutral = 4.0 * 1.66054e-27;
    double T = 300.0;
    double kB = 1.380649e-23;
    
    // Initial state: fast ion, thermal neutral
    Vec3 v_ion_init{1000.0, 0.0, 0.0};  // m/s
    
    std::mt19937_64 rng(12345);
    double sigma = std::sqrt(kB * T / m_neutral);
    std::normal_distribution<double> dist(0.0, sigma);
    
    Vec3 v_neutral{dist(rng), dist(rng), dist(rng)};
    
    std::cout << "Initial state:" << std::endl;
    std::cout << "  Ion velocity: (" << v_ion_init.x << ", " << v_ion_init.y << ", " << v_ion_init.z << ")" << std::endl;
    std::cout << "  Neutral velocity: (" << v_neutral.x << ", " << v_neutral.y << ", " << v_neutral.z << ")" << std::endl;
    
    // Initial energies
    double KE_ion_init = 0.5 * m_ion * (v_ion_init.x*v_ion_init.x + v_ion_init.y*v_ion_init.y + v_ion_init.z*v_ion_init.z);
    double KE_neutral_init = 0.5 * m_neutral * (v_neutral.x*v_neutral.x + v_neutral.y*v_neutral.y + v_neutral.z*v_neutral.z);
    double KE_total_init = KE_ion_init + KE_neutral_init;
    
    std::cout << "\nInitial energies:" << std::endl;
    std::cout << "  Ion KE: " << KE_ion_init / EV_TO_JOULES << " eV" << std::endl;
    std::cout << "  Neutral KE: " << KE_neutral_init / EV_TO_JOULES << " eV" << std::endl;
    std::cout << "  Total KE: " << KE_total_init / EV_TO_JOULES << " eV" << std::endl;
    
    // === STEP BY STEP TRANSFORMATION ===
    
    // 1. Relative velocity
    Vec3 v_rel{v_ion_init.x - v_neutral.x, v_ion_init.y - v_neutral.y, v_ion_init.z - v_neutral.z};
    double v_rel_mag = std::sqrt(v_rel.x*v_rel.x + v_rel.y*v_rel.y + v_rel.z*v_rel.z);
    std::cout << "\nStep 1: Relative velocity" << std::endl;
    std::cout << "  v_rel = (" << v_rel.x << ", " << v_rel.y << ", " << v_rel.z << ")" << std::endl;
    std::cout << "  |v_rel| = " << v_rel_mag << " m/s" << std::endl;
    
    // 2. Center-of-mass velocity
    double m_total = m_ion + m_neutral;
    Vec3 v_cm{
        (m_ion * v_ion_init.x + m_neutral * v_neutral.x) / m_total,
        (m_ion * v_ion_init.y + m_neutral * v_neutral.y) / m_total,
        (m_ion * v_ion_init.z + m_neutral * v_neutral.z) / m_total
    };
    std::cout << "\nStep 2: Center-of-mass velocity" << std::endl;
    std::cout << "  v_cm = (" << v_cm.x << ", " << v_cm.y << ", " << v_cm.z << ")" << std::endl;
    
    // Verify: v_ion - v_cm should equal (m_neutral/m_total) * v_rel
    Vec3 v_ion_cm{v_ion_init.x - v_cm.x, v_ion_init.y - v_cm.y, v_ion_init.z - v_cm.z};
    double expected_factor = m_neutral / m_total;
    std::cout << "  v_ion (CM frame) = (" << v_ion_cm.x << ", " << v_ion_cm.y << ", " << v_ion_cm.z << ")" << std::endl;
    std::cout << "  Expected: v_rel * (m_neutral/m_total) = v_rel * " << expected_factor << std::endl;
    std::cout << "  Check x: " << v_ion_cm.x << " vs " << v_rel.x * expected_factor << std::endl;
    
    // 3. Scatter: flip direction (simple test)
    Vec3 v_rel_scattered{-v_rel.x, -v_rel.y, -v_rel.z};
    std::cout << "\nStep 3: Scattered relative velocity (flipped)" << std::endl;
    std::cout << "  v_rel_scattered = (" << v_rel_scattered.x << ", " << v_rel_scattered.y << ", " << v_rel_scattered.z << ")" << std::endl;
    
    // 4. Transform back to lab frame
    Vec3 v_ion_final{
        v_cm.x + v_rel_scattered.x * expected_factor,
        v_cm.y + v_rel_scattered.y * expected_factor,
        v_cm.z + v_rel_scattered.z * expected_factor
    };
    std::cout << "\nStep 4: Transform back to lab frame" << std::endl;
    std::cout << "  v_ion_final = v_cm + v_rel_scattered * (m_neutral/m_total)" << std::endl;
    std::cout << "  v_ion_final = (" << v_ion_final.x << ", " << v_ion_final.y << ", " << v_ion_final.z << ")" << std::endl;
    
    // What should neutral velocity be?
    Vec3 v_neutral_final{
        v_cm.x - v_rel_scattered.x * (m_ion / m_total),
        v_cm.y - v_rel_scattered.y * (m_ion / m_total),
        v_cm.z - v_rel_scattered.z * (m_ion / m_total)
    };
    std::cout << "  v_neutral_final = (" << v_neutral_final.x << ", " << v_neutral_final.y << ", " << v_neutral_final.z << ")" << std::endl;
    
    // Final energies
    double KE_ion_final = 0.5 * m_ion * (v_ion_final.x*v_ion_final.x + v_ion_final.y*v_ion_final.y + v_ion_final.z*v_ion_final.z);
    double KE_neutral_final = 0.5 * m_neutral * (v_neutral_final.x*v_neutral_final.x + v_neutral_final.y*v_neutral_final.y + v_neutral_final.z*v_neutral_final.z);
    double KE_total_final = KE_ion_final + KE_neutral_final;
    
    std::cout << "\nFinal energies:" << std::endl;
    std::cout << "  Ion KE: " << KE_ion_final / EV_TO_JOULES << " eV" << std::endl;
    std::cout << "  Neutral KE: " << KE_neutral_final / EV_TO_JOULES << " eV" << std::endl;
    std::cout << "  Total KE: " << KE_total_final / EV_TO_JOULES << " eV" << std::endl;
    
    std::cout << "\nEnergy conservation check:" << std::endl;
    std::cout << "  ΔE_total: " << (KE_total_final - KE_total_init) / EV_TO_JOULES << " eV" << std::endl;
    std::cout << "  Conserved: " << (std::abs(KE_total_final - KE_total_init) < 1e-20 ? "YES ✓" : "NO ✗") << std::endl;
}

// Test 3: Compare CPU vs GPU formula application
void test_cpu_vs_gpu_formula() {
    std::cout << "\n=== TEST 3: CPU VS GPU TRANSFORMATION FORMULA ===" << std::endl;
    
    double m_ion = 29.0 * 1.66054e-27;
    double m_neutral = 4.0 * 1.66054e-27;
    double m_total = m_ion + m_neutral;
    
    // Test case: ion moving, neutral at rest
    Vec3 v_ion{1000.0, 0.0, 0.0};
    Vec3 v_neutral{0.0, 0.0, 0.0};
    
    // Relative velocity
    Vec3 v_rel{v_ion.x - v_neutral.x, v_ion.y - v_neutral.y, v_ion.z - v_neutral.z};
    
    // COM velocity
    Vec3 v_cm{
        (m_ion * v_ion.x + m_neutral * v_neutral.x) / m_total,
        (m_ion * v_ion.y + m_neutral * v_neutral.y) / m_total,
        (m_ion * v_ion.z + m_neutral * v_neutral.z) / m_total
    };
    
    std::cout << "Input:" << std::endl;
    std::cout << "  v_ion = (" << v_ion.x << ", 0, 0)" << std::endl;
    std::cout << "  v_neutral = (0, 0, 0)" << std::endl;
    std::cout << "  v_cm = (" << v_cm.x << ", 0, 0)" << std::endl;
    std::cout << "  v_rel = (" << v_rel.x << ", 0, 0)" << std::endl;
    
    // Scatter: flip direction
    Vec3 v_rel_scattered{-v_rel.x, 0.0, 0.0};
    
    // === CPU FORMULA (from CollisionKernels.cpp line 128) ===
    // return Vcom + vrel_scattered * (neutral_mass * inv_mt);
    double inv_mt = 1.0 / m_total;
    Vec3 v_ion_cpu{
        v_cm.x + v_rel_scattered.x * (m_neutral * inv_mt),
        v_cm.y + v_rel_scattered.y * (m_neutral * inv_mt),
        v_cm.z + v_rel_scattered.z * (m_neutral * inv_mt)
    };
    
    // === GPU FORMULA (from collision_kernels_gpu.cu line 226) ===
    // vx = vcm_x + vrel_scattered_x * reduced_mass_factor;
    // where reduced_mass_factor = m_neutral / m_total
    double reduced_mass_factor = m_neutral / m_total;
    Vec3 v_ion_gpu{
        v_cm.x + v_rel_scattered.x * reduced_mass_factor,
        v_cm.y + v_rel_scattered.y * reduced_mass_factor,
        v_cm.z + v_rel_scattered.z * reduced_mass_factor
    };
    
    std::cout << "\nCPU formula:" << std::endl;
    std::cout << "  v_ion_new = v_cm + v_rel_scattered * (m_neutral / m_total)" << std::endl;
    std::cout << "  v_ion_new = (" << v_ion_cpu.x << ", 0, 0)" << std::endl;
    
    std::cout << "\nGPU formula:" << std::endl;
    std::cout << "  v_ion_new = v_cm + v_rel_scattered * (m_neutral / m_total)" << std::endl;
    std::cout << "  v_ion_new = (" << v_ion_gpu.x << ", 0, 0)" << std::endl;
    
    std::cout << "\nDifference:" << std::endl;
    std::cout << "  Δvx = " << (v_ion_gpu.x - v_ion_cpu.x) << " m/s" << std::endl;
    
    if (std::abs(v_ion_gpu.x - v_ion_cpu.x) < 1e-10) {
        std::cout << "  ✓ IDENTICAL!" << std::endl;
    } else {
        std::cout << "  ✗ DIFFERENT!" << std::endl;
    }
    
    // Check energies
    double KE_ion_cpu = 0.5 * m_ion * (v_ion_cpu.x*v_ion_cpu.x);
    double KE_ion_gpu = 0.5 * m_ion * (v_ion_gpu.x*v_ion_gpu.x);
    
    std::cout << "\nKinetic energies:" << std::endl;
    std::cout << "  CPU: " << KE_ion_cpu / EV_TO_JOULES << " eV" << std::endl;
    std::cout << "  GPU: " << KE_ion_gpu / EV_TO_JOULES << " eV" << std::endl;
}

int main() {
    std::cout << std::setprecision(6) << std::fixed;
    
    test_neutral_sampling();
    test_single_collision_energy();
    test_cpu_vs_gpu_formula();
    
    std::cout << "\n=== ALL DIAGNOSTICS COMPLETE ===" << std::endl;
    return 0;
}

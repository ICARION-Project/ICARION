// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

// Quick test to verify HSS physics fix
#ifdef ICARION_USE_GPU

#include <iostream>
#include <cmath>
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/types/CollisionTypes.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"

using namespace ICARION;
using namespace ICARION::config;
using namespace ICARION::physics;

static void run_collision(HSSCollisionHandler& handler,
                          IonState& ion,
                          double dt,
                          PhysicsRng& rng,
                          const EnvironmentConfig& env) {
    auto ens = core::IonEnsemble::from_legacy({ion});
    auto view = ens.collision_data(0);
    handler.handle_collision(view, dt, rng, env);
    ion.vel = view.kin.vel();
}

int main() {
    // Environment
    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_mass_kg = 4.0026 * 1.66054e-27;
    env.particle_density_m_3 = env.pressure_Pa / (1.380649e-23 * env.temperature_K);
    
    // Ion properties
    const double mass = 29.0 * 1.66054e-27;
    const double CCS = 45e-20;
    const double dt = 1e-6;
    
    // Initial velocity (same for both)
    const double vx0 = 500.0, vy0 = 300.0, vz0 = 200.0;
    
    // CPU test - run 10000 ions
    std::cout << "CPU: Processing 10000 ions..." << std::endl;
    std::vector<IonState> cpu_ions(10000);
    for (auto& ion : cpu_ions) {
        ion.vel = Vec3(vx0, vy0, vz0);
        ion.mass_kg = mass;
        ion.CCS_m2 = CCS;
        ion.active = true;
    }
    
    PhysicsRng cpu_rng(42);
    HSSCollisionHandler cpu_handler(false, nullptr);
    for (auto& ion : cpu_ions) {
        run_collision(cpu_handler, ion, dt, cpu_rng, env);
    }
    
    // Compute mean
    double cpu_mean_vx = 0, cpu_mean_vy = 0, cpu_mean_vz = 0;
    for (const auto& ion : cpu_ions) {
        cpu_mean_vx += ion.vel.x;
        cpu_mean_vy += ion.vel.y;
        cpu_mean_vz += ion.vel.z;
    }
    cpu_mean_vx /= cpu_ions.size();
    cpu_mean_vy /= cpu_ions.size();
    cpu_mean_vz /= cpu_ions.size();
    double cpu_mean_speed = std::sqrt(cpu_mean_vx*cpu_mean_vx + cpu_mean_vy*cpu_mean_vy + cpu_mean_vz*cpu_mean_vz);
    
    std::cout << "CPU mean: vx=" << cpu_mean_vx << " vy=" << cpu_mean_vy << " vz=" << cpu_mean_vz << " speed=" << cpu_mean_speed << std::endl;
    
    // GPU test - run 10000 ions
    std::cout << "GPU: Processing 10000 ions..." << std::endl;
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    if (!gpu_ctx) {
        std::cerr << "GPU not available!" << std::endl;
        return 1;
    }
    
    std::vector<IonState> gpu_ions(10000);
    for (auto& ion : gpu_ions) {
        ion.vel = Vec3(vx0, vy0, vz0);
        ion.mass_kg = mass;
        ion.CCS_m2 = CCS;
        ion.active = true;
    }
    
    auto gpu_helper = icarion::gpu::GPUCollisionHelper::create(*gpu_ctx, 5000, "HSS", 42);
    gpu_helper->process_collisions_batch(gpu_ions, dt, env);
    
    // Compute mean
    double gpu_mean_vx = 0, gpu_mean_vy = 0, gpu_mean_vz = 0;
    for (const auto& ion : gpu_ions) {
        gpu_mean_vx += ion.vel.x;
        gpu_mean_vy += ion.vel.y;
        gpu_mean_vz += ion.vel.z;
    }
    gpu_mean_vx /= gpu_ions.size();
    gpu_mean_vy /= gpu_ions.size();
    gpu_mean_vz /= gpu_ions.size();
    double gpu_mean_speed = std::sqrt(gpu_mean_vx*gpu_mean_vx + gpu_mean_vy*gpu_mean_vy + gpu_mean_vz*gpu_mean_vz);
    
    std::cout << "GPU mean: vx=" << gpu_mean_vx << " vy=" << gpu_mean_vy << " vz=" << gpu_mean_vz << " speed=" << gpu_mean_speed << std::endl;
    
    // Compare
    double speed_diff_pct = 100.0 * std::abs(cpu_mean_speed - gpu_mean_speed) / cpu_mean_speed;
    std::cout << "\nSpeed difference: " << speed_diff_pct << "%" << std::endl;
    
    if (speed_diff_pct < 20.0) {
        std::cout << "✓ PASS: CPU and GPU are statistically similar (< 20% difference)" << std::endl;
        return 0;
    } else {
        std::cout << "✗ FAIL: CPU and GPU differ significantly (> 20% difference)" << std::endl;
        return 1;
    }
}

#else
int main() {
    std::cerr << "GPU not enabled!" << std::endl;
    return 1;
}
#endif

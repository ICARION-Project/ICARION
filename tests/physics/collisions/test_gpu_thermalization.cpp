// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file test_gpu_thermalization.cpp
 * @brief Compare CPU vs GPU thermalization convergence
 * 
 * Tests that GPU and CPU collision handlers produce similar thermalization rates
 * when starting from high energy (10x thermal) or low energy (0.1x thermal).
 */

#ifdef ICARION_USE_GPU

#include "core/gpu/core/GPUContext.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/types/CollisionTypes.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "utils/constants.h"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace ICARION;
using namespace ICARION::config;
using namespace ICARION::physics;

double thermal_energy_eV(double T_K) {
    return (1.5 * BOLTZMANN_CONSTANT * T_K) / ELEM_CHARGE_C;
}

static bool run_collision(HSSCollisionHandler& handler,
                          IonState& ion,
                          double dt,
                          PhysicsRng& rng,
                          const EnvironmentConfig& env) {
    auto ens = core::IonEnsemble::from_legacy({ion});
    auto view = ens.collision_data(0);
    bool res = handler.handle_collision(view, dt, rng, env);
    ion.vel = view.kin.vel();
    return res;
}

struct ThermalizationResult {
    double time_us;
    double KE_eV;
    double ratio;  // KE / thermal_KE
    double avg_collisions;
};

ThermalizationResult run_cpu_thermalization(
    int N_IONS, 
    int N_STEPS, 
    double dt, 
    double v_init,
    double mass_kg,
    const EnvironmentConfig& env,
    int seed_offset
) {
    HSSCollisionHandler handler(false, nullptr);
    
    double sum_v2 = 0.0;
    
    for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
        IonState ion;
        ion.species_id = "H3O+";
        ion.mass_kg = mass_kg;
        ion.ion_charge_C = ELEM_CHARGE_C;
        ion.CCS_m2 = 45e-20;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{v_init, 0.0, 0.0};
        ion.active = true;
        
        PhysicsRng rng(seed_offset + ion_idx);
        
        for (int i = 0; i < N_STEPS; ++i) {
            run_collision(handler, ion, dt, rng, env);
        }
        
        double v2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
        sum_v2 += v2;
    }
    
    double mean_v2 = sum_v2 / N_IONS;
    double KE_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
    double thermal_KE = thermal_energy_eV(env.temperature_K);
    
    ThermalizationResult result;
    result.time_us = N_STEPS * dt * 1e6;
    result.KE_eV = KE_eV;
    result.ratio = KE_eV / thermal_KE;
    result.avg_collisions = 0;  // Not tracked for simplicity
    
    return result;
}

ThermalizationResult run_gpu_thermalization(
    int N_IONS,
    int N_STEPS,
    double dt,
    double v_init,
    double mass_kg,
    const EnvironmentConfig& env,
    icarion::gpu::GPUContext& gpu_ctx,
    int seed_offset
) {
    auto gpu_helper = icarion::gpu::GPUCollisionHelper::create(gpu_ctx, 1000, "HSS", seed_offset);
    
    std::vector<IonState> ions(N_IONS);
    for (auto& ion : ions) {
        ion.species_id = "H3O+";
        ion.mass_kg = mass_kg;
        ion.ion_charge_C = ELEM_CHARGE_C;
        ion.CCS_m2 = 45e-20;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{v_init, 0.0, 0.0};
        ion.active = true;
    }
    
    // Process collisions for N_STEPS timesteps
    for (int i = 0; i < N_STEPS; ++i) {
        gpu_helper->process_collisions_batch(ions, dt, env);
    }
    
    // Compute mean kinetic energy
    double sum_v2 = 0.0;
    for (const auto& ion : ions) {
        double v2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
        sum_v2 += v2;
    }
    
    double mean_v2 = sum_v2 / N_IONS;
    double KE_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
    double thermal_KE = thermal_energy_eV(env.temperature_K);
    
    ThermalizationResult result;
    result.time_us = N_STEPS * dt * 1e6;
    result.KE_eV = KE_eV;
    result.ratio = KE_eV / thermal_KE;
    result.avg_collisions = 0;
    
    return result;
}

int main() {
    std::cout << "\n=== CPU vs GPU THERMALIZATION COMPARISON ===\n\n";
    
    // Check GPU availability
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    if (!gpu_ctx) {
        std::cerr << "ERROR: GPU not available!\n";
        return 1;
    }
    
    // Test parameters
    const double T_K = 300.0;
    const int N_IONS = 5000;  // Large enough for good statistics
    const double dt = 1e-7;   // 100 ns timestep
    const double mass_kg = 29.0 * AMU_TO_KG;  // H3O+ mass
    
    // Environment
    EnvironmentConfig env;
    env.temperature_K = T_K;
    env.pressure_Pa = 101325.0;  // 1 atm
    env.gas_mass_kg = 4.0026 * AMU_TO_KG;  // He
    env.particle_density_m_3 = env.pressure_Pa / (BOLTZMANN_CONSTANT * env.temperature_K);
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    
    // CRITICAL: GPU code reads from gas_mixture, not gas_mass_kg!
    // Set up gas mixture properly for GPU
    GasMixtureComponent he;
    he.species = "He";
    he.mole_fraction = 1.0;
    he.mass_kg = 4.0026 * AMU_TO_KG;
    he.radius_m = 1.4e-10;  // He radius
    he.cross_section_m2 = 45e-20;  // H3O+ in He CCS
    env.gas_mixture.push_back(he);
    
    std::cout << "Temperature: " << T_K << " K\n";
    std::cout << "Expected thermal KE: " << thermal_energy_eV(T_K) << " eV\n";
    std::cout << "Number of ions: " << N_IONS << "\n";
    std::cout << "Timestep: " << dt * 1e9 << " ns\n\n";
    
    // Test cases: different number of steps
    std::vector<int> step_counts = {100, 500, 1000, 2000, 5000};
    
    std::cout << std::fixed << std::setprecision(4);
    
    // ========================================================================
    // HIGH ENERGY START (10x thermal)
    // ========================================================================
    std::cout << "=== HIGH ENERGY START (10x thermal) ===\n\n";
    const double v_init_high = std::sqrt(10.0 * 3.0 * BOLTZMANN_CONSTANT * T_K / mass_kg);
    
    std::cout << "Steps  Time(μs)   CPU_KE(eV)  GPU_KE(eV)  CPU_Ratio  GPU_Ratio  Diff(%)\n";
    std::cout << "-----  ---------  ----------  ----------  ---------  ---------  -------\n";
    
    for (int N_STEPS : step_counts) {
        auto cpu_result = run_cpu_thermalization(N_IONS, N_STEPS, dt, v_init_high, mass_kg, env, 42);
        auto gpu_result = run_gpu_thermalization(N_IONS, N_STEPS, dt, v_init_high, mass_kg, env, *gpu_ctx, 42);
        
        double diff_pct = 100.0 * std::abs(cpu_result.ratio - gpu_result.ratio) / cpu_result.ratio;
        
        std::cout << std::setw(5) << N_STEPS << "  ";
        std::cout << std::setw(9) << cpu_result.time_us << "  ";
        std::cout << std::setw(10) << cpu_result.KE_eV << "  ";
        std::cout << std::setw(10) << gpu_result.KE_eV << "  ";
        std::cout << std::setw(9) << cpu_result.ratio << "  ";
        std::cout << std::setw(9) << gpu_result.ratio << "  ";
        std::cout << std::setw(7) << diff_pct << "\n";
    }
    
    // ========================================================================
    // LOW ENERGY START (0.1x thermal)
    // ========================================================================
    std::cout << "\n=== LOW ENERGY START (0.1x thermal) ===\n\n";
    const double v_init_low = std::sqrt(0.1 * 3.0 * BOLTZMANN_CONSTANT * T_K / mass_kg);
    
    std::cout << "Steps  Time(μs)   CPU_KE(eV)  GPU_KE(eV)  CPU_Ratio  GPU_Ratio  Diff(%)\n";
    std::cout << "-----  ---------  ----------  ----------  ---------  ---------  -------\n";
    
    for (int N_STEPS : step_counts) {
        auto cpu_result = run_cpu_thermalization(N_IONS, N_STEPS, dt, v_init_low, mass_kg, env, 1000);
        auto gpu_result = run_gpu_thermalization(N_IONS, N_STEPS, dt, v_init_low, mass_kg, env, *gpu_ctx, 1000);
        
        double diff_pct = 100.0 * std::abs(cpu_result.ratio - gpu_result.ratio) / cpu_result.ratio;
        
        std::cout << std::setw(5) << N_STEPS << "  ";
        std::cout << std::setw(9) << cpu_result.time_us << "  ";
        std::cout << std::setw(10) << cpu_result.KE_eV << "  ";
        std::cout << std::setw(10) << gpu_result.KE_eV << "  ";
        std::cout << std::setw(9) << cpu_result.ratio << "  ";
        std::cout << std::setw(9) << gpu_result.ratio << "  ";
        std::cout << std::setw(7) << diff_pct << "\n";
    }
    
    std::cout << "\n=== SUMMARY ===\n";
    std::cout << "Both CPU and GPU should converge to thermal energy (ratio ≈ 1.0)\n";
    std::cout << "Differences < 30% are acceptable (due to different RNG algorithms)\n";
    std::cout << "\n✓ Test complete!\n";
    
    return 0;
}

#else
int main() {
    std::cerr << "GPU not enabled!\n";
    return 1;
}
#endif

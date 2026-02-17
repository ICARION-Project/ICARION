// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#include <iostream>

#ifdef ICARION_USE_GPU
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#include "core/types/IonState.h"
#include "core/config/types/EnvironmentConfig.h"
#include <iomanip>
#include <cmath>
#include <vector>

using namespace ICARION;
using namespace ICARION::config;
#endif

int main() {
#ifndef ICARION_USE_GPU
    std::cout << "GPU build disabled; skipping GPU-only thermalization test.\n";
    return 0;
#else
    const int N_IONS = 5000;
    const double mass_kg = 29.0 * 1.66054e-27;
    const double T = 300.0;
    const double kB = 1.380649e-23;
    const double eV = 1.60218e-19;
    
    EnvironmentConfig env;
    env.temperature_K = T;
    env.pressure_Pa = 101325.0;
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    
    GasMixtureComponent he;
    he.species = "He";
    he.mole_fraction = 1.0;
    he.mass_kg = 4.0026 * 1.66054e-27;
    he.radius_m = 1.4e-10;
    env.gas_mixture.push_back(he);
    
    double v_thermal = std::sqrt(3 * kB * T / mass_kg);
    double v_init = v_thermal * std::sqrt(10.0);
    
    std::vector<IonState> ions(N_IONS);
    for (auto& ion : ions) {
        ion.species_id = "H3O+";
        ion.mass_kg = mass_kg;
        ion.CCS_m2 = 45e-20;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{v_init, 0.0, 0.0};
        ion.active = true;
    }
    
    auto gpu_ctx = icarion::gpu::GPUContext::create(0);
    auto gpu_helper = icarion::gpu::GPUCollisionHelper::create(*gpu_ctx, 1000, "HSS", 42);
    
    double dt = 1e-7;
    std::vector<int> steps = {100, 500, 1000, 2000, 5000};
    
    std::cout << std::fixed << std::setprecision(4);
    std::cout << "\nGPU Thermalization (He gas)\n";
    std::cout << "Expected thermal: " << 1.5 * kB * T / eV << " eV\n\n";
    std::cout << "Steps    KE(eV)   Ratio\n";
    std::cout << "-----  --------  ------\n";
    
    for (int N : steps) {
        for (auto& ion : ions) ion.vel = Vec3{v_init, 0.0, 0.0};
        
        for (int i = 0; i < N; ++i) {
            gpu_helper->process_collisions_batch(ions, dt, env);
        }
        
        double sum_v2 = 0.0;
        for (const auto& ion : ions) {
            sum_v2 += ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
        }
        double KE = 0.5 * mass_kg * (sum_v2 / N_IONS) / eV;
        double ratio = KE / (1.5 * kB * T / eV);
        
        std::cout << std::setw(5) << N << std::setw(10) << KE << std::setw(8) << ratio << "\n";
    }
    
    return 0;
#endif
}

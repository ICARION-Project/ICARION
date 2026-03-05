// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifdef ICARION_USE_GPU

#include "core/physics/collisions/EHSSCollisionHandler.h"
#include "core/gpu/core/GPUContext.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "utils/constants.h"
#include <vector>
#include <cmath>
#include <stdexcept>
#include <iostream>

using namespace ICARION;
using namespace ICARION::config;
using namespace ICARION::physics;

// Simple thermalization check: CPU vs GPU EHSS should both move toward thermal KE
int main() {
    auto ctx = icarion::gpu::GPUContext::create(0);
    if (!ctx) {
        std::cerr << "No GPU context available; skipping\n";
        return 0;
    }

    // Geometry: single-sphere placeholder (acts close to HSS)
    GeometryMap geom;
    std::vector<Vec3> pos = {Vec3{0.0, 0.0, 0.0}};
    std::vector<double> radii = {1.5e-10};
    geom["A"] = {pos, radii};

    // Environment (He buffer)
    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_mass_kg = MOLAR_MASS_HE_KG;
    env.gas_radius_m = 1.4e-10;
    env.particle_density_m_3 = env.pressure_Pa / (BOLTZMANN_CONSTANT * env.temperature_K);

    // Ion ensemble
    const size_t N = 4000;
    const double mass_kg = 20.0 * AMU_TO_KG;
    std::vector<IonState> ions_cpu(N), ions_gpu(N);
    const double kB = BOLTZMANN_CONSTANT;
    const double thermal_v = std::sqrt(3.0 * kB * env.temperature_K / mass_kg);
    const double v_init = thermal_v * std::sqrt(10.0);  // 10x thermal KE
    for (size_t i = 0; i < N; ++i) {
        IonState ion;
        ion.species_id = "A";
        ion.mass_kg = mass_kg;
        ion.CCS_m2 = M_PI * radii[0] * radii[0];
        ion.active = true;
        ion.born = true;
        ion.vel = Vec3{v_init, 0.0, 0.0};
        ions_cpu[i] = ion;
        ions_gpu[i] = ion;
    }

    // CPU handler
    EHSSCollisionHandler cpu_handler(geom, false, nullptr);
    PhysicsRng rng_cpu(1234);

    // GPU helper
    auto gpu_helper = icarion::gpu::GPUCollisionHelper::create(*ctx, 1, "EHSS", 1234);
    gpu_helper->set_geometry(geom);

    const double dt = 1e-8;
    const int steps = 500;

    core::IonEnsemble ensemble_cpu = core::IonEnsemble::from_legacy(ions_cpu);
    // CPU collisions
    for (int s = 0; s < steps; ++s) {
        for (size_t i = 0; i < N; ++i) {
            auto view = ensemble_cpu.collision_data(i);
            cpu_handler.handle_collision(view, dt, rng_cpu, env);
        }
    }
    ions_cpu = ensemble_cpu.to_legacy();

    // GPU collisions
    for (int s = 0; s < steps; ++s) {
        gpu_helper->process_collisions_batch(ions_gpu, dt, env);
    }

    auto mean_ke_ratio = [&](const std::vector<IonState>& ions) {
        double sum = 0.0;
        for (const auto& ion : ions) {
            double v2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
            sum += 0.5 * ion.mass_kg * v2;
        }
        double ke_mean = sum / static_cast<double>(ions.size());
        double ke_thermal = 1.5 * kB * env.temperature_K;
        return ke_mean / ke_thermal;
    };

    double cpu_ratio = mean_ke_ratio(ions_cpu);
    double gpu_ratio = mean_ke_ratio(ions_gpu);

    if (!(cpu_ratio > 0.3 && cpu_ratio < 2.5)) {
        throw std::runtime_error("CPU EHSS failed to thermalize: ratio=" + std::to_string(cpu_ratio));
    }
    if (!(gpu_ratio > 0.3 && gpu_ratio < 2.5)) {
        throw std::runtime_error("GPU EHSS failed to thermalize: ratio=" + std::to_string(gpu_ratio));
    }

    // Parity: expect similar thermalization level
    if (std::fabs(cpu_ratio - gpu_ratio) > 0.3) {
        throw std::runtime_error("GPU EHSS thermalization deviates from CPU (cpu=" +
                                 std::to_string(cpu_ratio) + ", gpu=" + std::to_string(gpu_ratio) + ")");
    }

    return 0;
}

#else
int main() { return 0; }
#endif

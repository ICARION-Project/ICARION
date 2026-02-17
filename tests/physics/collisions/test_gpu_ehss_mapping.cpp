// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
#ifdef ICARION_USE_GPU

#include "core/gpu/core/GPUContext.h"
#include "core/gpu/collisions/GPUCollisionHelper.h"
#include "core/types/IonState.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"
#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>

using namespace ICARION;
using namespace ICARION::config;

int main() {
    auto ctx = icarion::gpu::GPUContext::create(0);
    if (!ctx) {
        std::cerr << "No GPU context available\n";
        return 0; // Skip silently if no GPU
    }

    auto helper = icarion::gpu::GPUCollisionHelper::create(*ctx, 1, "EHSS", 123);
    if (!helper) {
        throw std::runtime_error("Failed to create GPUCollisionHelper");
    }

    // Two-species geometry map: A and B
    icarion::gpu::GPUCollisionHelper::GeometryMap geom;
    std::vector<Vec3> positions = {Vec3{0.0, 0.0, 0.0}};
    std::vector<double> radii = {1e-10};
    geom["A"] = {positions, radii};
    geom["B"] = {positions, radii};
    helper->set_geometry(geom);

    // Environment (single component)
    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_mixture.push_back({});
    auto& comp = env.gas_mixture.back();
    comp.species = "He";
    comp.mole_fraction = 1.0;
    comp.density_m3 = 2.5e25;
    comp.mass_kg = MOLAR_MASS_HE_KG;
    comp.radius_m = 1.4e-10;
    comp.participates_in_collisions = true;

    // Ions belonging to species B to exercise non-zero species index
    std::vector<IonState> ions(4);
    for (auto& ion : ions) {
        ion.species_id = "B";
        ion.mass_kg = 20.0 * AMU_TO_KG;
        ion.CCS_m2 = 20e-20;
        ion.active = true;
        ion.vel = Vec3{100.0, 0.0, 0.0};
    }

    const double dt = 1e-7;
    bool ok = helper->process_collisions_batch(ions, dt, env);
    if (!ok) {
        throw std::runtime_error("GPU EHSS batch processing failed");
    }

    // At least one ion velocity should change if collisions executed with correct mapping
    bool changed = false;
    for (const auto& ion : ions) {
        if (std::fabs(ion.vel.x - 100.0) > 1e-6 || std::fabs(ion.vel.y) > 1e-6 || std::fabs(ion.vel.z) > 1e-6) {
            changed = true;
            break;
        }
    }
    if (!changed) {
        throw std::runtime_error("GPU EHSS mapping did not modify velocities (species index mapping may be wrong)");
    }

    return 0;
}

#else
int main() { return 0; }
#endif

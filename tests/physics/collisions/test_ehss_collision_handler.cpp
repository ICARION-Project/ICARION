// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file test_ehss_collision_handler.cpp
 * @brief Thermalization test for EHSSCollisionHandler
 * 
 * Tests that EHSS collision handler correctly thermalizes ions to
 * the environment temperature over many collision steps.
 */

#include "core/physics/collisions/EHSSCollisionHandler.h"
#include "core/physics/collisions/geometryUtils.h"  // Phase 2E: SSOT geometry loading
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

using namespace ICARION::physics;
using namespace ICARION::config;
using namespace ICARION::core;
using Catch::Approx;

static bool run_collision(EHSSCollisionHandler& handler,
                          IonState& ion,
                          double dt,
                          PhysicsRng& rng,
                          const EnvironmentConfig& env) {
    auto ens = IonEnsemble::from_legacy({ion});
    auto view = ens.collision_data(0);
    bool res = handler.handle_collision(view, dt, rng, env);
    ion.vel = view.kin.vel();
    return res;
}

// Helper: Calculate kinetic energy from velocity
double kinetic_energy_eV(const Vec3& vel, double mass_kg) {
    double v2 = vel.x * vel.x + vel.y * vel.y + vel.z * vel.z;
    double KE_J = 0.5 * mass_kg * v2;
    return KE_J / ELEM_CHARGE_C;  // Convert to eV
}

// Helper: Calculate thermal energy (3/2 * kB * T)
double thermal_energy_eV(double T_K) {
    return (1.5 * BOLTZMANN_CONSTANT * T_K) / ELEM_CHARGE_C;
}

// Helper: Load H3O+ geometry using SSOT utility (Phase 2E)
GeometryMap load_h3o_geometry() {
    // SSOT: Use centralized geometry loading
    // Path relative to build/tests/physics/collisions/
    return ICARION::physics::load_geometry_map({"H3O+"}, "../../../../data/molecules");
}

TEST_CASE("EHSSCollisionHandler: Thermalization of H3O+", "[collision][ehss][thermalization]") {
    // ================================================================
    // Setup: 1000 H3O+ ions in He buffer gas at 300 K
    // ================================================================
    
    const int N_IONS = 1000;
    const double T_K = 300.0;
    
    // Environment: He at 300 K, 10000 Pa (high pressure for fast thermalization)
    EnvironmentConfig env;
    env.temperature_K = T_K;
    env.pressure_Pa = 2000.0;
    env.gas_species = "N2";
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.compute_derived_properties();
    
    // Load realistic H3O+ geometry from DFT-optimized structure
    GeometryMap geometry = load_h3o_geometry();
    
    // Create handler
    EHSSCollisionHandler handler(geometry, false);
    
    // ================================================================
    // Simulate 1000 ions with 10,000 steps each
    // ================================================================
    
    const int N_STEPS = 5000;
    const double dt = 1e-8;
    
    double sum_vx2 = 0.0;
    double sum_vy2 = 0.0;
    double sum_vz2 = 0.0;
    int total_collisions = 0;
    
    #pragma omp parallel for reduction(+:sum_vx2,sum_vy2,sum_vz2,total_collisions)
    for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
        // Create ion: H3O+ (mass = 19 amu)
        IonState ion;
        ion.species_id = "H3O+";
        const double mass_kg = 19.0 * AMU_TO_KG;
        ion.mass_kg = mass_kg;
        ion.ion_charge_C = ELEM_CHARGE_C;
        ion.CCS_m2 = 24.9 * ANGSTROM2_TO_M2;  // From DFT calculation
        ion.pos = Vec3{0.0, 0.0, 0.0};
        
        // Start with small initial velocity (triggers collisions)
        const double v_init = 100.0;  // m/s
        ion.vel = Vec3{v_init, 0.0, 0.0};
        
        // Use different RNG seed for each ion
        PhysicsRng rng(42 + ion_idx);
        
        int collision_count = 0;
        for (int i = 0; i < N_STEPS; ++i) {
            bool collided = run_collision(handler, ion, dt, rng, env);
            if (collided) collision_count++;
        }
        
        sum_vx2 += ion.vel.x * ion.vel.x;
        sum_vy2 += ion.vel.y * ion.vel.y;
        sum_vz2 += ion.vel.z * ion.vel.z;
        total_collisions += collision_count;
    }
    
    // ================================================================
    // Verify ensemble-averaged thermalization
    // ================================================================
    
    const double mass_kg = 19.0 * AMU_TO_KG;
    double mean_v2 = (sum_vx2 + sum_vy2 + sum_vz2) / N_IONS;
    double KE_avg_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
    double KE_thermal_eV = thermal_energy_eV(T_K);
    double collision_rate = (double)total_collisions / (N_IONS * N_STEPS);
    
    INFO("Average KE: " << KE_avg_eV << " eV");
    INFO("Thermal KE: " << KE_thermal_eV << " eV");
    INFO("Collision rate: " << collision_rate * 100 << "%");
    INFO("Total collisions: " << total_collisions);
    
    // Ensemble average should be much closer to thermal energy
    REQUIRE(KE_avg_eV > 0.0);
    REQUIRE(KE_avg_eV == Approx(KE_thermal_eV).margin(0.1 * KE_thermal_eV));
    REQUIRE(collision_rate > 0.5);  // Most steps should have collisions
}

TEST_CASE("EHSSCollisionHandler: Thermalization from high energy", "[collision][ehss][thermalization]") {
    // ================================================================
    // Setup: 1000 H3O+ ions starting with 10x thermal energy
    // ================================================================
    
    const int N_IONS = 1000;
    const double T_K = 300.0;
    const double mass_kg = 19.0 * AMU_TO_KG;
    const double v_thermal = std::sqrt(3.0 * BOLTZMANN_CONSTANT * T_K / mass_kg);
    
    EnvironmentConfig env;
    env.temperature_K = T_K;
    env.pressure_Pa = 2000.0;
    env.gas_species = "N2";
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.compute_derived_properties();
    
    // Load realistic H3O+ geometry
    GeometryMap geometry = load_h3o_geometry();
    EHSSCollisionHandler handler(geometry, false);
    
    // ================================================================
    // Simulate 1000 ions starting with high energy
    // ================================================================
    
    const int N_STEPS = 5000;
    const double dt = 1e-8;
    
    double sum_v2_initial = 0.0;
    double sum_v2_final = 0.0;
    
    #pragma omp parallel for reduction(+:sum_v2_initial,sum_v2_final)
    for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
        IonState ion;
        ion.species_id = "H3O+";
        ion.mass_kg = mass_kg;
        ion.ion_charge_C = ELEM_CHARGE_C;
        ion.CCS_m2 = 24.9 * ANGSTROM2_TO_M2;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{10.0 * v_thermal, 0.0, 0.0};  // 10x thermal in x-direction
        // domain_gas_velocity_m_s removed from IonState (passed via EnvironmentConfig)
        
        double v_init2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
        sum_v2_initial += v_init2;
        
        PhysicsRng rng(123 + ion_idx);
        
        for (int i = 0; i < N_STEPS; ++i) {
            run_collision(handler, ion, dt, rng, env);
        }
        
        double v_final2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
        sum_v2_final += v_final2;
    }
    
    // ================================================================
    // Verify ensemble-averaged energy decreased to thermal level
    // ================================================================
    
    double mean_v2_initial = sum_v2_initial / N_IONS;
    double mean_v2_final = sum_v2_final / N_IONS;
    double KE_initial_avg_eV = 0.5 * mass_kg * mean_v2_initial / ELEM_CHARGE_C;
    double KE_final_avg_eV = 0.5 * mass_kg * mean_v2_final / ELEM_CHARGE_C;
    double KE_thermal_eV = thermal_energy_eV(T_K);
    
    INFO("Initial KE (avg): " << KE_initial_avg_eV << " eV");
    INFO("Final KE (avg): " << KE_final_avg_eV << " eV");
    INFO("Thermal KE: " << KE_thermal_eV << " eV");
    
    REQUIRE(KE_final_avg_eV < KE_initial_avg_eV);  // Energy decreased
    REQUIRE(KE_final_avg_eV == Approx(KE_thermal_eV).margin(0.1 * KE_thermal_eV));
}

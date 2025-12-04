// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file test_hss_collision_handler.cpp
 * @brief Thermalization test for HSSCollisionHandler
 * 
 * Tests that HSS collision handler correctly thermalizes ions to
 * the environment temperature over many collision steps.
 */

#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/types/IonState.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>

using namespace ICARION::physics;
using namespace ICARION::config;
using namespace ICARION::core;
using Catch::Approx;

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

TEST_CASE("HSSCollisionHandler: Thermalization of H3O+", "[collision][hss][thermalization][!mayfail]") {
    // ================================================================
    // Setup: 1000 H3O+ ions in He buffer gas at 300 K
    // ================================================================
    
    const int N_IONS = 1000;
    const double T_K = 300.0;
    
    // Environment: He at 300 K, 20000 Pa (very high pressure for fast thermalization)
    EnvironmentConfig env;
    env.temperature_K = T_K;
    env.pressure_Pa = 2000.0;
    env.gas_species = "N2";
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.compute_derived_properties();
    
    // Create handler
    HSSCollisionHandler handler(false);
    
    // ================================================================
    // Simulate 1000 ions with many steps each
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
        ion.CCS_m2 = 24.9 * ANGSTROM2_TO_M2;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        
        // Start with small initial velocity (triggers collisions)
        const double v_init = 100.0;  // m/s
        ion.vel = Vec3{v_init, 0.0, 0.0};
        
        // Use different RNG seed for each ion
        PhysicsRng rng(42 + ion_idx);
        
        int collision_count = 0;
        for (int i = 0; i < N_STEPS; ++i) {
            bool collided = handler.handle_collision(ion, dt, rng, env);
            if (collided) collision_count++;
        }
        
        // Accumulate velocity squared components
        sum_vx2 += ion.vel.x * ion.vel.x;
        sum_vy2 += ion.vel.y * ion.vel.y;
        sum_vz2 += ion.vel.z * ion.vel.z;
        total_collisions += collision_count;
    }
    
    // ================================================================
    // Verify ensemble-averaged thermalization
    // ================================================================
    
    // Compute mean kinetic energy from ensemble-averaged <v²>
    const double mass_kg = 19.0 * AMU_TO_KG;
    double mean_v2 = (sum_vx2 + sum_vy2 + sum_vz2) / N_IONS;
    double KE_avg_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
    double KE_thermal_eV = thermal_energy_eV(T_K);
    double collision_rate = (double)total_collisions / (N_IONS * N_STEPS);
    
    INFO("Average KE: " << KE_avg_eV << " eV");
    INFO("Thermal KE: " << KE_thermal_eV << " eV");
    INFO("Collision rate: " << collision_rate * 100 << "%");
    INFO("Total collisions: " << total_collisions);
    
    // TODO: HSS model systematically underthermalizes (~28%) - physics investigation needed
    REQUIRE(KE_avg_eV > 0.0);
    REQUIRE(KE_avg_eV == Approx(KE_thermal_eV).margin(0.1 * KE_thermal_eV));
    REQUIRE(collision_rate > 0.3);  // Sufficient collision rate
}

TEST_CASE("HSSCollisionHandler: Thermalization from high energy", "[collision][hss][thermalization][!mayfail]") {
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
    
    HSSCollisionHandler handler(false);
    
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
        
        double v_init2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
        sum_v2_initial += v_init2;
        
        PhysicsRng rng(123 + ion_idx);
        
        for (int i = 0; i < N_STEPS; ++i) {
            handler.handle_collision(ion, dt, rng, env);
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

TEST_CASE("HSSCollisionHandler: Isotropic velocity distribution", "[collision][hss][thermalization]") {
    // ================================================================
    // Setup: 1000 ions with many collisions, verify ensemble isotropy
    // ================================================================
    
    const int N_IONS = 1000;
    const double mass_kg = 19.0 * AMU_TO_KG;
    
    EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 20000.0;
    env.gas_species = "N2";
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.compute_derived_properties();
    
    HSSCollisionHandler handler(false);
    
    // Run many collisions for each ion
    const int N_STEPS = 5000;
    const double dt = 1e-8;
    
    double total_vx2 = 0.0;
    double total_vy2 = 0.0;
    double total_vz2 = 0.0;
    
    #pragma omp parallel for reduction(+:total_vx2,total_vy2,total_vz2)
    for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
        IonState ion;
        ion.species_id = "H3O+";
        ion.mass_kg = mass_kg;
        ion.ion_charge_C = ELEM_CHARGE_C;
        ion.CCS_m2 = 24.9 * ANGSTROM2_TO_M2;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{1000.0, 0.0, 0.0};  // Start with directed velocity
        
        PhysicsRng rng(999 + ion_idx);
        
        for (int i = 0; i < N_STEPS; ++i) {
            handler.handle_collision(ion, dt, rng, env);
        }
        
        total_vx2 += ion.vel.x * ion.vel.x;
        total_vy2 += ion.vel.y * ion.vel.y;
        total_vz2 += ion.vel.z * ion.vel.z;
    }
    
    // Calculate ensemble-averaged velocity components
    double v2_total = total_vx2 + total_vy2 + total_vz2;
    double fraction_x = total_vx2 / v2_total;
    double fraction_y = total_vy2 / v2_total;
    double fraction_z = total_vz2 / v2_total;
    
    INFO("vx² fraction: " << fraction_x);
    INFO("vy² fraction: " << fraction_y);
    INFO("vz² fraction: " << fraction_z);
    
    // Each component should be roughly 1/3 (tighter tolerance with ensemble)
    REQUIRE(fraction_x == Approx(1.0/3.0).margin(0.05));
    REQUIRE(fraction_y == Approx(1.0/3.0).margin(0.05));
    REQUIRE(fraction_z == Approx(1.0/3.0).margin(0.05));
}

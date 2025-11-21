// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file test_thermalization_convergence.cpp
 * @brief Test if thermalization converges given enough time
 * 
 * This test tracks the time evolution of ion kinetic energy to check:
 * 1. Is equilibrium reached or is convergence incomplete?
 * 2. How many collisions are needed for equilibrium at different T?
 */

#include "core/physics/collisions/EHSSCollisionHandler.h"
#include "core/physics/collisions/geometryUtils.h"
#include "core/types/IonState.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"
#include <iostream>
#include <iomanip>
#include <vector>

using namespace ICARION::physics;
using namespace ICARION::config;
using namespace ICARION::core;

double thermal_energy_eV(double T_K) {
    return (1.5 * BOLTZMANN_CONSTANT * T_K) / ELEM_CHARGE_C;
}

GeometryMap load_h3o_geometry() {
    return ICARION::physics::load_geometry_map({"H3O+"}, "../../../../data/molecules");
}

int main() {
    std::cout << "\n=== THERMALIZATION CONVERGENCE TEST ===\n\n";
    
    const double T_K = 300.0;
    const int N_IONS = 1000;
    const double dt = 1e-8;
    const double mass_kg = 19.0 * AMU_TO_KG;
    
    // Test with different numbers of steps
    std::vector<int> step_counts = {1000, 2000, 5000, 10000, 20000, 50000};
    
    EnvironmentConfig env;
    env.temperature_K = T_K;
    env.pressure_Pa = 5000.0;
    env.gas_species = "N2";
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.compute_derived_properties();
    
    GeometryMap geometry = load_h3o_geometry();
    EHSSCollisionHandler handler(geometry, false);
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Temperature: " << T_K << " K\n";
    std::cout << "Expected thermal KE: " << thermal_energy_eV(T_K) << " eV\n";
    std::cout << "Starting from HIGH energy (10x thermal)\n\n";
    
    std::cout << "Steps    Time(μs)  KE(eV)    Ratio    Collisions\n";
    std::cout << "-------  --------  --------  -------  ----------\n";
    
    for (int N_STEPS : step_counts) {
        const double v_init = std::sqrt(10.0 * 3.0 * BOLTZMANN_CONSTANT * T_K / mass_kg);
        
        double sum_v2 = 0.0;
        long total_collisions = 0;
        
        #pragma omp parallel for reduction(+:sum_v2,total_collisions)
        for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
            IonState ion;
            ion.species_id = "H3O+";
            ion.mass_kg = mass_kg;
            ion.ion_charge_C = ELEM_CHARGE_C;
            ion.CCS_m2 = 24.9 * ANGSTROM2_TO_M2;
            ion.pos = Vec3{0.0, 0.0, 0.0};
            ion.vel = Vec3{v_init, 0.0, 0.0};
            
            EhssRng rng(42 + ion_idx);
            
            int collision_count = 0;
            for (int i = 0; i < N_STEPS; ++i) {
                bool collided = handler.handle_collision(ion, dt, rng, env);
                if (collided) collision_count++;
            }
            
            double v2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
            sum_v2 += v2;
            total_collisions += collision_count;
        }
        
        double mean_v2 = sum_v2 / N_IONS;
        double KE_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
        double ratio = KE_eV / thermal_energy_eV(T_K);
        double time_us = N_STEPS * dt * 1e6;
        double avg_collisions = (double)total_collisions / N_IONS;
        
        std::cout << std::setw(7) << N_STEPS << "  ";
        std::cout << std::setw(8) << time_us << "  ";
        std::cout << std::setw(8) << KE_eV << "  ";
        std::cout << std::setw(7) << ratio << "  ";
        std::cout << std::setw(10) << avg_collisions << "\n";
    }
    
    std::cout << "\n=== LOW ENERGY START ===\n\n";
    std::cout << "Starting from LOW energy (0.1x thermal)\n\n";
    std::cout << "Steps    Time(μs)  KE(eV)    Ratio    Collisions\n";
    std::cout << "-------  --------  --------  -------  ----------\n";
    
    for (int N_STEPS : step_counts) {
        const double v_init = std::sqrt(0.1 * 3.0 * BOLTZMANN_CONSTANT * T_K / mass_kg);
        
        double sum_v2 = 0.0;
        long total_collisions = 0;
        
        #pragma omp parallel for reduction(+:sum_v2,total_collisions)
        for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
            IonState ion;
            ion.species_id = "H3O+";
            ion.mass_kg = mass_kg;
            ion.ion_charge_C = ELEM_CHARGE_C;
            ion.CCS_m2 = 24.9 * ANGSTROM2_TO_M2;
            ion.pos = Vec3{0.0, 0.0, 0.0};
            ion.vel = Vec3{v_init, 0.0, 0.0};
            
            EhssRng rng(1000 + ion_idx);
            
            int collision_count = 0;
            for (int i = 0; i < N_STEPS; ++i) {
                bool collided = handler.handle_collision(ion, dt, rng, env);
                if (collided) collision_count++;
            }
            
            double v2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
            sum_v2 += v2;
            total_collisions += collision_count;
        }
        
        double mean_v2 = sum_v2 / N_IONS;
        double KE_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
        double ratio = KE_eV / thermal_energy_eV(T_K);
        double time_us = N_STEPS * dt * 1e6;
        double avg_collisions = (double)total_collisions / N_IONS;
        
        std::cout << std::setw(7) << N_STEPS << "  ";
        std::cout << std::setw(8) << time_us << "  ";
        std::cout << std::setw(8) << KE_eV << "  ";
        std::cout << std::setw(7) << ratio << "  ";
        std::cout << std::setw(10) << avg_collisions << "\n";
    }
    
    std::cout << "\n=== COMPARISON AT 600K (high T) ===\n\n";
    
    const double T_high = 600.0;
    env.temperature_K = T_high;
    env.compute_derived_properties();
    
    std::cout << "Temperature: " << T_high << " K\n";
    std::cout << "Expected thermal KE: " << thermal_energy_eV(T_high) << " eV\n\n";
    std::cout << "Steps    Time(μs)  KE(eV)    Ratio\n";
    std::cout << "-------  --------  --------  -------\n";
    
    for (int N_STEPS : step_counts) {
        const double v_init = std::sqrt(10.0 * 3.0 * BOLTZMANN_CONSTANT * T_high / mass_kg);
        
        double sum_v2 = 0.0;
        
        #pragma omp parallel for reduction(+:sum_v2)
        for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
            IonState ion;
            ion.species_id = "H3O+";
            ion.mass_kg = mass_kg;
            ion.ion_charge_C = ELEM_CHARGE_C;
            ion.CCS_m2 = 24.9 * ANGSTROM2_TO_M2;
            ion.pos = Vec3{0.0, 0.0, 0.0};
            ion.vel = Vec3{v_init, 0.0, 0.0};
            
            EhssRng rng(2000 + ion_idx);
            
            for (int i = 0; i < N_STEPS; ++i) {
                handler.handle_collision(ion, dt, rng, env);
            }
            
            double v2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
            sum_v2 += v2;
        }
        
        double mean_v2 = sum_v2 / N_IONS;
        double KE_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
        double ratio = KE_eV / thermal_energy_eV(T_high);
        double time_us = N_STEPS * dt * 1e6;
        
        std::cout << std::setw(7) << N_STEPS << "  ";
        std::cout << std::setw(8) << time_us << "  ";
        std::cout << std::setw(8) << KE_eV << "  ";
        std::cout << std::setw(7) << ratio << "\n";
    }
    
    std::cout << "\n✓ If ratio increases with more steps → convergence incomplete\n";
    std::cout << "✓ If ratio plateaus → equilibrium reached (but wrong value!)\n";
    
    return 0;
}

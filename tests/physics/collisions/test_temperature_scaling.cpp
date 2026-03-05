// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file test_temperature_scaling.cpp
 * @brief Test if equilibrium ion temperature scales correctly with gas temperature
 * 
 * This test checks if the thermalization bug is a constant factor (0.74x)
 * or if the neutral sampling is fundamentally broken.
 * 
 * Expected: Final ion KE should scale linearly with gas temperature T
 * If we get 74% at 300K, we should get 74% at all temperatures.
 */

#include "core/physics/collisions/EHSSCollisionHandler.h"
#include "core/physics/collisions/geometryUtils.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"
#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <cmath>
#include <iostream>
#include <iomanip>

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
    return KE_J / ELEM_CHARGE_C;
}

// Helper: Calculate thermal energy (3/2 * kB * T)
double thermal_energy_eV(double T_K) {
    return (1.5 * BOLTZMANN_CONSTANT * T_K) / ELEM_CHARGE_C;
}

// Helper: Load H3O+ geometry
GeometryMap load_h3o_geometry() {
    return ICARION::physics::load_geometry_map({"H3O+"}, "../../../../data/molecules");
}

TEST_CASE("Temperature scaling: Does equilibrium scale with T_gas?", "[collision][thermalization][temperature]") {
    
    std::cout << "\n=== TEMPERATURE SCALING TEST ===\n";
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "Testing if final ion KE scales linearly with gas temperature\n";
    std::cout << "If bug is constant factor: ratio should be same at all T\n\n";
    
    // Test at multiple temperatures
    std::vector<double> temperatures = {150.0, 300.0, 450.0, 600.0};
    std::vector<double> ratios;
    
    const int N_IONS = 500;
    const int N_STEPS = 5000;
    const double dt = 1e-8;
    const double mass_kg = 19.0 * AMU_TO_KG;
    
    GeometryMap geometry = load_h3o_geometry();
    
    for (double T_K : temperatures) {
        // Setup environment
        EnvironmentConfig env;
        env.temperature_K = T_K;
        env.pressure_Pa = 5000.0;  // Moderate pressure for fast equilibration
        env.gas_species = "N2";
        env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
        env.compute_derived_properties();
        
        EHSSCollisionHandler handler(geometry, false);
        
        // Start ions with high energy (10x thermal)
        const double v_init = std::sqrt(10.0 * 3.0 * BOLTZMANN_CONSTANT * T_K / mass_kg);
        
        double sum_v2 = 0.0;
        int total_collisions = 0;
        
        #pragma omp parallel for reduction(+:sum_v2,total_collisions)
        for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
            IonState ion;
            ion.species_id = "H3O+";
            ion.mass_kg = mass_kg;
            ion.ion_charge_C = ELEM_CHARGE_C;
            ion.CCS_m2 = 24.9 * ANGSTROM2_TO_M2;
            ion.pos = Vec3{0.0, 0.0, 0.0};
            ion.vel = Vec3{v_init, 0.0, 0.0};
            
            PhysicsRng rng(42 + ion_idx);
            
            int collision_count = 0;
            for (int i = 0; i < N_STEPS; ++i) {
                bool collided = run_collision(handler, ion, dt, rng, env);
                if (collided) collision_count++;
            }
            
            double v2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
            sum_v2 += v2;
            total_collisions += collision_count;
        }
        
        double mean_v2 = sum_v2 / N_IONS;
        double KE_final_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
        double KE_expected_eV = thermal_energy_eV(T_K);
        double ratio = KE_final_eV / KE_expected_eV;
        double collision_rate = (double)total_collisions / (N_IONS * N_STEPS);
        
        ratios.push_back(ratio);
        
        std::cout << "T = " << std::setw(5) << T_K << " K:  ";
        std::cout << "KE_final = " << std::setw(8) << KE_final_eV << " eV,  ";
        std::cout << "KE_expected = " << std::setw(8) << KE_expected_eV << " eV,  ";
        std::cout << "Ratio = " << std::setw(6) << ratio << "  (" << std::setw(5) << (ratio*100) << "%),  ";
        std::cout << "Collisions = " << std::setw(5) << (collision_rate*100) << "%\n";
    }
    
    std::cout << "\n=== ANALYSIS ===\n";
    
    // Calculate mean and std dev of ratios
    double mean_ratio = 0.0;
    for (double r : ratios) mean_ratio += r;
    mean_ratio /= ratios.size();
    
    double std_dev = 0.0;
    for (double r : ratios) {
        std_dev += (r - mean_ratio) * (r - mean_ratio);
    }
    std_dev = std::sqrt(std_dev / ratios.size());
    
    std::cout << "Mean ratio: " << mean_ratio << " ± " << std_dev << "\n";
    std::cout << "Coefficient of variation: " << (std_dev / mean_ratio * 100) << "%\n\n";
    
    if (std_dev / mean_ratio < 0.05) {
        std::cout << "✓ Ratio is CONSTANT across temperatures (< 5% variation)\n";
        std::cout << "  → Bug is a constant scaling factor (~" << mean_ratio << "x)\n";
        std::cout << "  → Neutral sampling temperature dependence is CORRECT\n";
        std::cout << "  → Problem is in energy transfer magnitude, not temperature scaling\n";
    } else {
        std::cout << "✗ Ratio VARIES with temperature (> 5% variation)\n";
        std::cout << "  → Neutral sampling may be broken\n";
        std::cout << "  → Temperature dependence is incorrect\n";
    }
    
    // Test that ratios are consistent (within 10% relative variation)
    REQUIRE(std_dev / mean_ratio < 0.10);
    
    // Test that mean ratio is approximately 1.0 (correct thermalization after bug fix)
    INFO("Mean ratio across temperatures: " << mean_ratio);
    REQUIRE(mean_ratio == Approx(1.0).margin(0.15));
}

TEST_CASE("Temperature scaling: Low energy start", "[collision][thermalization][temperature]") {
    
    std::cout << "\n=== TEMPERATURE SCALING TEST (LOW ENERGY START) ===\n";
    std::cout << "Testing from low energy (0.1x thermal) to equilibrium\n\n";
    
    std::vector<double> temperatures = {150.0, 300.0, 450.0, 600.0};
    std::vector<double> ratios;
    
    const int N_IONS = 500;
    const int N_STEPS = 5000;
    const double dt = 1e-8;
    const double mass_kg = 19.0 * AMU_TO_KG;
    
    GeometryMap geometry = load_h3o_geometry();
    
    for (double T_K : temperatures) {
        EnvironmentConfig env;
        env.temperature_K = T_K;
        env.pressure_Pa = 5000.0;
        env.gas_species = "N2";
        env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
        env.compute_derived_properties();
        
        EHSSCollisionHandler handler(geometry, false);
        
        // Start ions with LOW energy (0.1x thermal)
        const double v_init = std::sqrt(0.1 * 3.0 * BOLTZMANN_CONSTANT * T_K / mass_kg);
        
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
            
            PhysicsRng rng(1000 + ion_idx);
            
            for (int i = 0; i < N_STEPS; ++i) {
                run_collision(handler, ion, dt, rng, env);
            }
            
            double v2 = ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
            sum_v2 += v2;
        }
        
        double mean_v2 = sum_v2 / N_IONS;
        double KE_final_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
        double KE_expected_eV = thermal_energy_eV(T_K);
        double ratio = KE_final_eV / KE_expected_eV;
        
        ratios.push_back(ratio);
        
        std::cout << "T = " << std::setw(5) << T_K << " K:  ";
        std::cout << "Ratio = " << std::setw(6) << ratio << "  (" << std::setw(5) << (ratio*100) << "%)\n";
    }
    
    std::cout << "\n";
    
    // Check consistency
    double mean_ratio = 0.0;
    for (double r : ratios) mean_ratio += r;
    mean_ratio /= ratios.size();
    
    double std_dev = 0.0;
    for (double r : ratios) std_dev += (r - mean_ratio) * (r - mean_ratio);
    std_dev = std::sqrt(std_dev / ratios.size());
    
    std::cout << "LOW energy start: Mean ratio = " << mean_ratio << " ± " << std_dev << "\n";
    
    REQUIRE(std_dev / mean_ratio < 0.10);
}

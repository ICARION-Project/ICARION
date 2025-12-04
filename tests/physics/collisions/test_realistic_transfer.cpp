// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

// Test: Energy transfer with SAMPLED neutrals (realistic)

#include "core/physics/collisions/core/CollisionKernels.h"
#include "core/physics/collisions/core/VelocitySampling.h"
#include "core/types/CollisionTypes.h"
#include "utils/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ICARION;
using namespace ICARION::physics::collision_core;
using ICARION::physics::PhysicsRng;
using ICARION::physics::EHSSParams;

int main() {
    EHSSParams p;
    p.mi = 19.0 * AMU_TO_KG;
    p.mn = 28.0 * AMU_TO_KG;  // N2
    p.kB = BOLTZMANN_CONSTANT;
    p.Tn = 300.0;
    p.ubx = 0.0;
    p.uby = 0.0;
    p.ubz = 0.0;
    p.sigma_eff = 25e-20;
    p.Rn = std::sqrt(p.sigma_eff / M_PI);
    
    PhysicsRng rng(42);
    
    // Test: Ion with high energy, neutral SAMPLED from Maxwell-Boltzmann
    const int N_trials = 10000;
    const double v_ion_init = 3000.0;  // Fast ion
    
    double sum_E_ion_before = 0.0;
    double sum_E_ion_after = 0.0;
    double sum_E_neutral_before = 0.0;
    double sum_E_neutral_after = 0.0;
    
    for (int i = 0; i < N_trials; ++i) {
        Vec3 v_ion{v_ion_init, 0.0, 0.0};
        Vec3 v_neutral = VelocitySampling::sample_neutral_velocity(p.Tn, p.mn, Vec3{p.ubx, p.uby, p.ubz}, rng);
        
        // Lab frame energies
        double E_ion_before = 0.5 * p.mi * (v_ion.x*v_ion.x + v_ion.y*v_ion.y + v_ion.z*v_ion.z);
        double E_neutral_before = 0.5 * p.mn * (v_neutral.x*v_neutral.x + v_neutral.y*v_neutral.y + v_neutral.z*v_neutral.z);
        
        sum_E_ion_before += E_ion_before;
        sum_E_neutral_before += E_neutral_before;
        
        // Collision
        Vec3 v_ion_after = CollisionKernels::hss_collision(v_ion, v_neutral, p.mi, p.mn, rng);
        
        // Neutral velocity after (from momentum conservation)
        Vec3 momentum = v_ion * p.mi + v_neutral * p.mn;
        Vec3 v_neutral_after = (momentum - v_ion_after * p.mi) * (1.0 / p.mn);
        
        double E_ion_after = 0.5 * p.mi * (v_ion_after.x*v_ion_after.x + v_ion_after.y*v_ion_after.y + v_ion_after.z*v_ion_after.z);
        double E_neutral_after = 0.5 * p.mn * (v_neutral_after.x*v_neutral_after.x + v_neutral_after.y*v_neutral_after.y + v_neutral_after.z*v_neutral_after.z);
        
        sum_E_ion_after += E_ion_after;
        sum_E_neutral_after += E_neutral_after;
    }
    
    double avg_E_ion_before = sum_E_ion_before / N_trials;
    double avg_E_ion_after = sum_E_ion_after / N_trials;
    double avg_E_neutral_before = sum_E_neutral_before / N_trials;
    double avg_E_neutral_after = sum_E_neutral_after / N_trials;
    
    double E_thermal_neutral = 1.5 * p.kB * p.Tn;
    
    std::cout << std::fixed << std::setprecision(6);
    std::cout << "=== ENERGY TRANSFER TEST (SAMPLED NEUTRALS) ===\n\n";
    std::cout << "Ion: H3O+ (19 amu), Neutral: N2 (28 amu), T = 300 K\n";
    std::cout << "Ion starts at " << v_ion_init << " m/s\n";
    std::cout << "Neutrals sampled from Maxwell-Boltzmann at 300K\n";
    std::cout << N_trials << " collisions\n\n";
    
    std::cout << "BEFORE collisions (average):\n";
    std::cout << "  E_ion:     " << avg_E_ion_before / ELEM_CHARGE_C << " eV\n";
    std::cout << "  E_neutral: " << avg_E_neutral_before / ELEM_CHARGE_C << " eV  (expected: " << E_thermal_neutral / ELEM_CHARGE_C << " eV)\n\n";
    
    std::cout << "AFTER collisions (average):\n";
    std::cout << "  E_ion:     " << avg_E_ion_after / ELEM_CHARGE_C << " eV\n";
    std::cout << "  E_neutral: " << avg_E_neutral_after / ELEM_CHARGE_C << " eV\n\n";
    
    std::cout << "Energy change:\n";
    std::cout << "  ΔE_ion:     " << (avg_E_ion_after - avg_E_ion_before) / ELEM_CHARGE_C << " eV  (";
    std::cout << ((avg_E_ion_after - avg_E_ion_before) / avg_E_ion_before * 100) << "%)\n";
    std::cout << "  ΔE_neutral: " << (avg_E_neutral_after - avg_E_neutral_before) / ELEM_CHARGE_C << " eV  (";
    std::cout << ((avg_E_neutral_after - avg_E_neutral_before) / avg_E_neutral_before * 100) << "%)\n\n";
    
    // Check energy conservation
    double E_total_before = avg_E_ion_before + avg_E_neutral_before;
    double E_total_after = avg_E_ion_after + avg_E_neutral_after;
    std::cout << "Total energy:\n";
    std::cout << "  Before: " << E_total_before / ELEM_CHARGE_C << " eV\n";
    std::cout << "  After:  " << E_total_after / ELEM_CHARGE_C << " eV\n";
    std::cout << "  ΔE_total: " << (E_total_after - E_total_before) / ELEM_CHARGE_C << " eV  (should be ~0)\n\n";
    
    // What should happen in equilibrium?
    double E_thermal = 1.5 * p.kB * p.Tn;
    std::cout << "Expected equilibrium:\n";
    std::cout << "  Both ion and neutral should have E = (3/2)kT = " << E_thermal / ELEM_CHARGE_C << " eV\n";
    std::cout << "  Current ion energy: " << avg_E_ion_after / ELEM_CHARGE_C << " eV  (";
    std::cout << (avg_E_ion_after / E_thermal * 100) << "% of thermal)\n\n";
    
    // Test 2: Multiple collisions with many ions (ensemble average)
    std::cout << "=== ENSEMBLE AVERAGE: MANY IONS, MANY COLLISIONS ===\n\n";
    
    const int N_IONS_TEST = 1000;
    std::vector<int> collision_counts = {10, 50, 100, 500, 1000, 5000};
    
    std::cout << "Starting 1000 ions at " << v_ion_init << " m/s\n";
    std::cout << "Each ion collides N times with freshly sampled neutrals\n\n";
    
    std::cout << "N_collisions    <E_ion> (eV)    Ratio\n";
    std::cout << "------------    ------------    -----\n";
    
    for (int N_coll : collision_counts) {
        double sum_E = 0.0;
        
        #pragma omp parallel for reduction(+:sum_E)
        for (int ion_idx = 0; ion_idx < N_IONS_TEST; ++ion_idx) {
            Vec3 v_ion_local{v_ion_init, 0.0, 0.0};
            PhysicsRng rng_local(1000 + ion_idx);
            
            for (int i = 0; i < N_coll; ++i) {
                Vec3 v_neutral = VelocitySampling::sample_neutral_velocity(
                    p.Tn, p.mn, Vec3{p.ubx, p.uby, p.ubz}, rng_local
                );
                v_ion_local = CollisionKernels::hss_collision(v_ion_local, v_neutral, p.mi, p.mn, rng_local);
            }
            
            double E = 0.5 * p.mi * (v_ion_local.x*v_ion_local.x + v_ion_local.y*v_ion_local.y + v_ion_local.z*v_ion_local.z);
            sum_E += E;
        }
        
        double avg_E = sum_E / N_IONS_TEST;
        double ratio = avg_E / E_thermal;
        
        std::cout << std::setw(12) << N_coll << "    ";
        std::cout << std::setw(12) << avg_E / ELEM_CHARGE_C << "    ";
        std::cout << std::setw(5) << (ratio * 100) << "%\n";
    }
    
    std::cout << "\n✓ If ratio converges to 100% → collision physics is correct\n";
    std::cout << "✓ If ratio converges to <100% → systematic energy loss bug\n";
    std::cout << "✓ If ratio converges to >100% → systematic energy gain bug\n";
    
    return 0;
}

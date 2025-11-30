// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

// Diagnostic test to understand thermalization energy distribution

#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/types/IonState.h"
#include "core/config/types/EnvironmentConfig.h"
#include "utils/constants.h"
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace ICARION::physics;
using namespace ICARION::config;
using namespace ICARION::core;

int main() {
    const int N_IONS = 5000;
    const double T_K = 300.0;
    const double mass_kg = 19.0 * AMU_TO_KG;
    
    EnvironmentConfig env;
    env.temperature_K = T_K;
    env.pressure_Pa = 2000.0;
    env.gas_species = "Ar";  // Changed to Helium
    env.gas_velocity_m_s = Vec3{0.0, 0.0, 0.0};
    env.compute_derived_properties();
    
    HSSCollisionHandler handler(false);
    
    const int N_STEPS = 100000;
    const double dt = 1e-10;  // 0.1 ns (safe for P ≈ 10%)
    
    double sum_vx2 = 0.0, sum_vy2 = 0.0, sum_vz2 = 0.0;
    double sum_v2 = 0.0;
    
    #pragma omp parallel for reduction(+:sum_vx2,sum_vy2,sum_vz2,sum_v2)
    for (int ion_idx = 0; ion_idx < N_IONS; ++ion_idx) {
        IonState ion;
        ion.species_id = "H3O+";
        ion.mass_kg = mass_kg;
        ion.ion_charge_C = ELEM_CHARGE_C;
        ion.CCS_m2 = 24.9 * ANGSTROM2_TO_M2;
        ion.pos = Vec3{0.0, 0.0, 0.0};
        ion.vel = Vec3{6000.0, 0.0, 0.0};  // Start with high energy in x
        
        EhssRng rng(456 + ion_idx);
        
        for (int i = 0; i < N_STEPS; ++i) {
            handler.handle_collision(ion, dt, rng, env);
        }
        
        sum_vx2 += ion.vel.x * ion.vel.x;
        sum_vy2 += ion.vel.y * ion.vel.y;
        sum_vz2 += ion.vel.z * ion.vel.z;
        sum_v2 += ion.vel.x * ion.vel.x + ion.vel.y * ion.vel.y + ion.vel.z * ion.vel.z;
    }
    
    double mean_vx2 = sum_vx2 / N_IONS;
    double mean_vy2 = sum_vy2 / N_IONS;
    double mean_vz2 = sum_vz2 / N_IONS;
    double mean_v2 = sum_v2 / N_IONS;
    
    double KE_x_eV = 0.5 * mass_kg * mean_vx2 / ELEM_CHARGE_C;
    double KE_y_eV = 0.5 * mass_kg * mean_vy2 / ELEM_CHARGE_C;
    double KE_z_eV = 0.5 * mass_kg * mean_vz2 / ELEM_CHARGE_C;
    double KE_tot_eV = 0.5 * mass_kg * mean_v2 / ELEM_CHARGE_C;
    
    double KE_per_DOF_expected = 0.5 * BOLTZMANN_CONSTANT * T_K / ELEM_CHARGE_C;
    double KE_3DOF_expected = 1.5 * BOLTZMANN_CONSTANT * T_K / ELEM_CHARGE_C;
    
    std::cout << std::setprecision(6) << std::fixed;
    std::cout << "\n=== Thermalization Energy Diagnostic ===\n\n";
    std::cout << "T = " << T_K << " K, N = " << N_IONS << " ions, " << N_STEPS << " steps\n\n";
    
    std::cout << "Expected per DOF: " << KE_per_DOF_expected << " eV\n";
    std::cout << "Expected total:   " << KE_3DOF_expected << " eV\n\n";
    
    std::cout << "Measured:\n";
    std::cout << "  <KE_x> = " << KE_x_eV << " eV  (ratio = " << KE_x_eV/KE_per_DOF_expected << ")\n";
    std::cout << "  <KE_y> = " << KE_y_eV << " eV  (ratio = " << KE_y_eV/KE_per_DOF_expected << ")\n";
    std::cout << "  <KE_z> = " << KE_z_eV << " eV  (ratio = " << KE_z_eV/KE_per_DOF_expected << ")\n";
    std::cout << "  <KE_tot> = " << KE_tot_eV << " eV  (ratio = " << KE_tot_eV/KE_3DOF_expected << ")\n\n";
    
    std::cout << "Sum of components: " << (KE_x_eV + KE_y_eV + KE_z_eV) << " eV\n";
    std::cout << "Direct total:      " << KE_tot_eV << " eV\n\n";
    
    if (std::abs(KE_tot_eV / KE_3DOF_expected - 1.0) < 0.15) {
        std::cout << "✅ Energy matches expected 3D thermal distribution\n";
    } else {
        std::cout << "❌ Energy does NOT match expected thermal distribution!\n";
        std::cout << "   Deviation: " << (KE_tot_eV / KE_3DOF_expected - 1.0) * 100 << "%\n";
    }
    
    return 0;
}

// SPDX-License-Identifier: MIT
// Debug script: Track what happens in first HSS timestep

#include <iostream>
#include <iomanip>
#include "core/physics/collisions/HSSCollisionHandler.h"
#include "core/types/CollisionTypes.h"
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/IonState.h"
#include "utils/constants.h"

using namespace ICARION;
using ICARION::physics::EhssRng;

int main() {
    std::cout << std::fixed << std::setprecision(8);
    
    // Create H3O+ ion at entrance with small thermal velocity
    core::IonState ion;
    ion.species_id = "H3O+";
    ion.pos = {0.0, 0.0, 0.0};  // At entrance
    ion.mass_kg = 19.0 * AMU_TO_KG;
    ion.ion_charge_C = ELEM_CHARGE_C;
    ion.CCS_m2 = 24.9e-20;
    
    // Small thermal velocity
    double T_K = 300.0;
    double v_th = std::sqrt(BOLTZMANN_CONSTANT * T_K / ion.mass_kg);
    ion.vel = {0.0, 0.0, 0.1 * v_th};
    
    std::cout << "=== Initial State ===\n";
    std::cout << "Position: (" << ion.pos.x*1000 << ", " << ion.pos.y*1000 << ", " << ion.pos.z*1000 << ") mm\n";
    std::cout << "Velocity: (" << ion.vel.x << ", " << ion.vel.y << ", " << ion.vel.z << ") m/s\n";
    std::cout << "v_thermal: " << v_th << " m/s\n";
    std::cout << "\n";
    
    // Environment: 1 atm N2, 300 K
    config::EnvironmentConfig env;
    env.temperature_K = 300.0;
    env.pressure_Pa = 101325.0;
    env.gas_species = "N2";
    env.gas_velocity_m_s = {0.0, 0.0, 0.0};
    env.compute_derived_properties();
    
    std::cout << "=== Environment ===\n";
    std::cout << "Gas: " << env.gas_species << "\n";
    std::cout << "Pressure: " << env.pressure_Pa << " Pa\n";
    std::cout << "Temperature: " << env.temperature_K << " K\n";
    std::cout << "Number density: " << env.number_density_m3 << " m^-3\n";
    std::cout << "Mean free path: " << (1.0 / (env.number_density_m3 * ion.CCS_m2)) * 1e6 << " μm\n";
    std::cout << "\n";
    
    // HSS handler
    physics::HSSCollisionHandler handler(false);
    EhssRng rng(42);
    
    // Timestep
    double dt = 1e-9;  // 1 ns
    
    std::cout << "=== First Timestep (dt = " << dt*1e9 << " ns) ===\n";
    
    // Store initial state
    auto ion_before = ion;
    
    // Handle collision
    bool collision_occurred = handler.handle_collision(ion, dt, rng, env);
    
    std::cout << "Collision occurred: " << (collision_occurred ? "YES" : "NO") << "\n";
    std::cout << "\n";
    
    std::cout << "=== After Collision ===\n";
    std::cout << "Position: (" << ion.pos.x*1000 << ", " << ion.pos.y*1000 << ", " << ion.pos.z*1000 << ") mm\n";
    std::cout << "Velocity: (" << ion.vel.x << ", " << ion.vel.y << ", " << ion.vel.z << ") m/s\n";
    std::cout << "Delta position: (" << (ion.pos.x-ion_before.pos.x)*1e6 << ", " 
              << (ion.pos.y-ion_before.pos.y)*1e6 << ", "
              << (ion.pos.z-ion_before.pos.z)*1e6 << ") μm\n";
    std::cout << "Delta velocity: (" << (ion.vel.x-ion_before.vel.x) << ", "
              << (ion.vel.y-ion_before.vel.y) << ", "
              << (ion.vel.z-ion_before.vel.z) << ") m/s\n";
    std::cout << "\n";
    
    // Check if velocity is reasonable
    double speed = std::sqrt(ion.vel.x*ion.vel.x + ion.vel.y*ion.vel.y + ion.vel.z*ion.vel.z);
    std::cout << "Final speed: " << speed << " m/s\n";
    std::cout << "Ratio to thermal: " << speed / v_th << "\n";
    
    if (speed > 10.0 * v_th) {
        std::cout << "\n⚠️  WARNING: Velocity >> thermal velocity!\n";
        std::cout << "This suggests HSS collision imparts too much momentum.\n";
    }
    
    return 0;
}

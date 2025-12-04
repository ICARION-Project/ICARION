// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#include "OUCollisionHandler.h"
#include "core/physics/collisions/core/CollisionKernels.h"
#include <stdexcept>

namespace ICARION::physics {

OUCollisionHandler::OUCollisionHandler(double gamma_coefficient, bool apply_damping)
    : gamma_(gamma_coefficient)
    , apply_damping_(apply_damping)
{
    if (gamma_ <= 0.0) {
        throw std::invalid_argument("OUCollisionHandler: gamma_coefficient must be > 0!");
    }
}

bool OUCollisionHandler::handle_collision(
    IonState& ion,
    double dt,
    PhysicsRng& rng,
    const config::EnvironmentConfig& env
) {
    // ===================================================================
    // READ TEMPERATURE DIRECTLY FROM ENV (SSOT!)
    // ===================================================================
    const double T_K = env.temperature_K;
    
    // Apply Ornstein-Uhlenbeck velocity kick
    // Uses modern CollisionKernels module
    // apply_damping_ controls whether damping is applied (false when using DampingForce)
    collision_core::CollisionKernels::ou_velocity_update(
        ion, rng, dt, gamma_, T_K, env.gas_velocity_m_s, apply_damping_
    );
    
    return true;  // Always "collides" (continuous process)
}

bool OUCollisionHandler::handle_collision_soa(
    core::IonCollisionData& view,
    double dt,
    PhysicsRng& rng,
    const config::EnvironmentConfig& env
) {
    IonState ion;
    ion.vel = view.kin.vel();
    ion.mass_kg = view.kin.get_mass();
    ion.ion_charge_C = view.kin.get_charge();
    ion.CCS_m2 = view.get_CCS();
    
    bool occurred = handle_collision(ion, dt, rng, env);
    view.kin.set_vel(ion.vel);
    return occurred;
}

} // namespace ICARION::physics

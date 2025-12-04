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

} // namespace ICARION::physics

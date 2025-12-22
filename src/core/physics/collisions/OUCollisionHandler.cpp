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

    const double T_K = env.temperature_K;

    collision_core::CollisionKernels::ou_velocity_update(
        ion, rng, dt, gamma_, T_K, env.gas_velocity_m_s, apply_damping_
    );

    view.kin.set_vel(ion.vel);
    return true;  // Always "collides" (continuous process)
}

} // namespace ICARION::physics

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "OUCollisionHandler.h"
#include "collisionHelpers.h"
#include <stdexcept>

namespace ICARION::physics {

OUCollisionHandler::OUCollisionHandler(double gamma_coefficient)
    : gamma_(gamma_coefficient)
{
    if (gamma_ <= 0.0) {
        throw std::invalid_argument("OUCollisionHandler: gamma_coefficient must be > 0!");
    }
}

bool OUCollisionHandler::handle_collision(
    IonState& ion,
    double dt,
    EhssRng& rng,
    const config::EnvironmentConfig& env
) {
    // ===================================================================
    // READ TEMPERATURE DIRECTLY FROM ENV (SSOT!)
    // ===================================================================
    const double T_K = env.temperature_K;
    
    // Apply Ornstein-Uhlenbeck velocity kick
    // Uses existing helper function from collisionHelpers.h
    apply_ou_velocity_kick(ion, rng, dt, gamma_, T_K);
    
    return true;  // Always "collides" (continuous process)
}

} // namespace ICARION::physics

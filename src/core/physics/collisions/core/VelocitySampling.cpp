// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "VelocitySampling.h"
#include <cmath>

namespace ICARION::physics::collision_core {

Vec3 VelocitySampling::sample_neutral_velocity(
    double temperature_K,
    double mass_kg,
    const Vec3& flow_velocity,
    EhssRng& rng
) {
    // PHYSICS: Exact implementation from collisionHelpers.cpp::sample_neutral_velocity()
    // NO CHANGES to algorithm - only refactored for clarity
    
    // Compute thermal velocity width: sigma = sqrt(kB*T/m)
    double sigma = thermal_velocity_width(temperature_K, mass_kg);
    
    // Sample three independent thermal velocity components
    // Using Box-Muller transform for Gaussian distribution
    Vec3 v_thermal{
        box_muller_sample(sigma, rng),
        box_muller_sample(sigma, rng),
        box_muller_sample(sigma, rng)
    };
    
    // Add bulk flow velocity to get lab-frame velocity
    return Vec3{
        v_thermal.x + flow_velocity.x,
        v_thermal.y + flow_velocity.y,
        v_thermal.z + flow_velocity.z
    };
}

double VelocitySampling::sample_thermal_component(
    double temperature_K,
    double mass_kg,
    EhssRng& rng
) {
    double sigma = thermal_velocity_width(temperature_K, mass_kg);
    return box_muller_sample(sigma, rng);
}

double VelocitySampling::thermal_velocity_width(
    double temperature_K,
    double mass_kg
) {
    // Thermal velocity width: sqrt(kB*T/m)
    // This is the standard deviation of the Maxwell-Boltzmann distribution
    // for one velocity component
    return std::sqrt(BOLTZMANN_CONSTANT * temperature_K / mass_kg);
}

double VelocitySampling::box_muller_sample(double sigma, EhssRng& rng) {
    // PHYSICS: Box-Muller transform (exact copy from old implementation)
    // Generates N(0, sigma) from two uniform random numbers
    
    double u1 = rng.uniform01();
    double u2 = rng.uniform01();
    
    // Box-Muller: sqrt(-2*ln(u1)) * cos(2*pi*u2) gives N(0,1)
    // Multiply by sigma to get N(0, sigma)
    return sigma * std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * M_PI * u2);
}

} // namespace ICARION::physics::collision_core

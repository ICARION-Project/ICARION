// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file VelocitySampling.h
 * @brief Velocity distribution sampling for collision models
 * 
 * Generates velocities from thermal distributions (Maxwell-Boltzmann).
 * Pure stateless functions for sampling gas molecule velocities.
 * 
 * **Physics:**
 * - Maxwell-Boltzmann distribution: f(v) ∝ exp(-mv²/2kT)
 * - Thermal velocity width: σ = sqrt(kB*T/m)
 * - Box-Muller transform for Gaussian sampling
 * 
 * @note This module contains PHYSICS CORE code. Changes must be validated
 *       against thermalization tests and statistical distributions.
 */

#pragma once

#include "core/types/Vec3.h"
#include "core/types/CollisionTypes.h"  // For EhssRng
#include "utils/constants.h"

namespace ICARION::physics::collision_core {

/**
 * @brief Velocity distribution sampling for collision models
 * 
 * Generates velocities from thermal distributions (Maxwell-Boltzmann).
 * Pure static functions (no state, no side effects).
 */
class VelocitySampling {
public:
    /**
     * @brief Sample neutral molecule velocity from Maxwell-Boltzmann distribution
     * 
     * Generates thermal velocity using Box-Muller transform:
     * v_thermal ~ N(0, sqrt(kB*T/m)) for each component (x,y,z)
     * 
     * Then adds bulk flow velocity: v_lab = v_thermal + v_bulk
     * 
     * @param temperature_K Gas temperature [K]
     * @param mass_kg Neutral molecule mass [kg]
     * @param flow_velocity Bulk flow velocity [m/s] (lab frame)
     * @param rng Random number generator
     * 
     * @return Velocity [m/s] in lab frame
     * 
     * @pre temperature_K > 0, mass_kg > 0
     * @post result has thermal spread sqrt(kB*T/m) around flow_velocity
     * 
     * @note Equivalent to old `sample_neutral_velocity()` from collisionHelpers.h
     */
    static Vec3 sample_neutral_velocity(
        double temperature_K,
        double mass_kg,
        const Vec3& flow_velocity,
        EhssRng& rng
    );
    
    /**
     * @brief Sample thermal velocity component (1D Gaussian)
     * 
     * Uses Box-Muller transform to generate one component of thermal velocity
     * from N(0, sqrt(kB*T/m)) distribution.
     * 
     * @param temperature_K Temperature [K]
     * @param mass_kg Particle mass [kg]
     * @param rng Random number generator
     * 
     * @return Velocity component [m/s] from N(0, sqrt(kB*T/m))
     * 
     * @pre temperature_K > 0, mass_kg > 0
     */
    static double sample_thermal_component(
        double temperature_K,
        double mass_kg,
        EhssRng& rng
    );
    
    /**
     * @brief Compute thermal velocity width (1σ standard deviation)
     * 
     * Returns the characteristic thermal velocity: sqrt(kB*T/m)
     * This is the standard deviation of the Maxwell-Boltzmann distribution
     * for one velocity component.
     * 
     * @param temperature_K Temperature [K]
     * @param mass_kg Particle mass [kg]
     * 
     * @return sqrt(kB*T/m) [m/s]
     * 
     * @pre temperature_K > 0, mass_kg > 0
     */
    static double thermal_velocity_width(
        double temperature_K,
        double mass_kg
    );

private:
    /**
     * @brief Box-Muller transform for Gaussian sampling
     * 
     * Generates a random sample from N(0, sigma) using two uniform random numbers.
     * 
     * @param sigma Standard deviation
     * @param rng Random number generator
     * @return Random sample from N(0, sigma)
     */
    static double box_muller_sample(double sigma, EhssRng& rng);
};

} // namespace ICARION::physics::collision_core

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file HSSCollisionHandler.h
 * @brief Hard-Sphere Stochastic (HSS) collision handler
 * 
 * Implements isotropic hard-sphere scattering with random deflection angle.
 * Uses effective collision cross-section (single sphere approximation).
 * Faster than EHSS but less physically accurate for non-spherical molecules.
 * 
 * **Physics:**
 * - Isotropic scattering (random deflection angle in COM frame)
 * - Uses stored effective CCS (ion.CCS_m2)
 * - Conserves momentum and energy
 * 
 * **SSOT Design:**
 * - Reads gas properties directly from `EnvironmentConfig`
 * - No parameter copies or intermediate structs
 * 
 * @date 2025-11-21
 * @version 1.0
 */

#pragma once

#include "ICollisionHandler.h"
#include <unordered_map>

namespace ICARION::physics {

/**
 * @brief HSS (Hard-Sphere Stochastic) collision handler
 * 
 * Implements simple isotropic hard-sphere scattering model.
 * Uses effective collision cross-section without molecular geometry.
 * 
 * **Use cases:**
 * - Fast simulations with many ions
 * - Spherical or near-spherical molecules
 * - When geometry data is unavailable
 * 
 * **Performance:**
 * - Much faster than EHSS (no geometry sampling)
 * - O(1) collision detection
 * - Recommended for large ensembles (> 10k ions)
 * 
 * **SSOT Pattern:**
 * ```cpp
 * // Correct: No parameters needed
 * HSSCollisionHandler handler;
 * 
 * // Handler reads from EnvironmentConfig directly
 * config::EnvironmentConfig env;
 * env.temperature_K = 300.0;
 * handler.handle_collision(ion, dt, rng, env);  // SSOT!
 * ```
 * 
 * @see EHSSCollisionHandler for structure-resolved model
 * @see DampingForce for deterministic collision models
 */
class HSSCollisionHandler : public ICollisionHandler {
public:
    /**
     * @brief Construct HSS handler
     * @param enable_logging Enable debug logging (CSV output)
     */
    explicit HSSCollisionHandler(bool enable_logging = false);
    
    /**
     * @brief Handle HSS collision for single timestep
     * 
     * **Algorithm:**
     * 1. Compute collision probability from mean free path (using ion.CCS_m2)
     * 2. If collision occurs:
     *    - Sample neutral velocity from Maxwell-Boltzmann distribution
     *    - Transform to center-of-mass frame
     *    - Apply random isotropic deflection
     *    - Transform back to lab frame
     * 
     * **SSOT:** Reads gas properties directly from `env`:
     * - env.temperature_K → thermal velocity
     * - env.particle_density_m_3 → collision frequency
     * - env.neutral_mass_kg → reduced mass
     * - env.gas_velocity_m_s → bulk flow
     * 
     * @param[in,out] ion Ion state (velocity modified if collision occurs)
     * @param[in] dt Timestep [s]
     * @param[in,out] rng Random number generator
     * @param[in] env Environment configuration (SSOT!)
     * 
     * @return true if collision occurred, false otherwise
     */
    bool handle_collision(
        IonState& ion,
        double dt,
        EhssRng& rng,
        const config::EnvironmentConfig& env
    ) override;
    
    std::string name() const override { return "HSS"; }
    
    CollisionStats get_stats() const override { return stats_; }
    void reset_stats() override { stats_ = {}; collisions_by_species_.clear(); }
    const std::unordered_map<std::string, size_t>& collisions_by_species() const { return collisions_by_species_; }
    
private:
    bool enable_logging_;
    mutable CollisionStats stats_;

    bool handle_single_gas(IonState& ion, double dt, EhssRng& rng, const config::EnvironmentConfig& env);
    std::unordered_map<std::string, size_t> collisions_by_species_;
};

} // namespace ICARION::physics

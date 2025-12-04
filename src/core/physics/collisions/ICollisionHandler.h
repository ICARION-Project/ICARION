// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file ICollisionHandler.h
 * @brief Interface for stochastic collision models
 * 
 * Discrete collision handlers (EHSS, HSS, OU) update ion velocity stochastically
 * using environment parameters supplied at call time. Deterministic damping lives in
 * DampingForce.
 */

#pragma once

#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"  // For IonCollisionData view
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/CollisionTypes.h"  // PhysicsRng
#include <string>
#include <cstddef>

namespace ICARION::physics {

/**
 * @brief Collision statistics (for tracking/debugging)
 */
struct CollisionStats {
    size_t total_collisions = 0;        ///< Total number of collision events
    size_t rejected_collisions = 0;     ///< Collisions rejected (e.g., low probability)
    double average_collision_rate = 0.0;///< Mean collision rate [Hz]
};

/**
 * @brief Abstract interface for stochastic collision handlers
 * 
 * All collision models (EHSS, HSS, OU) implement this interface.
 * 
 * **IMPORTANT:** This is ONLY for **stochastic** collision models!
 * Deterministic models (Friction, Langevin, HardSphere) use `DampingForce` instead.
 * 
 * **SSOT Design:**
 * - Handler reads parameters directly from `EnvironmentConfig` (no copies!)
 * - No intermediate structs (e.g., CollisionContext)
 * - Config references must outlive handler instance
 * 
 * **Collision Models:**
 * - **EHSS** (Explicit Hard-Sphere Scattering): Structure-resolved, atom-centered spheres
 * - **HSS** (Hard-Sphere Stochastic): Isotropic scattering, single effective sphere
 * - **OU** (Ornstein-Uhlenbeck): Thermal velocity kicks (add-on for deterministic models)
 * 
 * @see DampingForce for deterministic collision models
 * @see CollisionHandlerFactory for creating handlers from config
 */
class ICollisionHandler {
public:
    virtual ~ICollisionHandler() = default;
    
    /**
     * @brief Handle collision for single timestep
     * 
     * Determines if collision occurs (probabilistic) and updates ion velocity accordingly.
     * 
     * @param[in,out] ion Ion state (velocity modified in-place if collision occurs)
     * @param[in] dt Timestep [s]
     * @param[in,out] rng Random number generator
     * @param[in] env Environment configuration (SSOT - direct reference!)
     * 
     * @return true if collision occurred, false otherwise
     * 
     * @note Thread-safety: Not thread-safe! Each thread needs separate handler + RNG instance.
     * 
     * @example
     * ```cpp
     * config::EnvironmentConfig env;
     * env.temperature_K = 300.0;
     * env.pressure_Pa = 101325.0;
     * 
     * PhysicsRng rng(12345);
     * IonState ion;
     * 
     * bool collision_occurred = handler->handle_collision(ion, 1e-9, rng, env);
     * ```
     */
    virtual bool handle_collision(
        IonState& ion,
        double dt,
        PhysicsRng& rng,
        const config::EnvironmentConfig& env
    ) = 0;
    
    /**
     * @brief Handle collision using SoA view (Phase 3 - cache-optimized)
     * 
     * Zero-copy access to ion data via view struct.
     * Default implementation converts to IonState and calls handle_collision().
     * 
     * @param[in,out] view Ion collision data view (velocity modified in-place)
     * @param[in] dt Timestep [s]
     * @param[in,out] rng Random number generator
     * @param[in] env Environment configuration
     * 
     * @return true if collision occurred, false otherwise
     * 
     * @note Override this for optimal SoA performance. Default wrapper provided for compatibility.
     */
    virtual bool handle_collision_soa(
        core::IonCollisionData& view,
        double dt,
        PhysicsRng& rng,
        const config::EnvironmentConfig& env
    ) {
        // Default: convert to IonState and call legacy method
        IonState ion;
        ion.pos = view.kin.pos();
        ion.vel = view.kin.vel();
        ion.mass_kg = view.kin.get_mass();
        ion.ion_charge_C = view.kin.get_charge();
        ion.CCS_m2 = view.get_CCS();
        
        bool result = handle_collision(ion, dt, rng, env);
        
        // Write back modified velocity
        view.kin.set_vel(ion.vel);
        return result;
    }
    
    /**
     * @brief Get collision model name
     * 
     * @return Model name (e.g., "EHSS", "HSS", "OU")
     */
    virtual std::string name() const = 0;
    
    /**
     * @brief Get collision statistics (for debugging/logging)
     * 
     * @return Collision statistics (total collisions, rejection rate, etc.)
     */
    virtual CollisionStats get_stats() const { return {}; }
    
    /**
     * @brief Reset statistics counters
     */
    virtual void reset_stats() {}
};

} // namespace ICARION::physics

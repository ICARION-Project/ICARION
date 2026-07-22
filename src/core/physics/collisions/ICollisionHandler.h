// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file ICollisionHandler.h
 * @brief Interface for stochastic collision models
 * 
 * Discrete collision handlers (EHSS, HSS, OU) update ion velocity stochastically
 * using environment parameters supplied at call time. Deterministic damping lives in
 * DampingForce.
 */

#pragma once

#include "core/types/IonEnsemble.h"  // For IonCollisionData view
#include "core/config/types/EnvironmentConfig.h"
#include "core/types/CollisionTypes.h"  // PhysicsRng
#include <string>
#include <cstddef>
#include <vector>

namespace ICARION::physics {

/**
 * @brief Collision statistics (for tracking/debugging)
 */
struct CollisionStats {
    size_t total_collisions = 0;        ///< Total number of collision events
    size_t rejected_collisions = 0;     ///< Collisions rejected (e.g., low probability)
    double average_collision_rate = 0.0;///< Mean collision rate [Hz]
};

struct CollisionEventDiagnostics {
    double v_rel_before_m_s = 0.0;  ///< Exact pre-collision relative speed [m/s]
    double sigma_mt_m2 = 0.0;       ///< Momentum-transfer cross section used for this event [m^2]
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
     * @brief Handle collision using SoA view (cache-optimized hot path)
     * 
     * @param[in,out] view Ion collision data view (velocity modified in-place)
     * @param[in] dt Timestep [s]
     * @param[in,out] rng Random number generator
     * @param[in] env Environment configuration
     * 
     * @return true if collision occurred, false otherwise
     */
    virtual bool handle_collision(
        core::IonCollisionData& view,
        double dt,
        PhysicsRng& rng,
        const config::EnvironmentConfig& env,
        CollisionEventDiagnostics* diagnostics = nullptr
    ) = 0;
    
    /**
     * @brief Get collision model name
     * 
     * @return Model name (e.g., "EHSS", "HSS", "OU")
     */
    virtual std::string name() const = 0;
    
    /**
     * @brief Get collision statistics (for debugging/logging)
     *
     * Call only after the simulation engine's collision-worker barrier.
     * Implementations may reduce worker-owned counters without synchronizing
     * concurrent handle_collision() writes.
     * 
     * @return Collision statistics (total collisions, rejection rate, etc.)
     */
    virtual CollisionStats get_stats() const { return {}; }
    
    /**
     * @brief Reset statistics counters
     *
     * Call only while no handle_collision() invocations are active.
     */
    virtual void reset_stats() {}

    /**
     * @brief Whether this handler exposes a batch path (e.g., GPU).
     */
    virtual bool supports_batch() const { return false; }

    /**
     * @brief Optional batch API for accelerators.
     *
     * @param ensemble Ion ensemble (SoA) to mutate.
     * @param ion_indices Indices participating in this batch (typically same domain).
     * @param dt Timestep.
     * @param env Domain environment.
     * @param rng_pool Per-ion RNGs for CPU fallback (GPU implementations may ignore).
     * @return true if the handler processed the batch (GPU, etc.), false to request CPU fallback.
     */
    virtual bool handle_batch(
        core::IonEnsemble& ensemble,
        const std::vector<size_t>& ion_indices,
        double dt,
        const config::EnvironmentConfig& env,
        std::vector<physics::PhysicsRng>& rng_pool
    ) {
        (void)ensemble;
        (void)ion_indices;
        (void)dt;
        (void)env;
        (void)rng_pool;
        return false;
    }
};

} // namespace ICARION::physics

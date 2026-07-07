// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_CONFIG_PHYSICS_CONFIG_H
#define ICARION_CONFIG_PHYSICS_CONFIG_H

#include "PhysicsEnums.h"
#include "../validation/ValidationResult.h"
#include <string>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief Global physics settings
 * 
 * Physics options that apply across all domains.
 * Replaces GlobalParams (physics-related parts only).
 */
struct PhysicsConfig {
    // === Collision model ===
    CollisionModel collision_model = CollisionModel::NoCollisions;
    
    // === Feature flags ===
    bool enable_reactions = false;
    bool enable_space_charge = false;
    bool enable_space_charge_gpu = false;
    SpaceChargeModel space_charge_model_type = SpaceChargeModel::Auto; ///< Canonical. Parsed from JSON "space_charge_model".
    std::string space_charge_model = "auto"; ///< Compatibility mirror — do not dispatch on this.

    // === Stochastic collision event handling ===
    bool collision_multi_event_mode = false;   ///< If true, enable multi-collision approximation via micro-subcycling inside each dt.
    int collision_max_events_per_step = 16;    ///< Upper bound for micro-subcycles/events per dt when collision_multi_event_mode is enabled.
    bool collision_time_centered = false;      ///< If true, apply collisions in two half-steps around force integration (reduces splitting bias)
    bool collision_time_randomized = false;    ///< If true, randomize collision timing within each step (per-ion random pre/post split)
    int collision_subcycles_per_step = 1;      ///< If >1, split each collision application into equal micro-steps (recompute rates each micro-step)

    // === InteractionPotentialModel model controls ===
    // --- Orientation mode (SSOT: ipm_orientation_mode_type) ---
    IPMOrientationMode ipm_orientation_mode_type = IPMOrientationMode::Random; ///< Canonical. Parsed from JSON "ipm_orientation_mode".
    std::string ipm_orientation_mode = "random"; ///< Compatibility mirror — do not dispatch on this.
    int ipm_fixed_orientation_index = 0;           ///< Fixed orientation index used when mode=fixed (>=0)
    std::string ipm_vrel_log_prefix = "";          ///< If non-empty, enable InteractionPotentialModel v_rel histogram logging and write to this file/prefix
    std::string ipm_momentum_log_prefix = "";      ///< If non-empty, enable InteractionPotentialModel momentum diagnostic logging and write to this file/prefix
    
    // === Thermalization ===
    bool enable_ou_thermalization = false;
    bool force_ou_for_stochastic = false;
    
    /**
     * @brief Validate physics configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        // OU thermalization with stochastic models is incompatible
        if (enable_ou_thermalization && 
            (collision_model == CollisionModel::HSS || 
             collision_model == CollisionModel::EHSS)) {
            result.add_error("enable_ou_thermalization cannot be true when using stochastic collision models (HSS or EHSS). "
                           "OU is only compatible with deterministic damping models (Friction, Langevin, HardSphere).");
        }
        
        if (force_ou_for_stochastic && 
            (collision_model == CollisionModel::HSS || 
             collision_model == CollisionModel::EHSS)) {
            result.add_error("force_ou_for_stochastic cannot be true when using stochastic collision models (HSS or EHSS)");
        }

        if (collision_model == CollisionModel::InteractionPotentialModel) {
            if (ipm_orientation_mode_type == IPMOrientationMode::Fixed &&
                ipm_fixed_orientation_index < 0) {
                result.add_error("ipm_fixed_orientation_index must be >= 0");
            }
        }

        if (collision_max_events_per_step < 1) {
            result.add_error("collision_max_events_per_step must be >= 1");
        }

        if (collision_subcycles_per_step < 1) {
            result.add_error("collision_subcycles_per_step must be >= 1");
        }
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_PHYSICS_CONFIG_H

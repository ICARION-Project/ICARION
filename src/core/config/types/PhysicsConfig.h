// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_PHYSICS_CONFIG_H
#define ICARION_CONFIG_PHYSICS_CONFIG_H

#include "PhysicsEnums.h"
#include "../validation/ValidationResult.h"
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
    
    // === Thermalization ===
    bool enable_ou_thermalization = false;
    bool force_ou_for_stochastic = false;
    
    /**
     * @brief Validate physics configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        // OU thermalization with stochastic models might be redundant
        if (force_ou_for_stochastic && 
            (collision_model == CollisionModel::HSS || 
             collision_model == CollisionModel::EHSS)) {
            result.add_error("force_ou_for_stochastic cannot be true when using stochastic collision models (HSS or EHSS)");
        }
        
        return result;
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_PHYSICS_CONFIG_H

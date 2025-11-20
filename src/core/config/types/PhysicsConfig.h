// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_PHYSICS_CONFIG_H
#define ICARION_CONFIG_PHYSICS_CONFIG_H

#include "PhysicsEnums.h"
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
     * 
     * @throws std::runtime_error if invalid
     */
    void validate() const {
        
        // OU thermalization with stochastic models is redundant
        if (force_ou_for_stochastic && 
            (collision_model == CollisionModel::HSMC || 
             collision_model == CollisionModel::EHSS)) {
            throw std::runtime_error("force_ou_for_stochastic cannot be true when using stochastic collision models (HSMC or EHSS)");
        }
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_PHYSICS_CONFIG_H

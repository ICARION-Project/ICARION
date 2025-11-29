// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file BoundaryAction.h
 * @brief Interface for ion-boundary interaction actions
 * 
 * When an ion hits a domain boundary, a BoundaryAction determines what happens:
 * - Absorption: ion deactivated (removed from simulation)
 * - Specular Reflection: elastic bounce (v' = v - 2(v·n)n)
 * - Diffuse Reflection: random cosine-weighted direction
 * - Thermal Reflection: Maxwell-Boltzmann re-emission at wall temperature
 * - Transmission: pass through (multi-domain transitions)
 * 
 * Used by DomainManager when boundary crossing detected.
 */

#pragma once

#include "core/types/IonState.h"
#include "core/types/Vec3.h"
#include <string>
#include <memory>

namespace ICARION {
namespace integrator {

/**
 * @brief Abstract interface for boundary interaction actions
 * 
 * Each concrete action implements apply() to modify ion state.
 */
class BoundaryAction {
public:
    virtual ~BoundaryAction() = default;
    
    /**
     * @brief Apply boundary action to ion
     * 
     * @param ion Ion state (modified in-place)
     * @param normal Surface normal vector (unit length, pointing inward)
     * @param boundary_pos Position where ion hit boundary [m]
     * @param temperature_K Wall temperature [K] (for thermal actions)
     * 
     * Implementation modifies:
     * - ion.vel (reflection, thermal re-emission)
     * - ion.pos (set to boundary_pos for accuracy)
     * - ion.active (false for absorption)
     */
    virtual void apply(
        IonState& ion,
        const Vec3& normal,
        const Vec3& boundary_pos,
        double temperature_K
    ) = 0;
    
    /**
     * @brief Get action name for logging
     */
    virtual std::string name() const = 0;
};

}  // namespace integrator
}  // namespace ICARION

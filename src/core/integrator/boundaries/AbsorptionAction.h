// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file AbsorptionAction.h
 * @brief Boundary action: absorb ion (remove from simulation)
 * 
 * Used for:
 * - Detectors (ions absorbed at detector surface)
 * - Lossy walls (no reflection, ion captured)
 * - Domain exits without transmission
 */

#pragma once

#include "BoundaryAction.h"

namespace ICARION {
namespace integrator {

/**
 * @brief Absorption: ion stops at boundary and deactivates
 * 
 * Physics:
 * - Ion position set to boundary intersection
 * - Velocity set to zero (stopped)
 * - Active flag set to false (removed from simulation)
 * 
 * Use cases:
 * - Detector plates (MCP, Faraday cup)
 * - Absorbing walls (lossy boundaries)
 * - Ion traps where ions escape
 */
class AbsorptionAction : public BoundaryAction {
public:
    void apply(
        IonState& ion,
        const Vec3& normal,
        const Vec3& boundary_pos,
        double temperature_K
    ) override {
        // Set position to boundary intersection (for accurate trajectory output)
        ion.pos = boundary_pos;
        
        // Stop ion (absorbed)
        ion.vel = Vec3{0.0, 0.0, 0.0};
        
        // Deactivate (remove from simulation)
        ion.active = false;
    }
    
    std::string name() const override {
        return "Absorption";
    }
};

}  // namespace integrator
}  // namespace ICARION

// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file BoundaryActionFactory.h
 * @brief Factory for creating boundary action instances
 */

#pragma once

#include "BoundaryAction.h"
#include "AbsorptionAction.h"
#include "ReflectionAction.h"
#include "core/config/types/BoundaryConfig.h"
#include <memory>
#include <random>

namespace ICARION {
namespace integrator {

/**
 * @brief Factory for creating boundary actions from configuration
 */
class BoundaryActionFactory {
public:
    /**
     * @brief Create boundary action from config
     * 
     * @param config Boundary configuration
     * @param rng Random number generator (required for diffuse/thermal reflection)
     * @return Unique pointer to boundary action
     */
    static std::unique_ptr<BoundaryAction> create(
        const config::BoundaryConfig& config,
        std::mt19937* rng = nullptr
    ) {
        using Type = config::BoundaryActionType;
        
        switch (config.type) {
            case Type::Absorption:
                return std::make_unique<AbsorptionAction>();
                
            case Type::SpecularReflection:
                return std::make_unique<ReflectionAction>(
                    ReflectionAction::Type::SPECULAR,
                    config.accommodation_coeff,
                    rng
                );
                
            case Type::DiffuseReflection:
                if (!rng) {
                    throw std::runtime_error("DiffuseReflection requires RNG");
                }
                return std::make_unique<ReflectionAction>(
                    ReflectionAction::Type::DIFFUSE,
                    config.accommodation_coeff,
                    rng
                );
                
            case Type::ThermalReflection:
                if (!rng) {
                    throw std::runtime_error("ThermalReflection requires RNG");
                }
                return std::make_unique<ReflectionAction>(
                    ReflectionAction::Type::THERMAL,
                    config.accommodation_coeff,
                    rng
                );
                
            default:
                throw std::runtime_error("Unknown boundary action type");
        }
    }
};

}  // namespace integrator
}  // namespace ICARION

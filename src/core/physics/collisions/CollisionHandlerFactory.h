// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
/**
 * @file CollisionHandlerFactory.h
 * @brief Factory for creating collision handlers based on collision model
 * 
 * Creates appropriate ICollisionHandler instances based on the collision model
 * specified in PhysicsConfig. Returns nullptr for deterministic models that
 * use DampingForce instead.
 * 
 * **SSOT Design:**
 * - Takes PhysicsConfig reference directly (no parameter copies)
 * - Factory validates requirements (e.g., EHSS needs geometry map)
 * - Handlers maintain SSOT principle by taking EnvironmentConfig references
 * 
 * **Collision Model Routing:**
 * - Stochastic (ICollisionHandler): EHSS, HSS
 * - Deterministic (DampingForce): Friction, Langevin, HSD
 * - NoCollisions → nullptr (no collision handling)
 * - OU thermalization handled separately via enable_ou_thermalization flag
 * 
 * **Usage:**
 * ```cpp
 * auto handler = CollisionHandlerFactory::create(
 *     config.physics,
 *     &geometry_map,
 *     gamma_coefficient
 * );
 * 
 * if (handler) {
 *     handler->handle_collision(ion, dt, rng, config.environment);
 * }
 * ```
 * 
 * @date 2025-11-21
 * @version 1.0
 */

#pragma once

#include "ICollisionHandler.h"
#include "EHSSCollisionHandler.h"
#include "HSSCollisionHandler.h"
#include "OUCollisionHandler.h"
#include "GPUCollisionHandler.h"
#include "core/config/types/PhysicsConfig.h"
#include "core/config/types/SpeciesConfig.h"
#include <memory>
#include <stdexcept>

namespace ICARION::physics {

/**
 * @brief Factory for creating collision handlers
 * 
 * Instantiates the appropriate ICollisionHandler subclass based on
 * the collision model specified in PhysicsConfig.
 * 
 * **Returns nullptr for:**
 * - NoCollisions (no collision handling needed)
 * - Friction, Langevin, HSD (use DampingForce instead)
 * 
 * **Validation:**
 * - EHSS requires geometry_map (throws if nullptr)
 * - OU requires gamma_for_ou > 0 (throws if invalid)
 */
class CollisionHandlerFactory {
public:
    /**
     * @brief Create collision handler based on configuration
     * 
     * @param config Physics configuration (SSOT reference)
     * @param geometry_map Molecular geometry map (required for EHSS, optional otherwise)
     * @param gamma_for_ou Damping coefficient for OU thermalization (required if enable_ou_thermalization)
     * @param enable_logging Enable collision logging (default: false)
     * 
     * @return Unique pointer to handler, or nullptr for deterministic/no-collision models
     * 
     * @throws std::invalid_argument if EHSS requested but geometry_map is nullptr
     * @throws std::invalid_argument if OU thermalization enabled but gamma_for_ou <= 0
     */
    static std::unique_ptr<ICollisionHandler> create(
        const config::PhysicsConfig& config,
        const GeometryMap* geometry_map = nullptr,
        double gamma_for_ou = 0.0,
        bool enable_logging = false,
        const config::SpeciesDatabase* species_db = nullptr,
        bool enable_gpu = false,
        unsigned long long gpu_seed = 42,
        size_t gpu_threshold = 5000
    );

private:
    CollisionHandlerFactory() = delete;  // Static class, no instantiation
};

} // namespace ICARION::physics

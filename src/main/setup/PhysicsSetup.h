// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file PhysicsSetup.h
 * @brief Physics module initialization helper for main.cpp
 * 
 * Creates force registries, integration strategy, collision/reaction handlers.
 */

#pragma once

#include "core/config/types/FullConfig.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/physics/forces/ForceRegistry.h"
#include "core/integrator/strategies/IIntegrationStrategy.h"
#include "core/physics/collisions/ICollisionHandler.h"
#include "core/physics/reactions/IReactionHandler.h"
#include "core/io/fieldArrayLoader.h"
#include "core/config/types/AnalyticalFieldModel.h"
#include "core/config/types/FieldProviderModel.h"
#include <vector>
#include <memory>

namespace ICARION::setup {

/**
 * @brief Physics module initialization result
 * 
 * Contains all physics dependencies needed by SimulationEngine.
 */
struct PhysicsModules {
    /// Force registries (one per domain)
    std::vector<std::shared_ptr<physics::ForceRegistry>> force_registries;
    
    /// Integration strategy (RK4, RK45, or Boris)
    std::shared_ptr<integrator::IIntegrationStrategy> integrator;
    
    /// Collision handler (EHSS, HSS, OU, or Friction)
    std::shared_ptr<physics::ICollisionHandler> collision_handler;
    
    /// Reaction handler (stochastic or deterministic)
    std::shared_ptr<physics::IReactionHandler> reaction_handler;
};

/**
 * @brief Physics setup helper
 * 
 * Handles physics module initialization with auto-configuration:
 * - Space charge method selection (Direct vs Grid based on N)
 * - Molecular geometry loading for EHSS
 * - Force registry per domain with appropriate forces
 */
class PhysicsSetup {
public:
    /**
     * @brief Initialize all physics modules
     * 
     * Creates:
     * - ForceRegistry per domain with:
 *   * ElectricFieldForce (always)
 *   * MagneticFieldForce (if B > 0)
 *   * DampingForce (if Friction model)
 *   * Space-charge model provided via SpaceChargeModelFactory (Direct/Grid/GPU)
     * - Integration strategy (RK4/RK45/Boris from config)
     * - Collision handler (EHSS/HSS/OU/Friction from config)
     * - Reaction handler (if reactions enabled)
     * 
     * @param config Simulation configuration (SSOT)
     * @param ions Initial ion ensemble (for space charge grid estimation)
     * @return Physics modules ready for SimulationEngine
     * 
     * @note Space charge method selection is delegated to SpaceChargeModelFactory,
     *       which evaluates domain geometry, ion count, and GPU availability.
     */
    static PhysicsModules initialize(
        const config::FullConfig& config,
        const core::IonEnsemble& ions
    );

private:
    /**
     * @brief Create force registries (one per domain)
     * 
     * Adds fundamental forces to each domain:
     * - ElectricFieldForce (always) + matching FieldModel (analytical or grid)
     * - MagneticFieldForce (if B > 0)
     * - DampingForce (if Friction model)
     *
     * Also loads field arrays from HDF5 if field_array_terms are present and
     * wraps them in FieldProviderModel for SSOT field access.
     */
    static std::vector<std::shared_ptr<physics::ForceRegistry>> create_force_registries(
        const config::FullConfig& config
    );
    
    /// Static storage for field arrays (must persist for simulation lifetime)
    /// Using unique_ptr to prevent pointer invalidation on vector reallocation
    static inline std::vector<std::unique_ptr<FieldArray>> field_arrays_storage_;
    /// Static storage for field models (must persist for simulation lifetime)
    static inline std::vector<std::unique_ptr<config::IFieldModel>> field_models_storage_;
    
    /**
     * @brief Add space charge models to registries
     * 
     * Delegates selection to SpaceChargeModelFactory for each domain. Direct models
     * are shared between registries, while grid/GPU models are instantiated per-domain.
     */
    static void add_space_charge_forces(
        std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
        const config::FullConfig& config,
        const core::IonEnsemble& ions
    );
    
    /**
     * @brief Create integration strategy from config
     * 
     * Supported: RK4, RK45, Boris
     * Defaults to RK45 if unknown.
     */
    static std::shared_ptr<integrator::IIntegrationStrategy> create_integrator(
        const config::FullConfig& config
    );
    
    /**
     * @brief Create collision handler from config
     * 
     * Loads molecular geometries for EHSS if needed.
     * Falls back to HSS if geometry loading fails.
     */
    static std::shared_ptr<physics::ICollisionHandler> create_collision_handler(
        const config::FullConfig& config
    );
    
    /**
     * @brief Create reaction handler from config
     */
    static std::shared_ptr<physics::IReactionHandler> create_reaction_handler(
        const config::FullConfig& config
    );
};

}  // namespace ICARION::setup

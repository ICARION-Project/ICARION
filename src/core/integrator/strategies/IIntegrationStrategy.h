// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#pragma once

#include <string>
#include <vector>
#include <memory>

#include "core/types/IonEnsemble.h"
#include "core/config/types/DomainConfig.h"
#include "core/physics/forces/ForceRegistry.h"

namespace ICARION {
namespace integrator {

/**
 * @brief Integration strategy interface (SoA-only)
 *
 * Defines the contract for numerical integration methods (RK4, RK45, Boris, etc.)
 * using the SoA ion layout. Legacy AoS entry points were removed to avoid
 * parallel code paths.
 * 
 * **Design Principles:**
 * - Strategy Pattern: Swap integrators without changing client code
 * - SSOT Compliance: Uses DomainConfig directly (no parameter conversion)
 * - Dependency Injection: ForceRegistry passed as parameter
 * - Zero-cost Abstraction: Virtual call overhead negligible vs. computation
 * 
 * **Example Usage:**
 * ```cpp
 * // Create strategy
 * auto strategy = std::make_unique<RK4Strategy>();
 * 
 * // Create force registry
 * physics::ForceRegistry registry;
 * registry.add_force(std::make_unique<ElectricFieldForce>(domain));
 * registry.add_force(std::make_unique<MagneticFieldForce>(domain.fields.magnetic));
 * 
 * // Integrate one timestep
 * strategy->step(ensemble, idx, t, dt, registry);
 * ```
 */
class IIntegrationStrategy {
public:
    virtual ~IIntegrationStrategy() = default;
    
    /**
     * @brief Advance ion state by one timestep (SoA)
     *
     * @param ensemble Ion ensemble (SoA)
     * @param ion_idx Index of ion to integrate
     * @param t Current simulation time [s]
     * @param dt Timestep size [s]
     * @param force_registry Force computation engine
     *
     * **Responsibilities:**
     * - Compute intermediate stages (k1, k2, ... for RK methods)
     * - Update ion position and velocity
     * - Ensure numerical stability
     *
     * **Not Responsible For:**
     * - Boundary checks (handled by SimulationEngine)
     * - Collision events (handled by CollisionHandler)
     * - Reaction events (handled by ReactionHandler)
     * - Output writing (handled by OutputManager)
     *
     * **Thread Safety:**
     * - Read-only access to force_registry and other ions
     * - Modifies only the selected ion (caller must ensure thread safety)
     */
    virtual void step(
        core::IonEnsemble& ensemble,
        size_t ion_idx,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry
    ) = 0;
    
    /**
     * @brief Get human-readable strategy name
     * 
     * @return Strategy identifier (e.g., "RK4", "RK45", "Boris")
     * 
     * Used for logging, debugging, and output file metadata.
     */
    virtual std::string name() const = 0;
    
    /**
     * @brief Check if strategy supports adaptive timestepping
     * 
     * @return true if adaptive (RK45, DOPRI5), false if fixed (RK4, Boris)
     * 
     * Adaptive strategies can implement IAdaptiveIntegrationStrategy
     * to provide error estimation and timestep control.
     */
    virtual bool is_adaptive() const = 0;

    /**
     * @brief Whether the strategy provides a batch integration shortcut.
     *
     * Default implementation returns false and signals that SimulationEngine
     * should call `step` for every ion. GPU-backed strategies override this
     * to run a domain-aware batch.
     */
    virtual bool supports_batch() const { return false; }

    /**
     * @brief Optional batch integration hook.
     *
     * @param ensemble Ion ensemble (SoA) to integrate in-place.
     * @param t Current time.
     * @param dt Timestep.
     * @param registries Force registries (per domain).
     * @param domain_indices Cached domain index per ion (-1 for skipped ions).
     * @return true when the batch path ran (GPU, etc.), false if the caller
     *         should fall back to per-ion integration.
     */
    virtual bool step_batch(
        core::IonEnsemble& ensemble,
        double t,
        double dt,
        const std::vector<std::shared_ptr<physics::ForceRegistry>>& registries,
        const std::vector<int>& domain_indices
    ) {
        (void)ensemble;
        (void)t;
        (void)dt;
        (void)registries;
        (void)domain_indices;
        return false;
    }

    /**
     * @brief Enable or disable parallel execution for batch paths.
     *
     * Default implementation is a no-op.
     */
    virtual void set_parallel_enabled(bool enabled) { (void)enabled; }
};

} // namespace integrator
} // namespace ICARION

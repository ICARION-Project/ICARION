// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <string>

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
};

} // namespace integrator
} // namespace ICARION

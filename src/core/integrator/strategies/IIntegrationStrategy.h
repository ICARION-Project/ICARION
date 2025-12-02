// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <string>
#include <vector>

#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"  // For SoA integration
#include "core/config/types/DomainConfig.h"
#include "core/physics/forces/ForceRegistry.h"

namespace ICARION {
namespace integrator {

/**
 * @brief Integration strategy interface
 * 
 * Defines contract for numerical integration methods (RK4, RK45, Boris, etc.).
 * Replaces legacy integrate_one_step() with modular, testable design.
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
 * strategy->step(ion, t, dt, registry, domain, all_ions);
 * ```
 */
class IIntegrationStrategy {
public:
    virtual ~IIntegrationStrategy() = default;
    
    /**
     * @brief Advance ion state by one timestep
     * 
     * @param ion Ion state (position, velocity, mass, charge) [in/out]
     *            Updated in-place with new position/velocity
     * @param t Current simulation time [s]
     * @param dt Timestep size [s]
     * @param force_registry Force computation engine (knows its domain via domain() method)
     * @param all_ions All ion states at current time (for space charge)
     * 
     * **SSOT Compliance (Phase 12 Enhancement):**
     * - `force_registry`: ForceRegistry now stores DomainConfig internally
     * - No need to pass domain separately (eliminates parameter duplication!)
     * - Domain accessible via force_registry.domain() if needed
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
     * - Read-only access to force_registry, all_ions
     * - Modifies only `ion` parameter (caller must ensure thread safety)
     */
    virtual void step(
        IonState& ion,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry,
        const std::vector<IonState>& all_ions
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
     * @brief Advance single ion using SoA (Structure of Arrays) data layout
     * 
     * Phase 3B: Cache-optimized integration using direct array access.
     * Default implementation converts to IonState and calls step().
     * 
     * @param ensemble Ion ensemble (SoA)
     * @param ion_idx Index of ion to integrate
     * @param t Current simulation time [s]
     * @param dt Timestep size [s]
     * @param force_registry Force computation engine
     * 
     * **Performance Notes:**
     * - Override for true SoA speedups; default wrapper converts the entire
     *   ensemble to AoS and back and is intended as a correctness fallback.
     * - Direct array access can reduce cache misses when implemented by a
     *   strategy that keeps computation in SoA form.
     */
    virtual void step_soa(
        core::IonEnsemble& ensemble,
        size_t ion_idx,
        double t,
        double dt,
        const physics::ForceRegistry& force_registry
    ) {
        // Default: convert to IonState and call legacy method
        auto ions_temp = ensemble.to_legacy();
        step(ions_temp[ion_idx], t, dt, force_registry, ions_temp);
        
        // Write back
        ensemble.set_pos(ion_idx, ions_temp[ion_idx].pos);
        ensemble.set_vel(ion_idx, ions_temp[ion_idx].vel);
    }
};

} // namespace integrator
} // namespace ICARION

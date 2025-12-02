// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#pragma once

#include <string>
#include "core/types/IonState.h"
#include "core/types/Vec3.h"

namespace ICARION::physics {

// Forward declarations
struct ForceContext;

/**
 * @brief Force contribution interface
 * 
 * All forces (electric, magnetic, damping, space charge) implement this interface.
 * Enables modular force composition via ForceRegistry.
 * 
 * Design Principles:
 * - Single Responsibility: Each force computes ONE physical contribution
 * - Const-correctness: Forces don't modify ion state (pure computation)
 * - Composability: Forces can be combined via ForceRegistry
 * - Testability: Each force can be unit tested in isolation
 * 
 * Example Usage:
 * @code
 * // Create electric field force
 * auto e_force = std::make_unique<ElectricFieldForce>(domain, field_provider);
 * 
 * // Compute force on ion
 * ForceContext ctx{};
 * Vec3 force = e_force->compute(ion, t, ctx);
 * @endcode
 * 
 * @see ForceRegistry for force composition
 * @see ForceContext for shared context data
 */
class IForce {
public:
    virtual ~IForce() = default;
    
    /**
     * @brief Compute force contribution for single ion
     * 
     * @param ion Current ion state (position, velocity, charge, mass)
     * @param t Current simulation time [s]
     * @param context Optional context for force computation (field provider, all ions, etc.)
     * @return Force vector [N] in Cartesian coordinates
     * 
     * @note Should be thread-safe if used under OpenMP parallel loops
     */
    virtual Vec3 compute(
        const IonState& ion,
        double t,
        const ForceContext& context
    ) const = 0;
    
    /**
     * @brief Check if this force applies to given ion
     * 
     * Allows conditional force application (e.g., magnetic force only if B-field enabled).
     * Default implementation returns true (force always applies).
     * 
     * @param ion Ion to check
     * @return true if force should be computed for this ion
     * 
     * @note This is called before compute() to avoid unnecessary calculations
     */
    virtual bool applies_to(const IonState& ion) const {
        (void)ion;  // Unused parameter
        return true;
    }
    
    /**
     * @brief Get force name for logging/debugging
     * 
     * @return Human-readable force name (e.g., "ElectricField", "MagneticField")
     */
    virtual std::string name() const = 0;
};

} // namespace ICARION::physics

// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#pragma once

#include <memory>
#include <vector>
#include "IForce.h"
#include "ForceContext.h"
#include "core/config/types/DomainConfig.h"

namespace ICARION::physics {

/**
 * @brief Manages and aggregates all active forces for simulation
 * 
 * ForceRegistry implements the Composite pattern to enable modular force composition.
 * Forces can be added/removed dynamically without modifying the integrator.
 * 
 * Key Features:
 * - Dynamic force registration (add/remove at runtime)
 * - Automatic force aggregation (sum of all contributions)
 * - Conditional force application (via IForce::applies_to())
 * - Zero-cost abstraction (inline-friendly, no virtual call overhead in hot path)
 * 
 * Typical Usage:
 * @code
 * // Setup forces
 * ForceRegistry registry;
 * registry.add_force(std::make_unique<ElectricFieldForce>(domain));
 * registry.add_force(std::make_unique<MagneticFieldForce>(mag_config));
 * registry.add_force(std::make_unique<DampingForce>(gamma));
 * 
 * // In integration loop
 * ForceContext ctx{field_provider, domain, all_ions};
 * Vec3 total_force = registry.compute_total_force(ion, t, ctx);
 * Vec3 acceleration = total_force / ion.mass_kg;
 * @endcode
 * 
 * Thread Safety:
 * - compute_total_force() is const and thread-safe (read-only)
 * - add_force() / clear() are NOT thread-safe (call during setup only)
 * 
 * @see IForce for force interface
 * @see ForceContext for shared context data
 */
class ForceRegistry {
public:
    /**
     * @brief Default constructor (empty registry, no domain)
     * @deprecated Use ForceRegistry(const config::DomainConfig&) instead
     */
    ForceRegistry() = default;
    
    /**
     * @brief Construct registry with domain context (RECOMMENDED)
     * 
     * @param domain Domain configuration (geometry, fields, environment)
     * 
     * This constructor allows forces to access domain-specific context
     * without needing it passed through every method call.
     * 
     * Example:
     * @code
     * ForceRegistry registry(domain_config);
     * registry.add_force(std::make_unique<ElectricFieldForce>());
     * // Force can now access domain via registry.domain()
     * @endcode
     */
    explicit ForceRegistry(const config::DomainConfig& domain)
        : domain_(&domain) {}
    
    /**
     * @brief Add force to registry
     * 
     * Transfers ownership of the force to the registry.
     * Forces are evaluated in the order they are added.
     * 
     * @param force Force implementation (ownership transferred)
     * 
     * Example:
     * @code
     * registry.add_force(std::make_unique<ElectricFieldForce>(domain));
     * @endcode
     * 
     * @note NOT thread-safe (call during setup phase only)
     */
    void add_force(std::unique_ptr<IForce> force);
    
    /**
     * @brief Compute total force on ion (sum of all registered forces)
     * 
     * This is the hot-path method called every integration step.
     * Optimized for performance:
     * - Early exit if no forces registered
     * - Skips forces that don't apply (via applies_to())
     * - Inline-friendly design
     * 
     * @param ion Current ion state
     * @param t Current simulation time [s]
     * @param context Optional context for force computation
     * @return Total force vector [N]
     * 
     * @note Thread-safe (const method, read-only access)
     */
    Vec3 compute_total_force(
        const IonState& ion,
        double t,
        const ForceContext& context = {}
    ) const;
    
    /**
     * @brief Get all registered forces (for inspection/debugging)
     * 
     * @return Vector of force pointers (const reference, no ownership transfer)
     */
    const std::vector<std::unique_ptr<IForce>>& forces() const {
        return forces_;
    }
    
    /**
     * @brief Get number of registered forces
     * 
     * @return Number of active forces
     */
    size_t size() const {
        return forces_.size();
    }
    
    /**
     * @brief Check if registry is empty
     * 
     * @return true if no forces registered
     */
    bool empty() const {
        return forces_.empty();
    }
    
    /**
     * @brief Clear all forces (removes all registered forces)
     * 
     * Useful for resetting simulation or switching force configurations.
     * 
     * @note NOT thread-safe (call during setup phase only)
     */
    void clear() {
        forces_.clear();
    }
    
    /**
     * @brief Get domain configuration (if available)
     * 
     * @return Pointer to domain config, or nullptr if not set
     * 
     * Allows forces to access domain context (geometry, fields, environment).
     */
    const config::DomainConfig* domain() const {
        return domain_;
    }
    
private:
    /**
     * @brief Vector of registered forces (owned by registry)
     */
    std::vector<std::unique_ptr<IForce>> forces_;
    
    /**
     * @brief Domain configuration (optional, nullptr if not set)
     * 
     * Non-owning pointer (registry does not own domain config).
     * Must remain valid for lifetime of registry.
     */
    const config::DomainConfig* domain_ = nullptr;
};

} // namespace ICARION::physics

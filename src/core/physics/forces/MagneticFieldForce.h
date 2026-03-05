// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
/**
 * @file MagneticFieldForce.h
 * @brief Magnetic field force (Lorentz force) implementation
 * 
 * Computes Lorentz magnetic force F = q·(v × B) for ions in magnetic fields.
 * Supports uniform/gradient fields from config or a field provider. Field providers
 * must supply B via the shared interface (see note below).
 */

#pragma once

#include "IForce.h"
#include "ForceContext.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"

#include <memory>

// Forward declaration (IFieldProvider is in global namespace)
class IFieldProvider;

namespace ICARION {

// Forward declaration for config types
namespace config {
    struct MagneticFieldConfig;
}

namespace physics {

/**
 * @class MagneticFieldForce
 * @brief Computes magnetic Lorentz force F = q·(v × B)
 * 
 * Two modes of operation:
 * 
 * 1. **Analytical Mode** (uniform + gradient):
 *    - B(r) = B₀ + ∇B·r (linear model)
 *    - Fast, suitable for FTICR, Orbitrap with uniform B-field
 * 
 * 2. **Field Provider Mode** (grid/BEM/FEM):
 *    - Uses IFieldProvider for spatially-varying B-fields
 *    - Required for complex magnet geometries
 * 
 * **Usage:**
 * ```cpp
 * // Analytical mode (SSOT-compliant)
 * config::DomainConfig domain = load_config("config.json");
 * domain.fields.magnetic.enabled = true;
 * domain.fields.magnetic.field_strength_T = 7.0;  // 7T along z-axis
 * auto force = std::make_unique<MagneticFieldForce>(domain.fields.magnetic);
 * 
 * // Field provider mode (spatially-varying field)
 * auto b_field_provider = std::make_shared<GridFieldProvider>(...);
 * auto force = std::make_unique<MagneticFieldForce>(b_field_provider);
 * ```
 * 
 * @note Magnetic force is velocity-dependent: F = q·(v × B)
 * @note For FTICR: radial confinement from B-field, axial from E-field
 * @note For Orbitrap: B-field typically negligible (electrostatic trap)
 */
class MagneticFieldForce : public IForce {
public:
    /**
     * @brief Construct from magnetic field configuration (SSOT)
     * @param magnetic_config Reference to magnetic field configuration
     * 
     * ⚠️ Config reference must outlive this object!
     */
    explicit MagneticFieldForce(const ICARION::config::MagneticFieldConfig& magnetic_config);
    
    /**
     * @brief Construct from field provider (grid/BEM/FEM)
     * @param field_provider Field provider for B-field evaluation
     * 
     * @note Field provider must implement get_E() for B-field
     *       (reusing IFieldProvider interface, E→B semantic mapping)
     */
    explicit MagneticFieldForce(std::shared_ptr<::IFieldProvider> field_provider);
    
    /**
     * @brief Compute magnetic force F = q·(v × B)
     * 
     * @param ion Ion state (velocity required!)
     * @param t Current simulation time [s]
     * @param ctx Force context (field provider, domain config)
     * @return Force vector [N]
     * 
     * Uses field provider if available, otherwise analytical formula.
     */
    Vec3 compute(const core::IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override;

    Vec3 compute_soa(const ForceState& state, double t,
                     const ForceContext& ctx) const override;
    
    /**
     * @brief Get force name
     * @return "MagneticField" or "MagneticField(FieldProvider)"
     */
    std::string name() const override;

private:
    /**
     * @brief Compute analytical magnetic field B(r)
     * @param pos Position [m]
     * @return Magnetic field vector [T]
     */
    Vec3 compute_analytical_field(const Vec3& pos) const;
    
    // Field calculation mode
    bool use_field_provider_ = false;
    std::shared_ptr<::IFieldProvider> field_provider_ = nullptr;
    const ICARION::config::MagneticFieldConfig* magnetic_config_ = nullptr;  // SSOT: config reference
};

} // namespace physics
} // namespace ICARION

// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file MagneticFieldForce.h
 * @brief Magnetic field force (Lorentz force) implementation
 * 
 * Computes Lorentz magnetic force F = q·(v × B) for ions in magnetic fields.
 * Supports both uniform fields and field providers (grid-based, gradient fields).
 */

#pragma once

#include "IForce.h"
#include "ForceContext.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"

#include <memory>

// Forward declaration (IFieldProvider is in global namespace)
class IFieldProvider;

namespace ICARION {
namespace physics {

/**
 * @brief Parameters for analytical magnetic field configuration
 */
struct MagneticFieldParams {
    Vec3 uniform_field_T = {0.0, 0.0, 0.0};     ///< Uniform B-field [Tesla]
    Vec3 gradient_T_per_m = {0.0, 0.0, 0.0};    ///< Linear gradient [T/m]
    bool enabled = false;                        ///< Enable magnetic force
};

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
 * // Analytical mode (uniform B-field for FTICR)
 * MagneticFieldParams params;
 * params.uniform_field_T = {0.0, 0.0, 7.0};  // 7T along z-axis
 * params.enabled = true;
 * auto force = std::make_unique<MagneticFieldForce>(params);
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
     * @brief Construct from analytical field parameters
     * @param params Magnetic field configuration (uniform + gradient)
     */
    explicit MagneticFieldForce(const MagneticFieldParams& params);
    
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
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override;
    
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
    MagneticFieldParams analytical_params_;
};

} // namespace physics
} // namespace ICARION

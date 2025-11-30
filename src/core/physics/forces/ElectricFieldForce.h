// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

/**
 * @file ElectricFieldForce.h
 * @brief Electric field force implementation for all instrument types
 * 
 * Computes Lorentz electric force F = q·E for ions in electric fields.
 * Supports both analytical field calculations (instrument-specific) and
 * field provider-based evaluation (interpolated from grid data).
 */

#pragma once

#include "IForce.h"
#include "ForceContext.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "core/config/types/DomainConfig.h"

#include <memory>

// Forward declaration (IFieldProvider is in global namespace)
class IFieldProvider;

namespace ICARION {

// Forward declaration for config types
namespace config {
    struct DomainConfig;
}

namespace physics {

/**
 * @class ElectricFieldForce
 * @brief Computes electric field force F = q·E
 * 
 * Two modes of operation:
 * 
 * 1. **Analytical Mode** (no field provider):
 *    - Uses instrument-specific analytical field formulas
 *    - Fast, exact for ideal geometries
 *    - Supports: IMS, LQIT, TOF, FTICR, Orbitrap, QuadrupoleRF
 * 
 * 2. **Field Provider Mode** (with field provider):
 *    - Uses IFieldProvider to evaluate E-field at ion position
 *    - Supports grid-based interpolation, BEM, FEM, etc.
 *    - More general but requires precomputed field data
 * 
 * **Usage:**
 * ```cpp
 * // Analytical mode (SSOT-compliant)
 * config::DomainConfig domain = load_config("config.json");
 * auto force = std::make_unique<ElectricFieldForce>(domain);
 * 
 * // Field provider mode
 * auto field_provider = std::make_unique<GridFieldProvider>(...);
 * auto force = std::make_unique<ElectricFieldForce>(std::move(field_provider));
 * ```
 * 
 * @note Reads fields directly from DomainConfig (SSOT).
 */
class ElectricFieldForce : public IForce {
public:
    /**
     * @brief Construct from domain configuration (SSOT)
     * @param domain Domain configuration (instrument, fields, geometry)
     * 
     * ⚠️ Config reference must outlive this object!
     */
    explicit ElectricFieldForce(const ICARION::config::DomainConfig& domain);
    
    /**
     * @brief Construct from field provider (grid/BEM/FEM)
     * @param field_provider Field provider for E-field evaluation
     */
    explicit ElectricFieldForce(std::shared_ptr<::IFieldProvider> field_provider);
    
    /**
     * @brief Compute electric force F = q·E
     * 
     * @param ion Ion state (position, charge)
     * @param t Current simulation time [s]
     * @param ctx Force context (field provider, domain config)
     * @return Force vector [N]
     * 
     * Uses field provider if available, otherwise analytical formulas.
     */
    Vec3 compute(const IonState& ion, double t, const ForceContext& ctx) const override;
    
    /**
     * @brief Get force name
     * @return "ElectricField" or "ElectricField(InstrumentType)"
     */
    std::string name() const override;

    /**
     * @brief Get field provider (for GPU field extraction)
     * @return Pointer to field provider or nullptr if using analytical fields
     * 
     * Used by GPU integration to extract field data for texture upload.
     */
    const ::IFieldProvider* get_field_provider() const {
        return field_provider_.get();
    }

private:
    /**
     * @brief Compute analytical electric field for instrument type
     * @param ion Ion state
     * @param t Current time [s]
     * @return Electric field vector [V/m]
     * 
     * Reads parameters directly from domain_ (SSOT).
     */
    Vec3 compute_analytical_field(const IonState& ion, double t) const;
    
    // Instrument-specific field calculations (all read from domain_)
    Vec3 compute_lqit_field(const IonState& ion, double t) const;
    Vec3 compute_ims_field(const IonState& ion, double t) const;  // Time-dependent (RF field)
    Vec3 compute_tof_field(const IonState& ion) const;
    Vec3 compute_fticr_field(const IonState& ion) const;
    Vec3 compute_orbitrap_field(const IonState& ion) const;
    Vec3 compute_quadrupole_rf_field(const IonState& ion, double t) const;
    
    // Field calculation mode
    bool use_field_provider_ = false;
    std::shared_ptr<::IFieldProvider> field_provider_ = nullptr;
    const ICARION::config::DomainConfig* domain_ = nullptr;  // SSOT: config reference
};

} // namespace physics
} // namespace ICARION

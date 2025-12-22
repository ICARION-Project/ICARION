// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

/**
 * @file ElectricFieldForce.h
 * @brief Electric field force implementation for all instrument types
 * 
 * Computes Lorentz electric force F = q·E for ions in electric fields.
 * Supports both analytical field calculations (instrument-specific) and
 * field provider-based evaluation (interpolated from grid data). Field
 * provider mode requires precomputed fields; analytical mode covers only the
 * implemented instrument types.
 */

#pragma once

#include "IForce.h"
#include "ForceContext.h"
#include "core/types/Vec3.h"
#include "core/types/IonState.h"
#include "core/types/IonEnsemble.h"
#include "core/config/types/DomainConfig.h"
#include "core/config/types/IFieldModel.h"
#include "core/config/types/AnalyticalFieldModel.h"

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
 * Prefers field sampling via IFieldModel (SSOT: injected by setup).
 * Falls back to a provider (legacy) or an internal AnalyticalFieldModel if
 * no model is available.
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
    Vec3 compute(const core::IonEnsemble& ensemble, size_t ion_idx, double t,
                 const ForceContext& ctx) const override;

    Vec3 compute_soa(const ForceState& state, double t,
                     const ForceContext& ctx) const override;
    
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
    // Field calculation mode
    bool use_field_provider_ = false;
    std::shared_ptr<::IFieldProvider> field_provider_ = nullptr;
    const ICARION::config::DomainConfig* domain_ = nullptr;  // SSOT: config reference

    // Local fallback model (analytical) if no model is injected
    std::unique_ptr<ICARION::config::IFieldModel> fallback_model_;
};

} // namespace physics
} // namespace ICARION

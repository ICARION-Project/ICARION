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
#include "instrument/InstrumentTypes.h"

#include <memory>

// Forward declaration (IFieldProvider is in global namespace)
class IFieldProvider;

namespace ICARION {
namespace physics {

/**
 * @brief Analytical field parameters for instrument-specific calculations
 * 
 * Different instruments use different subsets of these parameters.
 * Unused parameters can be left at default values (0.0).
 */
struct AnalyticalFieldParams {
    // Common parameters
    ICARION::instrument::InstrumentType instrument_type = 
        ICARION::instrument::InstrumentType::UnknownInstrument;
    double radius_m = 0.0;        ///< Trap/cell radius [m]
    double length_m = 0.0;        ///< Axial length [m]
    
    // RF field (LQIT, QuadrupoleRF)
    double rf_voltage_V = 0.0;           ///< RF amplitude [V]
    double rf_frequency_Hz = 0.0;        ///< RF frequency [Hz]
    double dc_quad_voltage_V = 0.0;      ///< DC quadrupole offset [V]
    
    // DC field (IMS, TOF, axial confinement)
    double dc_axial_voltage_V = 0.0;     ///< DC voltage along axis [V]
    
    // AC field (resonant excitation)
    double ac_voltage_V = 0.0;           ///< AC amplitude [V]
    double ac_frequency_Hz = 0.0;        ///< AC frequency [Hz]
    Vec3 ac_direction = {1.0, 0.0, 0.0}; ///< AC field direction (normalized internally)
    
    // Orbitrap-specific
    double orbitrap_k = 0.0;             ///< Field curvature [V/m²]
    double orbitrap_r_char = 0.0;        ///< Characteristic radius [m]
    
    // FTICR-specific
    double fticr_voltage_V = 0.0;        ///< Trapping voltage [V]
    double fticr_char_length_m = 0.0;    ///< Characteristic length [m]
};

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
 * // Analytical mode (LQIT)
 * AnalyticalFieldParams params;
 * params.instrument_type = InstrumentType::LQIT;
 * params.radius_m = 0.005;
 * params.rf_voltage_V = 1000.0;
 * params.rf_frequency_Hz = 1e6;
 * auto force = std::make_unique<ElectricFieldForce>(params);
 * 
 * // Field provider mode
 * auto field_provider = std::make_unique<GridFieldProvider>(...);
 * auto force = std::make_unique<ElectricFieldForce>(std::move(field_provider));
 * ```
 */
class ElectricFieldForce : public IForce {
public:
    /**
     * @brief Construct from analytical field parameters
     * @param params Instrument-specific field parameters
     */
    explicit ElectricFieldForce(const AnalyticalFieldParams& params);
    
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

private:
    /**
     * @brief Compute analytical electric field for instrument type
     * @param ion Ion state
     * @param t Current time [s]
     * @return Electric field vector [V/m]
     */
    Vec3 compute_analytical_field(const IonState& ion, double t) const;
    
    // Instrument-specific field calculations
    Vec3 compute_lqit_field(const IonState& ion, double t) const;
    Vec3 compute_ims_field(const IonState& ion) const;
    Vec3 compute_tof_field(const IonState& ion) const;
    Vec3 compute_fticr_field(const IonState& ion) const;
    Vec3 compute_orbitrap_field(const IonState& ion) const;
    Vec3 compute_quadrupole_rf_field(const IonState& ion, double t) const;
    
    // Field calculation mode
    bool use_field_provider_ = false;
    std::shared_ptr<::IFieldProvider> field_provider_ = nullptr;
    AnalyticalFieldParams analytical_params_;
};

} // namespace physics
} // namespace ICARION

// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "ElectricFieldForce.h"
#include "core/utils/mathUtils.h"
#include "fieldsolver/utils/IFieldProvider.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace ICARION {
namespace physics {

// ============================================================================
// Constructors
// ============================================================================

ElectricFieldForce::ElectricFieldForce(const AnalyticalFieldParams& params)
    : use_field_provider_(false)
    , analytical_params_(params)
{
    // Validate instrument type
    if (params.instrument_type == ICARION::instrument::InstrumentType::UnknownInstrument) {
        throw std::invalid_argument(
            "ElectricFieldForce: Unknown instrument type. "
            "Must specify valid instrument_type for analytical mode."
        );
    }
}

ElectricFieldForce::ElectricFieldForce(std::shared_ptr<::IFieldProvider> field_provider)
    : use_field_provider_(true)
    , field_provider_(std::move(field_provider))
{
    if (!field_provider_) {
        throw std::invalid_argument(
            "ElectricFieldForce: Field provider is null. "
            "Use analytical constructor if field provider unavailable."
        );
    }
}

// ============================================================================
// IForce Interface Implementation
// ============================================================================

Vec3 ElectricFieldForce::compute(const IonState& ion, double t, const ForceContext& ctx) const {
    Vec3 E_field{0.0, 0.0, 0.0};
    
    // Priority: context field provider > constructor field provider > analytical
    if (ctx.field_provider) {
        E_field = ctx.field_provider->get_E(ion.pos);
    } else if (use_field_provider_ && field_provider_) {
        E_field = field_provider_->get_E(ion.pos);
    } else {
        E_field = compute_analytical_field(ion, t);
    }
    
    // F = q * E
    return E_field * ion.ion_charge_C;
}

std::string ElectricFieldForce::name() const {
    if (use_field_provider_) {
        return "ElectricField(FieldProvider)";
    }
    
    using Instrument = ICARION::instrument::InstrumentType;
    switch (analytical_params_.instrument_type) {
        case Instrument::LQIT:         return "ElectricField(LQIT)";
        case Instrument::IMS:          return "ElectricField(IMS)";
        case Instrument::TOF:          return "ElectricField(TOF)";
        case Instrument::FTICR:        return "ElectricField(FTICR)";
        case Instrument::Orbitrap:     return "ElectricField(Orbitrap)";
        case Instrument::QuadrupoleRF: return "ElectricField(QuadrupoleRF)";
        default:                       return "ElectricField(Unknown)";
    }
}

// ============================================================================
// Analytical Field Calculations
// ============================================================================

Vec3 ElectricFieldForce::compute_analytical_field(const IonState& ion, double t) const {
    using Instrument = ICARION::instrument::InstrumentType;
    
    switch (analytical_params_.instrument_type) {
        case Instrument::LQIT:
            return compute_lqit_field(ion, t);
        
        case Instrument::IMS:
            return compute_ims_field(ion);
        
        case Instrument::TOF:
            return compute_tof_field(ion);
        
        case Instrument::FTICR:
            return compute_fticr_field(ion);
        
        case Instrument::Orbitrap:
            return compute_orbitrap_field(ion);
        
        case Instrument::QuadrupoleRF:
            return compute_quadrupole_rf_field(ion, t);
        
        case Instrument::NoFixedInstrument:
            // Generic instrument - no analytical field
            return Vec3{0.0, 0.0, 0.0};
        
        default:
            throw std::runtime_error(
                "ElectricFieldForce: Unsupported instrument type for analytical field calculation"
            );
    }
}

// ----------------------------------------------------------------------------
// LQIT: Linear Quadrupole Ion Trap
// ----------------------------------------------------------------------------
// Quadrupolar RF field: E_x = (U_DC + U_RF·cos(ωt)) · 2x/r₀²
//                        E_y = -(U_DC + U_RF·cos(ωt)) · 2y/r₀²
//                        E_z = 0 (radial confinement only)
// 
// DC endcap confinement (axial):
//   - Linear restoring field in outer 10% of trap length
//   - E_z ∝ distance from edge (pushes ions inward from both sides)
//   - Approximation for DC endcap potentials on both trap ends
// ----------------------------------------------------------------------------
Vec3 ElectricFieldForce::compute_lqit_field(const IonState& ion, double t) const {
    Vec3 E_total{0.0, 0.0, 0.0};
    
    // RF + DC quadrupole field (radial)
    if (std::fabs(analytical_params_.rf_voltage_V) > 1e-12) {
        const double r0_sq = analytical_params_.radius_m * analytical_params_.radius_m;
        const double omega = 2.0 * M_PI * analytical_params_.rf_frequency_Hz;
        const double U_eff = analytical_params_.dc_quad_voltage_V 
                           + analytical_params_.rf_voltage_V * std::cos(omega * t);
        
        E_total.x =  2.0 * ion.pos.x * U_eff / r0_sq;
        E_total.y = -2.0 * ion.pos.y * U_eff / r0_sq;
    }
    
    // DC endcap field (axial confinement)
    // Linear restoring field in outer 10% of trap length
    // Approximates DC endcap potentials pushing ions inward from both sides
    if (std::fabs(analytical_params_.dc_axial_voltage_V) > 1e-12) {
        const double L = analytical_params_.length_m;
        const double z = ion.pos.z;
        const double edge_region = 0.1 * L;  // Outer 10%
        
        double E_z = 0.0;
        
        // Left endcap region (z < 0.1*L): E_z increases from 0 to max as z → 0
        if (z < edge_region) {
            // E_z > 0 (pushes right, away from left endcap)
            E_z = analytical_params_.dc_axial_voltage_V * (edge_region - z) / (edge_region * L);
        }
        // Right endcap region (z > 0.9*L): E_z decreases from 0 to -max as z → L
        else if (z > L - edge_region) {
            // E_z < 0 (pushes left, away from right endcap)
            E_z = -analytical_params_.dc_axial_voltage_V * (z - (L - edge_region)) / (edge_region * L);
        }
        // Central region: field-free (E_z = 0)
        
        E_total.z = E_z;
    }
    
    // AC field (resonant excitation): FIXED in x-direction for v1.0
    // Multi-domain rotation handling deferred to v1.1+
    if (std::fabs(analytical_params_.ac_voltage_V) > 1e-12) {
        const double omega_ac = 2.0 * M_PI * analytical_params_.ac_frequency_Hz;
        const double mag = (analytical_params_.ac_voltage_V / analytical_params_.radius_m) 
                         * std::cos(omega_ac * t);
        E_total.x += mag;  // Fixed x-direction
    }
    
    return E_total;
}

// ----------------------------------------------------------------------------
// IMS: Ion Mobility Spectrometry
// ----------------------------------------------------------------------------
// Uniform drift field: E_z = -U_drift / L
// Pushes ions along drift tube with constant axial field.
// ----------------------------------------------------------------------------
Vec3 ElectricFieldForce::compute_ims_field(const IonState& ion) const {
    (void)ion;  // Position-independent
    
    return Vec3{
        0.0,
        0.0,
        analytical_params_.dc_axial_voltage_V / analytical_params_.length_m
    };
}

// ----------------------------------------------------------------------------
// TOF: Time-of-Flight
// ----------------------------------------------------------------------------
// Similar to IMS - uniform drift field for acceleration/drift region.
// Can be field-free (E=0) or have extraction field.
// ----------------------------------------------------------------------------
Vec3 ElectricFieldForce::compute_tof_field(const IonState& ion) const {
    (void)ion;  // Position-independent
    
    return Vec3{
        0.0,
        0.0,
        analytical_params_.dc_axial_voltage_V / analytical_params_.length_m
    };
}

// ----------------------------------------------------------------------------
// FTICR: Fourier Transform Ion Cyclotron Resonance
// ----------------------------------------------------------------------------
// Quadrupolar trapping potential for axial confinement:
//   φ(r,z) = V_trap/d² · (2z² - x² - y²)
//   E_x = -∂φ/∂x = V_trap/d² · 2x
//   E_y = -∂φ/∂y = V_trap/d² · 2y
//   E_z = -∂φ/∂z = V_trap/d² · (-4z)
// Radial confinement from magnetic field (handled separately).
// ----------------------------------------------------------------------------
Vec3 ElectricFieldForce::compute_fticr_field(const IonState& ion) const {
    const double factor = analytical_params_.fticr_voltage_V 
                        / (analytical_params_.fticr_char_length_m 
                         * analytical_params_.fticr_char_length_m);
    
    // Axial center (for z-dependence)
    const double z_center = ion.pos.z - 0.5 * analytical_params_.length_m;
    
    return Vec3{
         factor * ion.pos.x,
         factor * ion.pos.y,
        -2.0 * factor * z_center
    };
}

// ----------------------------------------------------------------------------
// Orbitrap: Hyperlogarithmic Field
// ----------------------------------------------------------------------------
// Potential: U(r,z) = k/2 · (z² - r²/2 + r_char²·ln(r/r_char))
// Field: E = -∇U
//   E_x = -k/2 · x · (1 - r_char²/r²)
//   E_y = -k/2 · y · (1 - r_char²/r²)
//   E_z = -k · z
// Creates axial harmonic oscillation (mass-dependent frequency).
// ----------------------------------------------------------------------------
Vec3 ElectricFieldForce::compute_orbitrap_field(const IonState& ion) const {
    const double r_sq = std::max(1e-18, ion.pos.x * ion.pos.x + ion.pos.y * ion.pos.y);
    const double r_char_sq = analytical_params_.orbitrap_r_char * analytical_params_.orbitrap_r_char;
    const double C = 1.0 - r_char_sq / r_sq;
    
    // Center z-coordinate
    const double z_center = ion.pos.z - 0.5 * analytical_params_.length_m;
    
    return Vec3{
        0.5 * analytical_params_.orbitrap_k * ion.pos.x * C,
        0.5 * analytical_params_.orbitrap_k * ion.pos.y * C,
       -analytical_params_.orbitrap_k * z_center
    };
}

// ----------------------------------------------------------------------------
// QuadrupoleRF: Generic Quadrupole RF 
// ----------------------------------------------------------------------------
// Same as LQIT but may have different typical parameters.
// ----------------------------------------------------------------------------
Vec3 ElectricFieldForce::compute_quadrupole_rf_field(const IonState& ion, double t) const {
    // Identical to LQIT for now
    return compute_lqit_field(ion, t);
}

} // namespace physics
} // namespace ICARION

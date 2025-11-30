// SPDX-License-Identifier: MIT
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#include "ElectricFieldForce.h"
#include "core/utils/mathUtils.h"
#include "fieldsolver/utils/IFieldProvider.h"
#include "core/config/types/DomainConfig.h"

#include <cmath>
#include <algorithm>
#include <stdexcept>

namespace {
    constexpr double MIN_VOLTAGE_THRESHOLD = 1e-12;  ///< Minimum voltage [V] to consider field active (numerical safety)
    
    /**
     * @brief Helper to evaluate ValueOrWaveform at current time
     * @param val ValueOrWaveform (static or waveform)
     * @param t Current time [s]
     * @param library Waveform library for references
     * @return Evaluated value at time t
     */
    inline double eval_value(const ICARION::config::ValueOrWaveform& val, double t, 
                             const std::map<std::string, ICARION::config::Waveform>& library) {
        if (val.constant_value.has_value()) {
            return val.constant_value.value();
        }
        return val.evaluate(t, library);
    }
}

namespace ICARION {
namespace physics {

// ============================================================================
// Constructors
// ============================================================================

ElectricFieldForce::ElectricFieldForce(const ICARION::config::DomainConfig& domain)
    : use_field_provider_(false)
    , domain_(&domain)  // SSOT: store config reference
{
    // No validation needed - config can be any instrument type
}

ElectricFieldForce::ElectricFieldForce(std::shared_ptr<::IFieldProvider> field_provider)
    : use_field_provider_(true)
    , field_provider_(std::move(field_provider))
    , domain_(nullptr)  // Field provider mode: no config needed
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
        E_field = ctx.field_provider->get_E(ion.pos, t);  // Time-dependent evaluation
    } else if (use_field_provider_ && field_provider_) {
        E_field = field_provider_->get_E(ion.pos, t);      // Time-dependent evaluation
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
    
    if (!domain_) {
        return "ElectricField(NoConfig)";
    }
    
    // SSOT: Read instrument from config
    using Inst = ICARION::config::Instrument;
    switch (domain_->instrument) {
        case Inst::LQIT:              return "ElectricField(LQIT)";
        case Inst::IMS:               return "ElectricField(IMS)";
        case Inst::TOF:               return "ElectricField(TOF)";
        case Inst::FTICR:             return "ElectricField(FTICR)";
        case Inst::Orbitrap:          return "ElectricField(Orbitrap)";
        case Inst::QuadrupoleRF:      return "ElectricField(QuadrupoleRF)";
        case Inst::NoFixedInstrument: return "ElectricField(NoFixedInstrument)";
        case Inst::UnknownInstrument: return "ElectricField(Unknown)";
        default:                      return "ElectricField(Unknown)";
    }
}

// ============================================================================
// Analytical Field Calculations
// ============================================================================

Vec3 ElectricFieldForce::compute_analytical_field(const IonState& ion, double t) const {
    if (!domain_) {
        return Vec3{0.0, 0.0, 0.0};  // No config available
    }
    
    // SSOT: Read instrument from config
    using Inst = ICARION::config::Instrument;
    switch (domain_->instrument) {
        case Inst::LQIT:
            return compute_lqit_field(ion, t);
        
        case Inst::IMS:
            return compute_ims_field(ion, t);
        
        case Inst::TOF:
            return compute_tof_field(ion);
        
        case Inst::FTICR:
            return compute_fticr_field(ion);
        
        case Inst::Orbitrap:
            return compute_orbitrap_field(ion);
        
        case Inst::QuadrupoleRF:
            return compute_quadrupole_rf_field(ion, t);
        
        case Inst::NoFixedInstrument:
        case Inst::UnknownInstrument:
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
// 
// DC endcap confinement (axial):
//   - Harmonic potential: Φ(z) = α·z² where α = U_DC / L²
//   - Electric field: E_z = -2αz (z measured from trap center)
//   - Creates parabolic well for axial confinement
// ----------------------------------------------------------------------------
Vec3 ElectricFieldForce::compute_lqit_field(const IonState& ion, double t) const {
    Vec3 E_total{0.0, 0.0, 0.0};
    
    // SSOT: Read from config
    const auto& rf = domain_->fields.rf;
    const auto& dc = domain_->fields.dc;
    const auto& ac = domain_->fields.ac;
    const auto& geom = domain_->geometry;
    const auto& lib = domain_->fields.waveform_library;
    
    // v1.0: Evaluate voltages/frequencies at current time (only if set)
    const double rf_voltage = (rf.voltage_V.constant_value.has_value() || rf.voltage_V.waveform_ref.has_value()) 
                             ? eval_value(rf.voltage_V, t, lib) : 0.0;
    const double rf_freq = (rf.frequency_Hz.constant_value.has_value() || rf.frequency_Hz.waveform_ref.has_value()) 
                          ? eval_value(rf.frequency_Hz, t, lib) : 0.0;
    const double ac_voltage = (ac.voltage_V.constant_value.has_value() || ac.voltage_V.waveform_ref.has_value()) 
                             ? eval_value(ac.voltage_V, t, lib) : 0.0;
    const double ac_freq = (ac.frequency_Hz.constant_value.has_value() || ac.frequency_Hz.waveform_ref.has_value()) 
                          ? eval_value(ac.frequency_Hz, t, lib) : 0.0;
    const double dc_quad = (dc.quad_V.constant_value.has_value() || dc.quad_V.waveform_ref.has_value()) 
                          ? eval_value(dc.quad_V, t, lib) : 0.0;
    const double dc_axial = (dc.axial_V.constant_value.has_value() || dc.axial_V.waveform_ref.has_value()) 
                           ? eval_value(dc.axial_V, t, lib) : 0.0;
    
    // (1) RF + DC quadrupole field (radial) - matches RFField()
    if (std::fabs(rf_voltage) > MIN_VOLTAGE_THRESHOLD) {
        const double r0_sq = geom.radius_m * geom.radius_m;
        const double omega = 2.0 * M_PI * rf_freq;  // Compute omega from evaluated frequency
        const double U_eff = dc_quad + rf_voltage * std::cos(omega * t);
        
        E_total.x =  2.0 * ion.pos.x * U_eff / r0_sq;
        E_total.y = -2.0 * ion.pos.y * U_eff / r0_sq;
    }
    
    // (2) AC excitation field - matches ACField()
    if (std::fabs(ac_voltage) > MIN_VOLTAGE_THRESHOLD) {
        const double omega_ac = 2.0 * M_PI * ac_freq;  // Compute omega from evaluated frequency
        const double mag = -ac_voltage / geom.radius_m * std::cos(omega_ac * t);
        E_total.x += mag;  // Fixed x-direction
    }
    
    // (3) DC endcap field (axial confinement) - harmonic potential approximation
    // E_z = -2αz (z measured from trap center), where α = U_DC / L²
    // This creates a parabolic potential well: Φ(z) = α·z²
    if (std::fabs(dc_axial) > MIN_VOLTAGE_THRESHOLD) {
        const double L = geom.length_m;
        const double alpha = dc_axial / (L * L);
        
        // z measured from trap center
        const double z_centered = ion.pos.z - 0.5 * L;
        E_total.z = -2.0 * alpha * z_centered;
    }
    
    return E_total;
}

// ----------------------------------------------------------------------------
// IMS: Ion Mobility Spectrometry
// ----------------------------------------------------------------------------
// Uniform drift field: E_z = U_drift / L
// Optional radial RF field for focusing.
// Matches computeAccelerations.cpp implementation.
// ----------------------------------------------------------------------------
Vec3 ElectricFieldForce::compute_ims_field(const IonState& ion, double t) const {
    // SSOT: Read from config
    Vec3 E_total{0.0, 0.0, 0.0};
    const auto& dc = domain_->fields.dc;
    const auto& rf = domain_->fields.rf;
    const auto& lib = domain_->fields.waveform_library;
    
    // v1.0: Evaluate voltages at current time (only if set)
    const double dc_axial = (dc.axial_V.constant_value.has_value() || dc.axial_V.waveform_ref.has_value()) 
                           ? eval_value(dc.axial_V, t, lib) : 0.0;
    const double dc_quad = (dc.quad_V.constant_value.has_value() || dc.quad_V.waveform_ref.has_value()) 
                          ? eval_value(dc.quad_V, t, lib) : 0.0;
    const double rf_voltage = (rf.voltage_V.constant_value.has_value() || rf.voltage_V.waveform_ref.has_value()) 
                             ? eval_value(rf.voltage_V, t, lib) : 0.0;
    const double rf_freq = (rf.frequency_Hz.constant_value.has_value() || rf.frequency_Hz.waveform_ref.has_value()) 
                          ? eval_value(rf.frequency_Hz, t, lib) : 0.0;
    
    // (1) Axial DC field
    E_total.z = dc_axial / domain_->geometry.length_m;
    
    // (2) Optional radial RF field (for ion focusing)
    if (std::fabs(rf_voltage) > MIN_VOLTAGE_THRESHOLD) {
        const double r0_sq = domain_->geometry.radius_m * domain_->geometry.radius_m;
        const double omega = 2.0 * M_PI * rf_freq;
        const double U_eff = dc_quad + rf_voltage * std::cos(omega * t);
        
        E_total.x =  2.0 * ion.pos.x * U_eff / r0_sq;
        E_total.y = -2.0 * ion.pos.y * U_eff / r0_sq;
    }
    
    return E_total;
}

// ----------------------------------------------------------------------------
// TOF: Time-of-Flight
// ----------------------------------------------------------------------------
// DC field only in acceleration region, field-free drift region.
// Matches computeAccelerations.cpp behavior.
// ----------------------------------------------------------------------------
Vec3 ElectricFieldForce::compute_tof_field(const IonState& ion) const {
    // SSOT: Position-dependent field (acceleration vs drift region)
    const double acc_length = domain_->geometry.acc_length_m;
    
    // Field only in acceleration region (z < acc_length_m)
    if (ion.pos.z < acc_length) {
        // v1.0: Evaluate axial voltage at t=0 (DC field)
        const double dc_axial = eval_value(domain_->fields.dc.axial_V, 0.0, domain_->fields.waveform_library);
        return Vec3{
            0.0,
            0.0,
            dc_axial / acc_length
        };
    } else {
        // Field-free drift region
        return Vec3{0.0, 0.0, 0.0};
    }
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
    // SSOT: FTICR uses radial_V for quadrupolar trapping
    // Quadrupolar trapping: φ(r,z) = V_trap/d² · (2z² - x² - y²)
    // Matches FTICRField() in defineFields.cpp
    // v1.0: Evaluate radial voltage at t=0 (DC field)
    const double voltage = eval_value(domain_->fields.dc.radial_V, 0.0, domain_->fields.waveform_library);
    
    // Characteristic distance d = sqrt(L²/8 + r²/4) - matches computeAccelerations
    const double L = domain_->geometry.length_m;
    const double r = domain_->geometry.radius_m;
    const double d = std::sqrt((L * L / 8.0) + (r * r / 4.0));
    const double factor = voltage / (d * d);
    
    // Axial center (for z-dependence)
    const double z_center = ion.pos.z - 0.5 * L;
    
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
    // Transform position from global to local coordinates
    const Vec3 pos_local = domain_->rotation_global_to_local * (ion.pos - domain_->geometry.origin_m);
    
    // SSOT: Read from config
    // Hyperlogarithmic field: U(r,z) = k/2 · (z² - r²/2 + r_char²·ln(r/r_char))
    const double r_sq = std::max(1e-18, pos_local.x * pos_local.x + pos_local.y * pos_local.y);
    const double r_char = domain_->geometry.radius_char_m;  // Characteristic radius
    const double r_char_sq = r_char * r_char;
    const double C = 1.0 - r_char_sq / r_sq;
    
    // SSOT: Correct k calculation from computeAccelerations.cpp
    // v1.0: Evaluate radial voltage at t=0 (DC field)
    const double voltage = eval_value(domain_->fields.dc.radial_V, 0.0, domain_->fields.waveform_library);  // Radial voltage
    const double r_in = domain_->geometry.radius_in_m;
    const double r_out = domain_->geometry.radius_out_m;
    const double k = 2.0 * voltage / (r_char_sq * std::log(r_out / r_in)
                                     - 0.5 * (r_out * r_out - r_in * r_in));
    
    // Compute field in local coordinates
    // Orbitrap field uses z directly (not centered) - hyperboloid equation is z² - r²/2 = C
    // The field is naturally symmetric around z=0 (domain origin)
    Vec3 E_local{
        0.5 * k * pos_local.x * C,
        0.5 * k * pos_local.y * C,
       -k * pos_local.z
    };
    
    // Transform field back to global coordinates
    return domain_->rotation_local_to_global * E_local;
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

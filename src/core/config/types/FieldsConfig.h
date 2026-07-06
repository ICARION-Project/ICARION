// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2026 Christoph Schaefer
#ifndef ICARION_CONFIG_FIELDS_CONFIG_H
#define ICARION_CONFIG_FIELDS_CONFIG_H

#include "core/utils/mathUtils.h"
#include "../validation/ValidationResult.h"
#include "WaveformConfig.h"
#include <vector>
#include <string>
#include <map>
#include <cmath>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief DC field configuration
 * 
 * Supports both voltage specification and field strength (Townsend).
 * Voltages support time-varying waveforms.
 */
struct DCFieldConfig {
    // === Direct voltage specification (v1.0.0: static or waveform) ===
    ValueOrWaveform axial_V{0.0};       ///< Axial DC voltage [V]
    ValueOrWaveform quad_V{0.0};        ///< Quadrupole DC voltage [V]
    ValueOrWaveform radial_V{0.0};      ///< Radial DC voltage [V]
    
    // === Field strength specification (alternative to voltage) ===
    ValueOrWaveform EN_Td{0.0};         ///< Reduced field strength [Td]
    double EN_Vm2 = 0.0;                ///< E/N [V·m²] (computed from EN_Td)
    
    /**
     * @brief Validate DC field configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        // No hard constraints - DC fields can be zero or negative
        // Field strength and voltage are alternatives (both can be specified)
        return result;
    }
};

/**
 * @brief RF field configuration
 * 
 * Voltage and frequency support time-varying waveforms (chirps, modulation).
 */
struct RFFieldConfig {
    ValueOrWaveform voltage_V{0.0};     ///< RF amplitude [V] (0-to-peak, static or waveform)
    ValueOrWaveform frequency_Hz{0.0};  ///< RF frequency [Hz] (static or waveform)
    double phase_rad = 0.0;             ///< Initial phase [rad]
    
    // Derived (computed after load, only valid for static frequency)
    double angular_frequency_rad_s = 0.0;  ///< ω = 2π·f [rad/s] (static only)
    
    /**
     * @brief Compute derived quantities (for static frequency only)
     */
    void compute_derived() {
        // Only compute if frequency is static (not a waveform)
        if (frequency_Hz.constant_value.has_value()) {
            angular_frequency_rad_s = 2.0 * M_PI * frequency_Hz.constant_value.value();
        }
    }
    
    /**
     * @brief Validate RF field configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        // Validation deferred to runtime for waveforms
        // Static values checked here
        if (voltage_V.constant_value.has_value() && voltage_V.constant_value.value() < 0.0) {
            result.add_error("RF voltage cannot be negative");
        }
        if (frequency_Hz.constant_value.has_value() && frequency_Hz.constant_value.value() < 0.0) {
            result.add_error("RF frequency cannot be negative");
        }
        if (voltage_V.constant_value.has_value() && frequency_Hz.constant_value.has_value()) {
            if (voltage_V.constant_value.value() > 0.0 && frequency_Hz.constant_value.value() == 0.0) {
                result.add_error("RF voltage specified but frequency is zero");
            }
        }
        
        return result;
    }
};

/**
 * @brief One radial dipolar AC excitation axis.
 *
 * The instantaneous drive is amplitude_V(t) * ramp(t) *
 * cos(2*pi*frequency_Hz(t)*t + phase_rad).
 */
struct ACDipolarAxisConfig {
    bool enabled = false;
    ValueOrWaveform amplitude_V{0.0};   ///< AC amplitude [V] (static or waveform)
    ValueOrWaveform frequency_Hz{0.0};  ///< AC frequency [Hz] (static or waveform)
    double phase_rad = 0.0;             ///< Initial phase [rad]
    ValueOrWaveform ramp{1.0};          ///< Dimensionless amplitude multiplier

    double angular_frequency_rad_s = 0.0;  ///< Derived only for static frequency

    void compute_derived() {
        if (frequency_Hz.constant_value.has_value()) {
            angular_frequency_rad_s = 2.0 * M_PI * frequency_Hz.constant_value.value();
        }
    }

    ValidationResult validate(const std::string& axis_name) const {
        ValidationResult result;

        if (!amplitude_V.is_valid()) {
            result.add_error("AC dipolar axis '" + axis_name + "' amplitude_V must be a valid ValueOrWaveform");
        }
        if (!frequency_Hz.is_valid()) {
            result.add_error("AC dipolar axis '" + axis_name + "' frequency_Hz must be a valid ValueOrWaveform");
        }
        if (!ramp.is_valid()) {
            result.add_error("AC dipolar axis '" + axis_name + "' ramp must be a valid ValueOrWaveform");
        }
        if (amplitude_V.constant_value.has_value() && amplitude_V.constant_value.value() < 0.0) {
            result.add_error("AC dipolar axis '" + axis_name + "' amplitude_V cannot be negative");
        }
        if (frequency_Hz.constant_value.has_value() && frequency_Hz.constant_value.value() < 0.0) {
            result.add_error("AC dipolar axis '" + axis_name + "' frequency_Hz cannot be negative");
        }
        if (enabled &&
            amplitude_V.constant_value.has_value() &&
            frequency_Hz.constant_value.has_value() &&
            amplitude_V.constant_value.value() > 0.0 &&
            frequency_Hz.constant_value.value() == 0.0) {
            result.add_error("AC dipolar axis '" + axis_name + "' has nonzero amplitude but zero frequency");
        }

        return result;
    }
};

/**
 * @brief AC excitation field configuration (primarily for LQIT)
 *
 * Legacy voltage/frequency fields drive the local x-axis. The optional
 * dipolar axes allow independent x/y amplitudes, frequencies, phases, and ramps.
 */
struct ACFieldConfig {
    ValueOrWaveform voltage_V{0.0};     ///< AC amplitude [V] (static or waveform)
    ValueOrWaveform frequency_Hz{0.0};  ///< AC frequency [Hz] (static or waveform)
    double phase_rad = 0.0;             ///< Legacy single-axis phase [rad]
    
    // Derived (only valid for static frequency)
    double angular_frequency_rad_s = 0.0;  ///< ω = 2π·f [rad/s] (static only)

    // === LQIT two-axis dipolar excitation ===
    bool dipolar_excitation_defined = false;
    ACDipolarAxisConfig dipolar_x;
    ACDipolarAxisConfig dipolar_y;
    
    // === LQIT phase locking to RF ===
    bool lqit_lock_enable = false;
    double lqit_lock_phase_rad = 0.0;   ///< Phase offset to RF [rad]
    double lqit_lock_bandwidth_Hz = 0.0; ///< Lock bandwidth [Hz]
    
    /**
     * @brief Compute derived quantities (for static frequency only)
     */
    void compute_derived() {
        if (frequency_Hz.constant_value.has_value()) {
            angular_frequency_rad_s = 2.0 * M_PI * frequency_Hz.constant_value.value();
        }
        dipolar_x.compute_derived();
        dipolar_y.compute_derived();
    }
    
    /**
     * @brief Validate AC field configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        // Static value validation
        if (voltage_V.constant_value.has_value() && voltage_V.constant_value.value() < 0.0) {
            result.add_error("AC voltage cannot be negative");
        }
        if (frequency_Hz.constant_value.has_value() && frequency_Hz.constant_value.value() < 0.0) {
            result.add_error("AC frequency cannot be negative");
        }
        if (voltage_V.constant_value.has_value() && frequency_Hz.constant_value.has_value()) {
            if (voltage_V.constant_value.value() > 0.0 && frequency_Hz.constant_value.value() == 0.0) {
                result.add_error("AC voltage specified but frequency is zero");
            }
        }
        if (dipolar_excitation_defined) {
            result.merge(dipolar_x.validate("x"));
            result.merge(dipolar_y.validate("y"));
        }
        
        return result;
    }
};

/**
 * @brief Magnetic field configuration
 * 
 * Supports uniform field + optional linear gradient.
 */
struct MagneticFieldConfig {
    bool enabled = false;
    Vec3 field_strength_T = {0.0, 0.0, 0.0};    ///< Uniform B field [T]
    Vec3 field_gradient_T_m = {0.0, 0.0, 0.0};  ///< Linear gradient dB/dz [T/m]
    
    /**
     * @brief Validate magnetic field configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        // No hard constraints - B can be arbitrary
        if (enabled) {
            double B_mag = std::sqrt(
                field_strength_T.x * field_strength_T.x +
                field_strength_T.y * field_strength_T.y +
                field_strength_T.z * field_strength_T.z
            );
            if (B_mag > 100.0) {
                result.add_warning("Magnetic field > 100 T seems unrealistic");
            }
        }
        
        return result;
    }
};

/**
 * @brief TIMS axial field program.
 *
 * Evaluates E_z(z,t) = (1 - f(t)) * E_initial(z) + f(t) * E_final(z),
 * with optional 1D profiles along local z. RF/DC radial confinement still
 * comes from standard IMS RF/DC fields.
 */
struct TIMSFieldConfig {
    bool enabled = false;

    double axial_field_initial_uniform_V_m = 0.0; ///< Initial axial field before the ramp [V/m]
    double axial_field_final_uniform_V_m = 0.0;   ///< Final axial field after the ramp [V/m]

    std::vector<double> z_positions_m;          ///< Local z grid for profiles [m]
    std::vector<double> axial_field_initial_profile_V_m; ///< Initial axial field profile [V/m]
    std::vector<double> axial_field_final_profile_V_m;   ///< Final axial field profile [V/m]

    double ramp_start_s = 0.0;
    double ramp_end_s = 0.0;
    std::string ramp_mode = "linear";  ///< "linear" or "exponential"
    double ramp_tau_s = 1e-3;          ///< Exponential time constant [s]

    bool ramp_fraction_defined = false;
    ValueOrWaveform ramp_fraction{0.0};

    ValidationResult validate() const {
        ValidationResult result;
        if (!enabled) {
            return result;
        }

        if (ramp_end_s < ramp_start_s) {
            result.add_error("TIMS ramp_end_s must be >= ramp_start_s");
        }
        if (ramp_mode != "linear" && ramp_mode != "exponential") {
            result.add_error("TIMS ramp_mode must be 'linear' or 'exponential'");
        }
        if (ramp_mode == "exponential" && ramp_tau_s <= 0.0) {
            result.add_error("TIMS ramp_tau_s must be > 0 for exponential mode");
        }
        if (ramp_fraction_defined && !ramp_fraction.is_valid()) {
            result.add_error("TIMS ramp_fraction must be a valid ValueOrWaveform");
        }

        if (!z_positions_m.empty()) {
            if (z_positions_m.size() < 2) {
                result.add_error("TIMS z_positions_m must contain at least two points");
            }
            if (axial_field_initial_profile_V_m.size() != z_positions_m.size()) {
                result.add_error("TIMS axial_field_initial_profile_V_m size must match z_positions_m");
            }
            if (!axial_field_final_profile_V_m.empty() && axial_field_final_profile_V_m.size() != z_positions_m.size()) {
                result.add_error("TIMS axial_field_final_profile_V_m size must match z_positions_m when provided");
            }
            for (size_t i = 1; i < z_positions_m.size(); ++i) {
                if (z_positions_m[i] < z_positions_m[i - 1]) {
                    result.add_error("TIMS z_positions_m must be non-decreasing");
                    break;
                }
            }
        }

        return result;
    }
};

/**
 * @brief Complete field configuration for a domain
 * 
 * Aggregates DC, RF, AC, and magnetic fields.
 * Also includes precomputed field arrays (BEM/FEM results).
 * v1.0.0: Includes waveform library for named waveforms.
 */
struct FieldsConfig {
    DCFieldConfig dc;
    RFFieldConfig rf;
    ACFieldConfig ac;
    TIMSFieldConfig tims;
    MagneticFieldConfig magnetic;
    
    // === Waveform library (v1.0.0) ===
    std::map<std::string, Waveform> waveform_library;  ///< Named waveforms for @references
    
    // === Precomputed field arrays (BEM/FEM) ===
    
    /**
     * @brief Field array term for superposition
     * 
     * Each term represents a unit-voltage field from BEM/FEM,
     * scaled by a time-dependent factor.
     */
    struct FieldArrayTerm {
        std::string file;               ///< Path to HDF5 field file
        
        enum class ScaleKind {
            Constant,                   ///< Scale by constant factor
            DC_Axial,                   ///< Scale by DC.axial_V
            DC_Quad,                    ///< Scale by DC.quad_V
            DC_Radial,                  ///< Scale by DC.radial_V
            RF                          ///< Scale by RF.voltage_V * cos(ωt + φ)
        } kind = ScaleKind::Constant;
        
        double constant = 1.0;          ///< Constant scaling factor [V]
        double phase_rad = 0.0;         ///< Phase offset for RF terms [rad]
        double frequency_Hz = 0.0;      ///< Frequency override for RF (0 = use domain RF freq)
    };
    
    std::vector<FieldArrayTerm> field_array_terms;  ///< Superposition of field arrays
    
    // Legacy single-file support (deprecated, use field_array_terms)
    std::string legacy_field_array_file = "";
    
    /**
     * @brief Validate complete fields configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        result.merge(dc.validate());
        result.merge(rf.validate());
        result.merge(ac.validate());
        result.merge(tims.validate());
        result.merge(magnetic.validate());
        
        // Field array validation
        for (const auto& term : field_array_terms) {
            if (term.file.empty()) {
                result.add_error("Field array term has empty file path");
            }
            if (term.kind == FieldArrayTerm::ScaleKind::RF && term.frequency_Hz < 0.0) {
                result.add_error("Field array RF term has negative frequency");
            }
        }
        
        return result;
    }
    
    /**
     * @brief Compute derived quantities for all sub-configs
     */
    void compute_derived() {
        rf.compute_derived();
        ac.compute_derived();
    }
};

} // namespace ICARION::config

#endif // ICARION_CONFIG_FIELDS_CONFIG_H

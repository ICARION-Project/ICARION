// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 ICARION Project Contributors
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
 * @brief AC excitation field configuration (primarily for LQIT)
 * 
 * Voltage and frequency support time-varying waveforms.
 */
struct ACFieldConfig {
    ValueOrWaveform voltage_V{0.0};     ///< AC amplitude [V] (static or waveform)
    ValueOrWaveform frequency_Hz{0.0};  ///< AC frequency [Hz] (static or waveform)
    
    // Derived (only valid for static frequency)
    double angular_frequency_rad_s = 0.0;  ///< ω = 2π·f [rad/s] (static only)
    
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

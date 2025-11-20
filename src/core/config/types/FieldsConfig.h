// SPDX-License-Identifier: Apache-2.0
// SPDX-FileCopyrightText: 2025 ICARION Project Contributors

#ifndef ICARION_CONFIG_FIELDS_CONFIG_H
#define ICARION_CONFIG_FIELDS_CONFIG_H

#include "core/utils/mathUtils.h"
#include "../validation/ValidationResult.h"
#include <vector>
#include <string>
#include <cmath>
#include <stdexcept>

namespace ICARION::config {

/**
 * @brief DC field configuration
 * 
 * Supports both voltage specification and field strength (Townsend).
 * Future: Time-varying voltages (ramping, pulses).
 */
struct DCFieldConfig {
    // === Direct voltage specification ===
    double axial_V = 0.0;               ///< Axial DC voltage [V]
    double quad_V = 0.0;                ///< Quadrupole DC voltage [V]
    double radial_V = 0.0;              ///< Radial DC voltage [V]
    
    // === Field strength specification (alternative to voltage) ===
    double EN_Td = 0.0;                 ///< Reduced field strength [Td]
    double EN_Vm2 = 0.0;                ///< E/N [V·m²]
    
    // Future: Time-varying support (voltage ramping)
    // bool enable_ramp = false;
    // double ramp_start_s = 0.0;
    // double ramp_rate_V_s = 0.0;
    // std::vector<std::pair<double, double>> voltage_schedule;  // {time, voltage}
    
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
 * Future: Frequency chirps, amplitude modulation.
 */
struct RFFieldConfig {
    double voltage_V = 0.0;             ///< RF amplitude [V] (0-to-peak)
    double frequency_Hz = 0.0;          ///< RF frequency [Hz]
    double phase_rad = 0.0;             ///< Initial phase [rad]
    
    // Derived (computed after load)
    double angular_frequency_rad_s = 0.0;  ///< ω = 2π·f [rad/s]
    
    /**
     * @brief Compute derived quantities
     */
    void compute_derived() {
        angular_frequency_rad_s = 2.0 * M_PI * frequency_Hz;
    }
    
    /**
     * @brief Validate RF field configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        if (voltage_V < 0.0) {
            result.add_error("RF voltage cannot be negative");
        }
        if (frequency_Hz < 0.0) {
            result.add_error("RF frequency cannot be negative");
        }
        if (voltage_V > 0.0 && frequency_Hz == 0.0) {
            result.add_error("RF voltage specified but frequency is zero");
        }
        
        return result;
    }
};

/**
 * @brief AC excitation field configuration (primarily for LQIT)
 * 
 * Supports voltage and frequency sweeps.
 * Future: Replace sweeps with general waveform system.
 */
struct ACFieldConfig {
    double voltage_V = 0.0;             ///< AC amplitude [V]
    double frequency_Hz = 0.0;          ///< AC frequency [Hz]
    
    // Derived
    double angular_frequency_rad_s = 0.0;  ///< ω = 2π·f [rad/s]
    
    // === Voltage sweep (linear) ===
    bool enable_voltage_sweep = false;
    double amplitude_slope_V_s = 0.0;   ///< Voltage ramp rate [V/s]
    double start_time_s = 0.0;          ///< Sweep start time [s]
    double rise_time_s = 0.0;           ///< Sweep duration [s]
    
    // === Frequency sweep (linear) ===
    bool enable_frequency_sweep = false;
    double frequency_start_Hz = 0.0;    ///< Initial frequency [Hz]
    double frequency_sweep_slope_Hz_s = 0.0; ///< Frequency ramp rate [Hz/s]
    
    // === LQIT phase locking to RF ===
    bool lqit_lock_enable = false;
    double lqit_lock_phase_rad = 0.0;   ///< Phase offset to RF [rad]
    double lqit_lock_bandwidth_Hz = 0.0; ///< Lock bandwidth [Hz]
    
    // === Arbitrary waveforms (future) ===
    std::vector<std::pair<double, double>> voltage_time_table; ///< {time_s, voltage_V}
    
    /**
     * @brief Compute derived quantities
     */
    void compute_derived() {
        angular_frequency_rad_s = 2.0 * M_PI * frequency_Hz;
    }
    
    /**
     * @brief Validate AC field configuration
     */
    ValidationResult validate() const {
        ValidationResult result;
        
        if (voltage_V < 0.0) {
            result.add_error("AC voltage cannot be negative");
        }
        if (frequency_Hz < 0.0) {
            result.add_error("AC frequency cannot be negative");
        }
        
        // Sweep validation
        if (enable_voltage_sweep) {
            if (rise_time_s <= 0.0) {
                result.add_error("AC voltage sweep rise_time_s must be positive");
            }
            if (start_time_s < 0.0) {
                result.add_error("AC voltage sweep start_time_s cannot be negative");
            }
        }
        
        if (enable_frequency_sweep) {
            if (frequency_start_Hz < 0.0) {
                result.add_error("AC frequency sweep start frequency cannot be negative");
            }
        }
        
        // Voltage time table validation
        if (!voltage_time_table.empty()) {
            for (size_t i = 1; i < voltage_time_table.size(); ++i) {
                if (voltage_time_table[i].first <= voltage_time_table[i-1].first) {
                    result.add_error("AC voltage_time_table must be sorted by time");
                }
            }
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
 */
struct FieldsConfig {
    DCFieldConfig dc;
    RFFieldConfig rf;
    ACFieldConfig ac;
    MagneticFieldConfig magnetic;
    
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

// ICARION: Ion Collision And Reaction IntegratiON
// MIT License - Copyright (c) 2025 ICARION Project Contributors

#ifndef ICARION_WAVEFORM_CONFIG_H
#define ICARION_WAVEFORM_CONFIG_H

#include <variant>
#include <vector>
#include <string>
#include <optional>
#include <map>
#include <cmath>

namespace ICARION::config {

/**
 * @brief Waveform types for time-varying parameters
 */
enum class WaveformType {
    Constant,       ///< y = value
    Linear,         ///< y = start + (end - start) * t_norm
    Quadratic,      ///< y = a + b*t + c*t²
    Sinusoidal,     ///< y = offset + amplitude * sin(2π*f*t + φ)
    Pulsed,         ///< y = high (during pulse), else low
    Arbitrary       ///< Interpolated from time table
};

/**
 * @brief Constant waveform (static value)
 */
struct ConstantWaveform {
    double value = 0.0;
    
    double evaluate([[maybe_unused]] double t) const { return value; }
};

/**
 * @brief Linear ramp/sweep waveform
 */
struct LinearWaveform {
    double start_value = 0.0;
    double end_value = 0.0;
    double start_time_s = 0.0;
    double end_time_s = 1e9;
    bool clamp = true;  ///< Hold end_value after end_time
    
    double evaluate(double t) const {
        if (t < start_time_s) return start_value;
        if (t >= end_time_s) {
            return clamp ? end_value : start_value;
        }
        double duration = end_time_s - start_time_s;
        double t_norm = (t - start_time_s) / duration;
        return start_value + (end_value - start_value) * t_norm;
    }
};

/**
 * @brief Quadratic ramp (acceleration profile)
 * y = a + b·t + c·t²  (within [start_time_s, end_time_s])
 */
struct QuadraticWaveform {
    double a = 0.0;  ///< Constant term
    double b = 0.0;  ///< Linear term
    double c = 0.0;  ///< Quadratic term
    double start_time_s = 0.0;
    double end_time_s = 1e9;
    
    double evaluate(double t) const {
        if (t < start_time_s || t > end_time_s) return a;
        return a + b * t + c * t * t;
    }
};

/**
 * @brief Sinusoidal modulation waveform
 */
struct SinusoidalWaveform {
    double offset = 0.0;
    double amplitude = 0.0;
    double frequency_Hz = 0.0;
    double phase_rad = 0.0;
    
    double evaluate(double t) const {
        double omega = 2.0 * M_PI * frequency_Hz;
        return offset + amplitude * std::sin(omega * t + phase_rad);
    }
};

/**
 * @brief Pulsed waveform (single pulse)
 */
struct PulsedWaveform {
    double low_value = 0.0;
    double high_value = 0.0;
    double pulse_start_s = 0.0;
    double pulse_width_s = 0.0;
    
    double evaluate(double t) const {
        if (t >= pulse_start_s && t < pulse_start_s + pulse_width_s) {
            return high_value;
        }
        return low_value;
    }
};

/**
 * @brief Arbitrary waveform (interpolated from time table)
 */
struct ArbitraryWaveform {
    std::vector<double> times_s;
    std::vector<double> values;
    
    enum class Interpolation {
        Linear,     ///< Linear interpolation
        Step,       ///< Step function (hold previous)
        Cubic       ///< Cubic spline (smooth)
    };
    Interpolation interp = Interpolation::Linear;
    
    double evaluate(double t) const;
};

/**
 * @brief Waveform variant (type-safe union of all waveform types)
 */
using WaveformVariant = std::variant<
    ConstantWaveform,
    LinearWaveform,
    QuadraticWaveform,
    SinusoidalWaveform,
    PulsedWaveform,
    ArbitraryWaveform
>;

/**
 * @brief Generic waveform container
 */
struct Waveform {
    std::string id;  ///< Optional ID for named waveforms
    WaveformVariant data;
    
    /**
     * @brief Evaluate waveform at time t
     * @param t Time in seconds
     * @return Waveform value at time t
     */
    double evaluate(double t) const {
        return std::visit([t](const auto& w) { return w.evaluate(t); }, data);
    }
};

/**
 * @brief Flexible field value: static OR time-varying
 * 
 * Exactly one option must be set:
 * - constant_value: Static value (most common)
 * - waveform: Inline waveform definition
 * - waveform_ref: Reference to named waveform
 */
struct ValueOrWaveform {
    std::optional<double> constant_value;
    std::optional<Waveform> waveform;
    std::optional<std::string> waveform_ref;
    
    // Default constructor (all nullopt)
    ValueOrWaveform() = default;
    
    // Convenience constructor for constant value
    explicit ValueOrWaveform(double value) : constant_value(value) {}
    
    /**
     * @brief Check if exactly one option is set
     */
    bool is_valid() const {
        int count = (constant_value.has_value() ? 1 : 0) +
                    (waveform.has_value() ? 1 : 0) +
                    (waveform_ref.has_value() ? 1 : 0);
        return count == 1;
    }
    
    /**
     * @brief Check if time-varying (not constant)
     */
    bool is_time_varying() const {
        return waveform.has_value() || waveform_ref.has_value();
    }
    
    /**
     * @brief Evaluate at time t
     * @param t Time in seconds
     * @param library Waveform library (for resolving references)
     * @return Value at time t
     */
    double evaluate(double t, const std::map<std::string, Waveform>& library) const;
};

} // namespace ICARION::config

#endif // ICARION_WAVEFORM_CONFIG_H

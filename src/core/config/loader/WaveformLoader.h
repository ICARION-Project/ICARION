// ICARION: Ion Collision And Reaction IntegratiON
// SPDX-License-Identifier: GPL-3.0-only
// Copyright (c) 2025 Christoph Schaefer
#ifndef ICARION_WAVEFORM_LOADER_H
#define ICARION_WAVEFORM_LOADER_H

#include "core/config/types/WaveformConfig.h"
#include <json/json.h>
#include <map>
#include <string>

namespace ICARION::config {

/**
 * @brief JSON loader for waveform configurations
 */
class WaveformLoader {
public:
    /**
     * @brief Load a single waveform from JSON object
     * @param json JSON object with "type" field and type-specific parameters
     * @return Waveform object
     * @throws std::runtime_error if JSON is invalid or type is unknown
     */
    static Waveform load(const Json::Value& json);
    
    /**
     * @brief Load waveform library (map of named waveforms)
     * @param json JSON object where keys are waveform IDs
     * @return Map of waveform ID → Waveform
     * @throws std::runtime_error if JSON is invalid
     */
    static std::map<std::string, Waveform> load_library(const Json::Value& json);
    
    /**
     * @brief Load ValueOrWaveform (handles static value, inline, or reference)
     * @param json JSON value (number, object, or string)
     * @param local_library Domain-local waveform library
     * @param global_library Global waveform library (SSOT)
     * @return ValueOrWaveform object
     * @throws std::runtime_error if JSON format is invalid
     */
    static ValueOrWaveform load_value_or_waveform(
        const Json::Value& json,
        const std::map<std::string, Waveform>& local_library,
        const std::map<std::string, Waveform>& global_library = {}
    );
    
private:
    /**
     * @brief Load specific waveform types
     */
    static ConstantWaveform load_constant(const Json::Value& json);
    static LinearWaveform load_linear(const Json::Value& json);
    static QuadraticWaveform load_quadratic(const Json::Value& json);
    static ExponentialWaveform load_exponential(const Json::Value& json);
    static SinusoidalWaveform load_sinusoidal(const Json::Value& json);
    static PWMWaveform load_pwm(const Json::Value& json);
    static PulsedWaveform load_pulsed(const Json::Value& json);
    static ArbitraryWaveform load_arbitrary(const Json::Value& json);
    
    /**
     * @brief Helper: Get required field with type checking
     */
    template<typename T>
    static T get_required(const Json::Value& json, const std::string& field, const std::string& type_name);
    
    /**
     * @brief Helper: Get optional field with default value
     */
    template<typename T>
    static T get_optional(const Json::Value& json, const std::string& field, T default_value);
};

} // namespace ICARION::config

#endif // ICARION_WAVEFORM_LOADER_H

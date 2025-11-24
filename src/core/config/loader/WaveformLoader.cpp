#include "WaveformLoader.h"
#include <stdexcept>
#include <sstream>

namespace ICARION::config {

// Helper to get required field
template<typename T>
T WaveformLoader::get_required(const Json::Value& json, const std::string& field, const std::string& type_name) {
    if (!json.isMember(field)) {
        throw std::runtime_error("Missing required field '" + field + "' in " + type_name + " waveform");
    }
    
    const Json::Value& val = json[field];
    
    if constexpr (std::is_same_v<T, double>) {
        if (!val.isNumeric()) {
            throw std::runtime_error("Field '" + field + "' in " + type_name + " must be numeric");
        }
        return val.asDouble();
    } else if constexpr (std::is_same_v<T, bool>) {
        if (!val.isBool()) {
            throw std::runtime_error("Field '" + field + "' in " + type_name + " must be boolean");
        }
        return val.asBool();
    } else if constexpr (std::is_same_v<T, std::string>) {
        if (!val.isString()) {
            throw std::runtime_error("Field '" + field + "' in " + type_name + " must be string");
        }
        return val.asString();
    }
    
    throw std::runtime_error("Unsupported type in get_required");
}

// Helper to get optional field
template<typename T>
T WaveformLoader::get_optional(const Json::Value& json, const std::string& field, T default_value) {
    if (!json.isMember(field)) {
        return default_value;
    }
    
    const Json::Value& val = json[field];
    
    if constexpr (std::is_same_v<T, double>) {
        if (!val.isNumeric()) return default_value;
        return val.asDouble();
    } else if constexpr (std::is_same_v<T, bool>) {
        if (!val.isBool()) return default_value;
        return val.asBool();
    } else if constexpr (std::is_same_v<T, std::string>) {
        if (!val.isString()) return default_value;
        return val.asString();
    }
    
    return default_value;
}

ConstantWaveform WaveformLoader::load_constant(const Json::Value& json) {
    ConstantWaveform w;
    w.value = get_required<double>(json, "value", "constant");
    return w;
}

LinearWaveform WaveformLoader::load_linear(const Json::Value& json) {
    LinearWaveform w;
    w.start_value = get_required<double>(json, "start", "linear");
    w.end_value = get_required<double>(json, "end", "linear");
    w.duration_s = get_required<double>(json, "duration_s", "linear");
    w.start_time_s = get_optional<double>(json, "start_time_s", 0.0);
    w.clamp = get_optional<bool>(json, "clamp", true);
    
    if (w.duration_s <= 0.0) {
        throw std::runtime_error("LinearWaveform: duration_s must be positive");
    }
    
    return w;
}

QuadraticWaveform WaveformLoader::load_quadratic(const Json::Value& json) {
    QuadraticWaveform w;
    w.a = get_required<double>(json, "a", "quadratic");
    w.b = get_required<double>(json, "b", "quadratic");
    w.c = get_required<double>(json, "c", "quadratic");
    w.start_time_s = get_optional<double>(json, "start_time_s", 0.0);
    w.end_time_s = get_optional<double>(json, "end_time_s", 1e9);
    
    return w;
}

SinusoidalWaveform WaveformLoader::load_sinusoidal(const Json::Value& json) {
    SinusoidalWaveform w;
    w.amplitude = get_required<double>(json, "amplitude", "sinusoidal");
    w.frequency_Hz = get_required<double>(json, "frequency_Hz", "sinusoidal");
    w.offset = get_optional<double>(json, "offset", 0.0);
    w.phase_rad = get_optional<double>(json, "phase_rad", 0.0);
    
    if (w.frequency_Hz <= 0.0) {
        throw std::runtime_error("SinusoidalWaveform: frequency_Hz must be positive");
    }
    
    return w;
}

PulsedWaveform WaveformLoader::load_pulsed(const Json::Value& json) {
    PulsedWaveform w;
    w.low_value = get_required<double>(json, "low", "pulsed");
    w.high_value = get_required<double>(json, "high", "pulsed");
    w.pulse_start_s = get_required<double>(json, "pulse_start_s", "pulsed");
    w.pulse_width_s = get_required<double>(json, "pulse_width_s", "pulsed");
    
    if (w.pulse_width_s <= 0.0) {
        throw std::runtime_error("PulsedWaveform: pulse_width_s must be positive");
    }
    if (w.pulse_start_s < 0.0) {
        throw std::runtime_error("PulsedWaveform: pulse_start_s must be non-negative");
    }
    
    return w;
}

ArbitraryWaveform WaveformLoader::load_arbitrary(const Json::Value& json) {
    ArbitraryWaveform w;
    
    // Load times array
    if (!json.isMember("times")) {
        throw std::runtime_error("ArbitraryWaveform: missing required field 'times'");
    }
    const Json::Value& times = json["times"];
    if (!times.isArray()) {
        throw std::runtime_error("ArbitraryWaveform: 'times' must be an array");
    }
    
    for (const auto& t : times) {
        if (!t.isNumeric()) {
            throw std::runtime_error("ArbitraryWaveform: all elements in 'times' must be numeric");
        }
        w.times_s.push_back(t.asDouble());
    }
    
    // Load values array
    if (!json.isMember("values")) {
        throw std::runtime_error("ArbitraryWaveform: missing required field 'values'");
    }
    const Json::Value& values = json["values"];
    if (!values.isArray()) {
        throw std::runtime_error("ArbitraryWaveform: 'values' must be an array");
    }
    
    for (const auto& v : values) {
        if (!v.isNumeric()) {
            throw std::runtime_error("ArbitraryWaveform: all elements in 'values' must be numeric");
        }
        w.values.push_back(v.asDouble());
    }
    
    // Validate array lengths
    if (w.times_s.size() != w.values.size()) {
        throw std::runtime_error("ArbitraryWaveform: 'times' and 'values' must have same length");
    }
    if (w.times_s.size() < 2) {
        throw std::runtime_error("ArbitraryWaveform: need at least 2 points");
    }
    
    // Validate times are sorted
    for (size_t i = 1; i < w.times_s.size(); ++i) {
        if (w.times_s[i] <= w.times_s[i-1]) {
            throw std::runtime_error("ArbitraryWaveform: times must be strictly increasing");
        }
    }
    
    // Load interpolation mode
    std::string interp_str = get_optional<std::string>(json, "interpolation", "linear");
    if (interp_str == "linear") {
        w.interp = ArbitraryWaveform::Interpolation::Linear;
    } else if (interp_str == "step") {
        w.interp = ArbitraryWaveform::Interpolation::Step;
    } else if (interp_str == "cubic") {
        w.interp = ArbitraryWaveform::Interpolation::Cubic;
    } else {
        throw std::runtime_error("ArbitraryWaveform: unknown interpolation mode '" + interp_str + "'");
    }
    
    return w;
}

Waveform WaveformLoader::load(const Json::Value& json) {
    if (!json.isObject()) {
        throw std::runtime_error("Waveform must be a JSON object");
    }
    
    if (!json.isMember("type")) {
        throw std::runtime_error("Waveform must have 'type' field");
    }
    
    std::string type = json["type"].asString();
    
    Waveform waveform;
    
    if (type == "constant") {
        waveform.data = load_constant(json);
    } else if (type == "linear") {
        waveform.data = load_linear(json);
    } else if (type == "quadratic") {
        waveform.data = load_quadratic(json);
    } else if (type == "sinusoidal") {
        waveform.data = load_sinusoidal(json);
    } else if (type == "pulsed") {
        waveform.data = load_pulsed(json);
    } else if (type == "arbitrary") {
        waveform.data = load_arbitrary(json);
    } else {
        throw std::runtime_error("Unknown waveform type: " + type);
    }
    
    // Optional: Load waveform ID
    if (json.isMember("id")) {
        waveform.id = json["id"].asString();
    }
    
    return waveform;
}

std::map<std::string, Waveform> WaveformLoader::load_library(const Json::Value& json) {
    std::map<std::string, Waveform> library;
    
    if (!json.isObject()) {
        throw std::runtime_error("Waveform library must be a JSON object");
    }
    
    for (const auto& id : json.getMemberNames()) {
        Waveform w = load(json[id]);
        w.id = id;  // Set ID from key
        library[id] = w;
    }
    
    return library;
}

ValueOrWaveform WaveformLoader::load_value_or_waveform(
    const Json::Value& json,
    const std::map<std::string, Waveform>& library
) {
    ValueOrWaveform result;
    
    if (json.isNumeric()) {
        // OLD FORMAT: Static value
        result.constant_value = json.asDouble();
    } else if (json.isString()) {
        // NEW FORMAT: Waveform reference (@waveform_id)
        std::string ref = json.asString();
        if (ref.empty() || ref[0] != '@') {
            throw std::runtime_error("Waveform reference must start with '@' (got: '" + ref + "')");
        }
        result.waveform_ref = ref.substr(1);  // Remove '@'
        
        // Validate reference exists
        if (library.find(*result.waveform_ref) == library.end()) {
            throw std::runtime_error("Waveform reference not found in library: " + *result.waveform_ref);
        }
    } else if (json.isObject()) {
        // NEW FORMAT: Inline waveform definition
        result.waveform = load(json);
    } else {
        throw std::runtime_error("ValueOrWaveform: invalid JSON type (must be number, string, or object)");
    }
    
    if (!result.is_valid()) {
        throw std::runtime_error("ValueOrWaveform: internal error (invalid state after loading)");
    }
    
    return result;
}

} // namespace ICARION::config

# Waveform System Design Document

**Author:** Architecture Team  
**Date:** November 24, 2025  
**Status:** PROPOSED (awaiting approval)  
**Target Release:** v1.1 (post Phase 6)  
**Estimated Implementation:** 11 hours

---

## Executive Summary

**Problem:** Current field configuration uses separate boolean flags + parameters for each sweep type, leading to:
- **Boolean soup:** 7+ fields per parameter for voltage sweeps
- **Not extensible:** Adding new sweep types requires new booleans
- **Inconsistent:** Different approaches for voltage (sweep) vs arbitrary (time_table)

**Solution:** Unified **Waveform System** that replaces all sweep flags with a single flexible type:
```cpp
// OLD (7 fields for voltage sweep):
double voltage_V;
bool enable_voltage_sweep;
double amplitude_slope_V_s;
double start_time_s;
double rise_time_s;

// NEW (1 field):
ValueOrWaveform voltage_V;  // Static OR time-varying
```

**Benefits:**
- ✅ **75% field reduction** (8 → 2 fields in ACFieldConfig)
- ✅ **Extensible:** Add new waveform types without API changes
- ✅ **Consistent:** Same pattern for all parameters
- ✅ **Reusable:** Named waveforms shareable across domains
- ✅ **Backward compatible:** Old JSON configs still work

---

## Table of Contents

1. [Motivation & Current Problems](#motivation--current-problems)
2. [Design Goals](#design-goals)
3. [Proposed Architecture](#proposed-architecture)
4. [JSON Schema & Examples](#json-schema--examples)
5. [Implementation Plan](#implementation-plan)
6. [Migration Strategy](#migration-strategy)
7. [Testing Strategy](#testing-strategy)
8. [Performance Considerations](#performance-considerations)
9. [Future Extensions](#future-extensions)
10. [Decision Log](#decision-log)

---

## Motivation & Current Problems

### Current Implementation (v1.0)

**File:** `src/core/config/types/FieldsConfig.h:95-120`

```cpp
struct ACFieldConfig {
    double voltage_V = 0.0;
    double frequency_Hz = 0.0;
    
    // === Voltage sweep (linear) ===
    bool enable_voltage_sweep = false;
    double amplitude_slope_V_s = 0.0;
    double start_time_s = 0.0;
    double rise_time_s = 0.0;
    
    // === Frequency sweep (linear) ===
    bool enable_frequency_sweep = false;
    double frequency_start_Hz = 0.0;
    double frequency_sweep_slope_Hz_s = 0.0;
    
    // === Arbitrary waveforms (inconsistent) ===
    std::vector<std::pair<double, double>> voltage_time_table;
};
```

### Problems Identified

#### Problem 1: Boolean Soup Anti-Pattern
```cpp
// Want to add sinusoidal modulation?
bool enable_sinusoidal_modulation;
double sinusoidal_amplitude;
double sinusoidal_frequency_Hz;
double sinusoidal_phase_rad;
// = 4 MORE fields per parameter!
```

**Impact:**
- ACFieldConfig: 12 fields (2 parameters × 6 fields each)
- DCFieldConfig potential: 18 fields (3 parameters × 6 fields)
- **Validation hell:** Check enable_X → check param_1 → check param_2 → ...

#### Problem 2: Not Extensible
New waveform types require:
1. Add new bool flag → API change
2. Add parameters → API change
3. Update loader → Logic change
4. Update validator → Validation change
5. Update docs → Documentation change

**Cost per new waveform type:** ~4 hours of changes

#### Problem 3: Inconsistency
```cpp
// Voltage: Linear sweep (enable_voltage_sweep)
// Voltage: Arbitrary (voltage_time_table)
// Frequency: Linear sweep (enable_frequency_sweep)
// Frequency: Arbitrary (??) - NOT SUPPORTED
```

**Why is voltage_time_table separate from sweep flags?**

#### Problem 4: Copy-Paste Configuration
```json
// Want same sweep for 3 domains? Copy-paste 7 fields × 3!
{
  "domain_1": {
    "AC": {
      "voltage_V": 0,
      "enable_voltage_sweep": true,
      "amplitude_slope_V_s": 500000,
      "start_time_s": 0,
      "rise_time_s": 0.001
    }
  },
  "domain_2": { /* same 5 fields repeated */ },
  "domain_3": { /* same 5 fields repeated */ }
}
```

**Maintenance nightmare:** Change sweep → update 3 places

---

## Design Goals

### Primary Goals

1. **Simplicity:** One field instead of 5-7 fields per parameter
2. **Extensibility:** Add new waveform types without API changes
3. **Consistency:** Same pattern for voltage, frequency, DC fields
4. **Reusability:** Define waveform once, reference multiple times
5. **Backward Compatibility:** Old JSON configs still work

### Secondary Goals

6. **Performance:** Waveform evaluation < 10ns overhead
7. **Type Safety:** Compile-time guarantees for waveform validity
8. **Validation:** Clear error messages for invalid waveforms
9. **Testability:** Easy to unit test individual waveform types
10. **Documentation:** Self-documenting JSON (type field explicit)

### Non-Goals

- **GPU waveforms:** Deferred to Phase 9 (GPU module)
- **Multi-dimensional waveforms:** Not needed for v1.1
- **Real-time waveform editing:** Static config only
- **Waveform DSL:** No custom language, JSON sufficient

---

## Proposed Architecture

### Architecture Overview

```
┌─────────────────────────────────────────────────────────────┐
│                    FieldsConfig (Domain)                    │
├─────────────────────────────────────────────────────────────┤
│  DCFieldConfig:                                             │
│    axial_V: ValueOrWaveform                                 │
│    quad_V: ValueOrWaveform                                  │
│    radial_V: ValueOrWaveform                                │
│                                                              │
│  RFFieldConfig:                                             │
│    voltage_V: ValueOrWaveform                               │
│    frequency_Hz: ValueOrWaveform                            │
│                                                              │
│  ACFieldConfig:                                             │
│    voltage_V: ValueOrWaveform                               │
│    frequency_Hz: ValueOrWaveform                            │
│                                                              │
│  Waveform Library:                                          │
│    std::map<string, Waveform> waveforms                     │
│    ("ac_ramp" → LinearWaveform{...})                        │
│    ("rf_chirp" → LinearWaveform{...})                       │
└─────────────────────────────────────────────────────────────┘
```

### Core Types

#### 1. WaveformType (Enum)

```cpp
// src/core/config/types/WaveformConfig.h

namespace ICARION::config {

enum class WaveformType {
    Constant,       ///< y = value
    Linear,         ///< y = start + (end - start) * t_norm
    Quadratic,      ///< y = a + b*t + c*t²
    Sinusoidal,     ///< y = offset + amplitude * sin(2π*f*t + φ)
    Pulsed,         ///< y = high (during pulse), else low
    Arbitrary       ///< Interpolated from time table
};

} // namespace ICARION::config
```

#### 2. Waveform Structs (Variants)

```cpp
// Constant waveform (default behavior)
struct ConstantWaveform {
    double value = 0.0;
};

// Linear ramp/sweep
struct LinearWaveform {
    double start_value = 0.0;
    double end_value = 0.0;
    double start_time_s = 0.0;
    double duration_s = 0.0;
    bool clamp = true;  ///< Hold end_value after duration
    
    double evaluate(double t) const {
        if (t < start_time_s) return start_value;
        if (t >= start_time_s + duration_s) {
            return clamp ? end_value : start_value;
        }
        double t_norm = (t - start_time_s) / duration_s;
        return start_value + (end_value - start_value) * t_norm;
    }
};

// Quadratic ramp (acceleration profile)
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

// Sinusoidal modulation
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

// Pulsed waveform (single pulse)
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

// Arbitrary waveform (interpolated)
struct ArbitraryWaveform {
    std::vector<double> times_s;
    std::vector<double> values;
    
    enum class Interpolation {
        Linear,     ///< Linear interpolation
        Step,       ///< Step function (hold previous)
        Cubic       ///< Cubic spline (smooth)
    } interp = Interpolation::Linear;
    
    double evaluate(double t) const;  // Implementation in .cpp
};
```

#### 3. Waveform (Variant Container)

```cpp
// Generic waveform using std::variant
using WaveformVariant = std::variant<
    ConstantWaveform,
    LinearWaveform,
    QuadraticWaveform,
    SinusoidalWaveform,
    PulsedWaveform,
    ArbitraryWaveform
>;

struct Waveform {
    std::string id;  ///< Optional ID for named waveforms
    WaveformVariant data;
    
    // Evaluate waveform at time t using std::visit
    double evaluate(double t) const {
        return std::visit([t](const auto& w) { return w.evaluate(t); }, data);
    }
};
```

**Why std::variant?**
- Type-safe (no dynamic_cast)
- Zero overhead (no virtual functions)
- Pattern matching with std::visit
- Compile-time exhaustiveness checking

#### 4. ValueOrWaveform (Flexible Field Type)

```cpp
struct ValueOrWaveform {
    // Option 1: Static value (most common case)
    std::optional<double> constant_value;
    
    // Option 2: Inline waveform definition
    std::optional<Waveform> waveform;
    
    // Option 3: Reference to named waveform
    std::optional<std::string> waveform_ref;
    
    // Validation: Exactly one option set
    bool is_valid() const {
        int count = (constant_value.has_value() ? 1 : 0) +
                    (waveform.has_value() ? 1 : 0) +
                    (waveform_ref.has_value() ? 1 : 0);
        return count == 1;
    }
    
    // Check if time-varying
    bool is_time_varying() const {
        return waveform.has_value() || waveform_ref.has_value();
    }
    
    // Evaluate at time t (requires waveform library for refs)
    double evaluate(double t, const std::map<std::string, Waveform>& library) const;
};
```

### Refactored Field Configs

#### DCFieldConfig (Before & After)

**BEFORE (v1.0):**
```cpp
struct DCFieldConfig {
    double axial_V = 0.0;
    double quad_V = 0.0;
    double radial_V = 0.0;
    double EN_Td = 0.0;
    
    // Future: Would need 12 more fields for sweeps!
    // bool enable_axial_sweep;
    // double axial_slope_V_s;
    // double axial_start_time_s;
    // double axial_rise_time_s;
    // ... same for quad_V and radial_V
};
```

**AFTER (v1.1):**
```cpp
struct DCFieldConfig {
    ValueOrWaveform axial_V;    ///< Can be static OR time-varying
    ValueOrWaveform quad_V;
    ValueOrWaveform radial_V;
    ValueOrWaveform EN_Td;      ///< Alternative specification
    
    ValidationResult validate() const;
};
```

**Reduction:** 4 fields (same) but extensible to sweeps without adding 12 fields!

#### ACFieldConfig (Before & After)

**BEFORE (v1.0):**
```cpp
struct ACFieldConfig {
    double voltage_V = 0.0;
    double frequency_Hz = 0.0;
    
    // Voltage sweep
    bool enable_voltage_sweep = false;
    double amplitude_slope_V_s = 0.0;
    double start_time_s = 0.0;
    double rise_time_s = 0.0;
    
    // Frequency sweep
    bool enable_frequency_sweep = false;
    double frequency_start_Hz = 0.0;
    double frequency_sweep_slope_Hz_s = 0.0;
    
    // Arbitrary (inconsistent)
    std::vector<std::pair<double, double>> voltage_time_table;
    
    // = 12 fields total
};
```

**AFTER (v1.1):**
```cpp
struct ACFieldConfig {
    ValueOrWaveform voltage_V;     ///< Static or time-varying
    ValueOrWaveform frequency_Hz;  ///< Static or chirped
    
    // LQIT phase locking (domain-specific)
    bool lqit_lock_enable = false;
    double lqit_lock_phase_rad = 0.0;
    double lqit_lock_bandwidth_Hz = 0.0;
    
    ValidationResult validate() const;
    
    // = 5 fields (vs 12 before) - 58% reduction!
};
```

**Reduction:** 12 → 5 fields (58% fewer fields!)

#### FieldsConfig (Top-Level)

```cpp
struct FieldsConfig {
    DCFieldConfig dc;
    RFFieldConfig rf;
    ACFieldConfig ac;
    MagneticFieldConfig magnetic;
    
    // NEW: Named waveform library (reusable across fields)
    std::map<std::string, Waveform> waveforms;
    
    // Field arrays (unchanged)
    std::vector<FieldArrayTerm> field_array_terms;
    
    ValidationResult validate() const;
    void compute_derived();
};
```

---

## JSON Schema & Examples

### JSON Schema (schema/waveform.schema.json)

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "$id": "https://icarion.io/schemas/v1.1/waveform.schema.json",
  "title": "ICARION Waveform Schema v1.1",
  
  "definitions": {
    "ConstantWaveform": {
      "type": "object",
      "required": ["type", "value"],
      "properties": {
        "type": {"const": "constant"},
        "value": {"type": "number", "description": "Static value"}
      }
    },
    
    "LinearWaveform": {
      "type": "object",
      "required": ["type", "start", "end", "duration_s"],
      "properties": {
        "type": {"const": "linear"},
        "start": {"type": "number", "description": "Starting value"},
        "end": {"type": "number", "description": "Ending value"},
        "start_time_s": {"type": "number", "default": 0, "minimum": 0},
        "duration_s": {"type": "number", "exclusiveMinimum": 0},
        "clamp": {"type": "boolean", "default": true, "description": "Hold end value after duration"}
      }
    },
    
    "QuadraticWaveform": {
      "type": "object",
      "required": ["type", "a", "b", "c"],
      "properties": {
        "type": {"const": "quadratic"},
        "a": {"type": "number", "description": "Constant term"},
        "b": {"type": "number", "description": "Linear coefficient"},
        "c": {"type": "number", "description": "Quadratic coefficient"},
        "start_time_s": {"type": "number", "default": 0},
        "end_time_s": {"type": "number", "default": 1e9}
      }
    },
    
    "SinusoidalWaveform": {
      "type": "object",
      "required": ["type", "amplitude", "frequency_Hz"],
      "properties": {
        "type": {"const": "sinusoidal"},
        "offset": {"type": "number", "default": 0},
        "amplitude": {"type": "number", "description": "Oscillation amplitude"},
        "frequency_Hz": {"type": "number", "exclusiveMinimum": 0},
        "phase_rad": {"type": "number", "default": 0}
      }
    },
    
    "PulsedWaveform": {
      "type": "object",
      "required": ["type", "low", "high", "pulse_start_s", "pulse_width_s"],
      "properties": {
        "type": {"const": "pulsed"},
        "low": {"type": "number", "description": "Baseline value"},
        "high": {"type": "number", "description": "Pulse value"},
        "pulse_start_s": {"type": "number", "minimum": 0},
        "pulse_width_s": {"type": "number", "exclusiveMinimum": 0}
      }
    },
    
    "ArbitraryWaveform": {
      "type": "object",
      "required": ["type", "times", "values"],
      "properties": {
        "type": {"const": "arbitrary"},
        "times": {
          "type": "array",
          "items": {"type": "number"},
          "minItems": 2,
          "description": "Time points [s], must be sorted"
        },
        "values": {
          "type": "array",
          "items": {"type": "number"},
          "minItems": 2,
          "description": "Values at time points (same length as times)"
        },
        "interpolation": {
          "enum": ["linear", "step", "cubic"],
          "default": "linear"
        }
      }
    },
    
    "Waveform": {
      "oneOf": [
        {"$ref": "#/definitions/ConstantWaveform"},
        {"$ref": "#/definitions/LinearWaveform"},
        {"$ref": "#/definitions/QuadraticWaveform"},
        {"$ref": "#/definitions/SinusoidalWaveform"},
        {"$ref": "#/definitions/PulsedWaveform"},
        {"$ref": "#/definitions/ArbitraryWaveform"}
      ]
    },
    
    "ValueOrWaveform": {
      "oneOf": [
        {
          "type": "number",
          "description": "Static value (shorthand for constant waveform)"
        },
        {
          "$ref": "#/definitions/Waveform",
          "description": "Inline waveform definition"
        },
        {
          "type": "string",
          "pattern": "^@[A-Za-z0-9_]+$",
          "description": "Reference to named waveform (@waveform_id)"
        }
      ]
    }
  }
}
```

### JSON Examples

#### Example 1: Backward Compatible (Static Values)

```json
{
  "fields": {
    "AC": {
      "voltage_V": 100.0,
      "frequency_Hz": 1000000.0
    }
  }
}
```

**Loader behavior:** Converts `100.0` → `ValueOrWaveform{constant_value=100.0}`

#### Example 2: Linear Voltage Ramp (Inline)

```json
{
  "fields": {
    "AC": {
      "voltage_V": {
        "type": "linear",
        "start": 0.0,
        "end": 500.0,
        "start_time_s": 0.0,
        "duration_s": 0.001,
        "clamp": true
      },
      "frequency_Hz": 1000000.0
    }
  }
}
```

**Replaces OLD config:**
```json
{
  "AC": {
    "voltage_V": 0.0,
    "enable_voltage_sweep": true,
    "amplitude_slope_V_s": 500000.0,
    "start_time_s": 0.0,
    "rise_time_s": 0.001
  }
}
```

**Reduction:** 5 fields → 1 field ✅

#### Example 3: Named Waveforms (Reusable)

```json
{
  "waveforms": {
    "ac_voltage_ramp": {
      "type": "linear",
      "start": 0,
      "end": 500,
      "duration_s": 0.001
    },
    "rf_chirp": {
      "type": "linear",
      "start": 1000000,
      "end": 2000000,
      "duration_s": 0.01
    },
    "sinusoidal_modulation": {
      "type": "sinusoidal",
      "offset": 250,
      "amplitude": 50,
      "frequency_Hz": 100
    }
  },
  
  "domains": [
    {
      "name": "domain_1",
      "fields": {
        "AC": {
          "voltage_V": "@ac_voltage_ramp",
          "frequency_Hz": 1000000
        },
        "RF": {
          "voltage_V": "@sinusoidal_modulation",
          "frequency_Hz": "@rf_chirp"
        }
      }
    },
    {
      "name": "domain_2",
      "fields": {
        "AC": {
          "voltage_V": "@ac_voltage_ramp",  
          "frequency_Hz": 500000
        }
      }
    }
  ]
}
```

**Benefits:**
- ✅ Define waveform once, reuse multiple times
- ✅ Named waveforms are self-documenting
- ✅ Change waveform in one place → updates all references

#### Example 4: Arbitrary Waveform (Time Table)

```json
{
  "fields": {
    "AC": {
      "voltage_V": {
        "type": "arbitrary",
        "times": [0, 0.001, 0.002, 0.003, 0.004],
        "values": [0, 100, 50, 200, 0],
        "interpolation": "linear"
      }
    }
  }
}
```

**Backward compatible with:**
```cpp
std::vector<std::pair<double, double>> voltage_time_table = {
    {0.0, 0.0}, {0.001, 100.0}, {0.002, 50.0}, {0.003, 200.0}, {0.004, 0.0}
};
```

#### Example 5: Quadratic Ramp (Acceleration Profile)

```json
{
  "fields": {
    "DC": {
      "axial_V": {
        "type": "quadratic",
        "a": 0,
        "b": 100,
        "c": 50000,
        "start_time_s": 0,
        "end_time_s": 0.01
      }
    }
  }
}
```

**Result:** V(t) = 100t + 50000t² (accelerating ramp)

---

## Implementation Plan

### Phase 1: Waveform Core (4 hours)

**Files to create:**
- `src/core/config/types/WaveformConfig.h` (waveform types & structs)
- `src/core/config/types/WaveformConfig.cpp` (evaluate() implementations)
- `schema/waveform.schema.json` (JSON schema)

**Deliverables:**
- All waveform types (Constant, Linear, Quadratic, Sinusoidal, Pulsed, Arbitrary)
- `Waveform::evaluate(double t)` using std::visit
- `ValueOrWaveform` with validation
- Unit tests: `tests/config/test_waveform_types.cpp` (20 tests)

**Test Coverage:**
- ConstantWaveform::evaluate() (1 test)
- LinearWaveform::evaluate() (5 tests: before start, during, after, clamped, unclamped)
- QuadraticWaveform::evaluate() (3 tests: before, during, after)
- SinusoidalWaveform::evaluate() (4 tests: t=0, T/4, T/2, T)
- PulsedWaveform::evaluate() (3 tests: before, during, after pulse)
- ArbitraryWaveform::evaluate() (4 tests: linear, step, cubic, out of bounds)

**Timeline:** Day 1 (4 hours)

---

### Phase 2: JSON Loader (2 hours)

**Files to create:**
- `src/core/config/loader/WaveformLoader.h` (JSON → Waveform)
- `src/core/config/loader/WaveformLoader.cpp` (implementation)

**Key Functions:**
```cpp
class WaveformLoader {
public:
    // Load single waveform from JSON object
    static Waveform load(const Json::Value& json);
    
    // Load waveform library (map of named waveforms)
    static std::map<std::string, Waveform> load_library(const Json::Value& json);
    
    // Load ValueOrWaveform (handles static value, inline, or reference)
    static ValueOrWaveform load_value_or_waveform(
        const Json::Value& json,
        const std::map<std::string, Waveform>& library
    );
};
```

**Test Coverage:**
- `tests/config/test_waveform_loader.cpp` (15 tests)
  - Load each waveform type from JSON (6 tests)
  - Load ValueOrWaveform (3 tests: static, inline, reference)
  - Load waveform library (2 tests: valid, duplicate IDs)
  - Error handling (4 tests: invalid type, missing fields, bad reference)

**Timeline:** Day 1-2 (2 hours)

---

### Phase 3: Refactor Field Configs (3 hours)

**Files to modify:**
- `src/core/config/types/FieldsConfig.h` (replace booleans with ValueOrWaveform)
- `src/core/config/loader/DomainConfigLoader.cpp` (update parsing)
- `src/core/config/loader/DomainConfigLoader.h` (add waveform library parameter)

**Key Changes:**

**FieldsConfig.h:**
```cpp
struct DCFieldConfig {
    ValueOrWaveform axial_V;  // Was: double axial_V
    ValueOrWaveform quad_V;
    ValueOrWaveform radial_V;
    ValueOrWaveform EN_Td;
};

struct RFFieldConfig {
    ValueOrWaveform voltage_V;
    ValueOrWaveform frequency_Hz;
    double phase_rad = 0.0;
};

struct ACFieldConfig {
    ValueOrWaveform voltage_V;
    ValueOrWaveform frequency_Hz;
    
    // REMOVE:
    // bool enable_voltage_sweep;
    // double amplitude_slope_V_s;
    // ... (all sweep flags deleted)
};
```

**DomainConfigLoader.cpp:**
```cpp
// OLD:
if (json["AC"]["voltage_V"].isNumeric()) {
    ac.voltage_V = json["AC"]["voltage_V"].asDouble();
}

// NEW:
ac.voltage_V = WaveformLoader::load_value_or_waveform(
    json["AC"]["voltage_V"],
    fields_config.waveforms  // Library from top-level
);
```

**Backward Compatibility:**
```cpp
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
        if (ref[0] != '@') throw std::runtime_error("Waveform ref must start with @");
        result.waveform_ref = ref.substr(1);  // Remove '@'
    } else if (json.isObject()) {
        // NEW FORMAT: Inline waveform
        result.waveform = WaveformLoader::load(json);
    } else {
        throw std::runtime_error("Invalid ValueOrWaveform format");
    }
    
    return result;
}
```

**Test Coverage:**
- Update existing config loader tests (backward compatibility)
- Add new tests for waveform loading
- `tests/config/test_domain_config_waveforms.cpp` (10 tests)

**Timeline:** Day 2-3 (3 hours)

---

### Phase 4: Update Physics Code (2 hours)

**Files to modify:**
- `src/core/physics/forces/ElectricFieldForce.cpp`
- `src/core/physics/forces/RFFieldForce.cpp`
- `src/core/physics/forces/ACFieldForce.cpp` (if exists)

**Key Change: Evaluate Waveforms at Runtime**

**OLD:**
```cpp
// In ElectricFieldForce::compute(t)
double voltage = domain_.fields.ac.voltage_V;

if (domain_.fields.ac.enable_voltage_sweep) {
    double t_rel = t - domain_.fields.ac.start_time_s;
    if (t_rel >= 0 && t_rel <= domain_.fields.ac.rise_time_s) {
        voltage += domain_.fields.ac.amplitude_slope_V_s * t_rel;
    }
}
```

**NEW:**
```cpp
// Helper function (in FieldsConfig.h or utils)
double evaluate_value_or_waveform(
    const ValueOrWaveform& val,
    double t,
    const std::map<std::string, Waveform>& library
) {
    if (val.constant_value.has_value()) {
        return *val.constant_value;
    } else if (val.waveform.has_value()) {
        return val.waveform->evaluate(t);
    } else if (val.waveform_ref.has_value()) {
        auto it = library.find(*val.waveform_ref);
        if (it == library.end()) {
            throw std::runtime_error("Waveform reference not found: " + *val.waveform_ref);
        }
        return it->second.evaluate(t);
    }
    return 0.0;  // Fallback
}

// In ElectricFieldForce::compute(t)
double voltage = evaluate_value_or_waveform(
    domain_.fields.ac.voltage_V,
    t,
    domain_.fields.waveforms
);
```

**Test Coverage:**
- `tests/physics/forces/test_field_force_waveforms.cpp` (8 tests)
  - Static voltage (backward compat)
  - Linear voltage sweep
  - Sinusoidal modulation
  - Waveform reference
  - Verify force vs time curve

**Timeline:** Day 3 (2 hours)

---

### Phase 5: Migration Tools (2 hours)

**Files to create:**
- `scripts/migrate_to_waveforms.py` (auto-convert old configs)
- `docs/MIGRATION_WAVEFORMS.md` (migration guide)

**Migration Script:**
```python
#!/usr/bin/env python3
"""
Convert v1.0 configs (sweep flags) to v1.1 configs (waveforms)
"""

import json
import sys

def migrate_ac_field(ac_config):
    """Convert AC field with sweep flags to waveform format"""
    result = {}
    
    # Voltage
    if ac_config.get("enable_voltage_sweep", False):
        # Convert to linear waveform
        start = ac_config.get("voltage_V", 0.0)
        slope = ac_config.get("amplitude_slope_V_s", 0.0)
        rise_time = ac_config.get("rise_time_s", 0.0)
        start_time = ac_config.get("start_time_s", 0.0)
        
        end = start + slope * rise_time
        
        result["voltage_V"] = {
            "type": "linear",
            "start": start,
            "end": end,
            "start_time_s": start_time,
            "duration_s": rise_time,
            "clamp": True
        }
    elif "voltage_time_table" in ac_config:
        # Convert to arbitrary waveform
        table = ac_config["voltage_time_table"]
        result["voltage_V"] = {
            "type": "arbitrary",
            "times": [t for t, v in table],
            "values": [v for t, v in table],
            "interpolation": "linear"
        }
    else:
        # Static value
        result["voltage_V"] = ac_config.get("voltage_V", 0.0)
    
    # Frequency
    if ac_config.get("enable_frequency_sweep", False):
        start = ac_config.get("frequency_start_Hz", 0.0)
        slope = ac_config.get("frequency_sweep_slope_Hz_s", 0.0)
        # ... similar conversion
    else:
        result["frequency_Hz"] = ac_config.get("frequency_Hz", 0.0)
    
    return result

def migrate_config(old_config):
    """Main migration function"""
    new_config = old_config.copy()
    
    for domain in new_config.get("domains", []):
        if "fields" in domain:
            fields = domain["fields"]
            
            if "AC" in fields:
                fields["AC"] = migrate_ac_field(fields["AC"])
            
            # Similar for RF, DC if they have sweeps
    
    return new_config

if __name__ == "__main__":
    with open(sys.argv[1], "r") as f:
        old = json.load(f)
    
    new = migrate_config(old)
    
    with open(sys.argv[2], "w") as f:
        json.dump(new, f, indent=2)
    
    print(f"✓ Migrated {sys.argv[1]} → {sys.argv[2]}")
```

**Migration Guide (docs/MIGRATION_WAVEFORMS.md):**
```md
# Waveform System Migration Guide (v1.0 → v1.1)

## Automatic Migration

```bash
python3 scripts/migrate_to_waveforms.py config_v1.0.json config_v1.1.json
```

## Manual Migration

### OLD (v1.0): Voltage Sweep with Flags
```json
{
  "AC": {
    "voltage_V": 0,
    "enable_voltage_sweep": true,
    "amplitude_slope_V_s": 500000,
    "start_time_s": 0,
    "rise_time_s": 0.001
  }
}
```

### NEW (v1.1): Linear Waveform
```json
{
  "AC": {
    "voltage_V": {
      "type": "linear",
      "start": 0,
      "end": 500,
      "duration_s": 0.001
    }
  }
}
```

## Backward Compatibility

Old configs (v1.0) will continue to work in v1.1:
- Static values: `"voltage_V": 100` → auto-converted to ConstantWaveform
- Sweep flags: Will be supported until v2.0 (deprecated warning)
```

**Timeline:** Day 3-4 (2 hours)

---

### Phase 6: Documentation (2 hours)

**Files to update:**
- `docs/CONFIG_GUIDE.md` (add Waveform section)
- `docs/ARCHITECTURE.md` (update Field Config section)
- `examples/waveforms/` (example configs)

**CONFIG_GUIDE.md additions:**
```md
### Waveforms (Time-Varying Parameters)

ICARION supports time-varying field parameters through waveforms.

#### Static Values (Backward Compatible)

```json
{"voltage_V": 100.0}
```

#### Linear Ramp

```json
{
  "voltage_V": {
    "type": "linear",
    "start": 0,
    "end": 500,
    "duration_s": 0.001
  }
}
```

#### Named Waveforms (Reusable)

```json
{
  "waveforms": {
    "my_ramp": {"type": "linear", "start": 0, "end": 500, "duration_s": 0.001}
  },
  "fields": {
    "AC": {"voltage_V": "@my_ramp"}
  }
}
```

#### Available Waveform Types

| Type | Use Case | Parameters |
|------|----------|------------|
| `constant` | Static value | `value` |
| `linear` | Voltage/frequency sweep | `start`, `end`, `duration_s` |
| `quadratic` | Acceleration profile | `a`, `b`, `c` |
| `sinusoidal` | Amplitude modulation | `amplitude`, `frequency_Hz` |
| `pulsed` | Single pulse | `low`, `high`, `pulse_start_s`, `pulse_width_s` |
| `arbitrary` | Custom curve | `times[]`, `values[]` |
```

**Example Configs:**
- `examples/waveforms/linear_voltage_ramp.json`
- `examples/waveforms/frequency_chirp.json`
- `examples/waveforms/amplitude_modulation.json`
- `examples/waveforms/reusable_waveforms.json`

**Timeline:** Day 4-5 (2 hours)

---

## Total Implementation Time: 11 hours

| Phase | Task | Hours |
|-------|------|-------|
| 1 | Waveform Core Types | 4 |
| 2 | JSON Loader | 2 |
| 3 | Refactor Field Configs | 3 |
| 4 | Update Physics Code | 2 |
| 5 | Migration Tools | 2 |
| 6 | Documentation | 2 |
| **Total** | **Full Implementation** | **11** |

**Timeline:** 5 working days (spread over 1-2 weeks)

---

## Migration Strategy

### Backward Compatibility Guarantee

**Principle:** All v1.0 configs work in v1.1 without changes

**Implementation:**
1. **Static values auto-convert:**
   ```cpp
   // JSON: "voltage_V": 100
   // Loader creates: ValueOrWaveform{constant_value=100}
   ```

2. **Sweep flags supported (deprecated):**
   - Keep `enable_voltage_sweep` logic in loader
   - Emit deprecation warning: "Sweep flags deprecated, use waveforms"
   - Convert to LinearWaveform internally

3. **Gradual migration:**
   - v1.1: Both formats supported
   - v1.2: Sweep flags deprecated (warnings)
   - v2.0: Sweep flags removed

### Migration Workflow

**Step 1: Automatic Migration (Recommended)**
```bash
cd /path/to/configs
python3 $ICARION_ROOT/scripts/migrate_to_waveforms.py old.json new.json
```

**Step 2: Manual Review**
- Check generated waveforms match intent
- Consider named waveforms for reusability
- Test with `--validate-config` flag

**Step 3: Gradual Rollout**
- Migrate examples/ first
- Migrate tests/ configs
- User configs at their convenience

---

## Testing Strategy

### Unit Tests (20 tests)

**File:** `tests/config/test_waveform_types.cpp`

```cpp
TEST_CASE("ConstantWaveform evaluates to constant", "[waveform]") {
    ConstantWaveform w{100.0};
    CHECK(w.evaluate(0.0) == 100.0);
    CHECK(w.evaluate(1.0) == 100.0);
}

TEST_CASE("LinearWaveform ramps correctly", "[waveform]") {
    LinearWaveform w{0.0, 100.0, 0.0, 0.1, true};
    CHECK(w.evaluate(0.0) == Approx(0.0));
    CHECK(w.evaluate(0.05) == Approx(50.0));
    CHECK(w.evaluate(0.1) == Approx(100.0));
    CHECK(w.evaluate(0.2) == Approx(100.0));  // Clamped
}

TEST_CASE("SinusoidalWaveform oscillates", "[waveform]") {
    SinusoidalWaveform w{250.0, 50.0, 100.0, 0.0};
    CHECK(w.evaluate(0.0) == Approx(250.0));
    CHECK(w.evaluate(0.0025) == Approx(300.0).margin(1.0));  // Peak
    CHECK(w.evaluate(0.005) == Approx(250.0).margin(1.0));   // Zero
}

// ... 17 more tests
```

### Loader Tests (15 tests)

**File:** `tests/config/test_waveform_loader.cpp`

```cpp
TEST_CASE("Load linear waveform from JSON", "[loader]") {
    Json::Value json;
    json["type"] = "linear";
    json["start"] = 0.0;
    json["end"] = 100.0;
    json["duration_s"] = 0.1;
    
    Waveform w = WaveformLoader::load(json);
    CHECK(std::holds_alternative<LinearWaveform>(w.data));
}

TEST_CASE("Load ValueOrWaveform: static value", "[loader]") {
    Json::Value json = 100.0;
    std::map<std::string, Waveform> lib;
    
    ValueOrWaveform val = WaveformLoader::load_value_or_waveform(json, lib);
    CHECK(val.constant_value.has_value());
    CHECK(*val.constant_value == 100.0);
}

TEST_CASE("Load ValueOrWaveform: waveform reference", "[loader]") {
    Json::Value json = "@my_waveform";
    std::map<std::string, Waveform> lib;
    lib["my_waveform"] = /* ... */;
    
    ValueOrWaveform val = WaveformLoader::load_value_or_waveform(json, lib);
    CHECK(val.waveform_ref.has_value());
    CHECK(*val.waveform_ref == "my_waveform");
}

// ... 12 more tests
```

### Integration Tests (8 tests)

**File:** `tests/integration/test_field_force_waveforms.cpp`

```cpp
TEST_CASE("AC voltage linear ramp produces expected force", "[integration]") {
    // Setup config with linear voltage waveform
    // Run simulation for 1ms
    // Verify force increases linearly
}

TEST_CASE("RF frequency chirp changes field frequency", "[integration]") {
    // Setup RF with frequency waveform
    // Sample E-field at multiple times
    // Verify frequency increases as expected
}
```

### E2E Tests (5 tests)

**File:** `tests/e2e/test_waveform_simulations.cpp`

```cpp
TEST_CASE("IMS with voltage ramp: ions arrive later", "[e2e]") {
    // Run IMS with static 100V
    // Run IMS with 0→100V ramp
    // Compare arrival times (ramped should be later)
}
```

### Total Test Coverage

| Category | Count | File |
|----------|-------|------|
| Unit tests | 20 | test_waveform_types.cpp |
| Loader tests | 15 | test_waveform_loader.cpp |
| Integration tests | 8 | test_field_force_waveforms.cpp |
| E2E tests | 5 | test_waveform_simulations.cpp |
| **Total** | **48** | |

---

## Performance Considerations

### Waveform Evaluation Overhead

**Target:** < 10ns per waveform evaluation

**Benchmark:**
```cpp
// tests/benchmark/bench_waveform_eval.cpp
BENCHMARK("ConstantWaveform evaluate") {
    ConstantWaveform w{100.0};
    return w.evaluate(0.5);
};  // Expected: ~2ns

BENCHMARK("LinearWaveform evaluate") {
    LinearWaveform w{0, 100, 0, 1, true};
    return w.evaluate(0.5);
};  // Expected: ~5ns

BENCHMARK("SinusoidalWaveform evaluate") {
    SinusoidalWaveform w{250, 50, 100, 0};
    return w.evaluate(0.5);
};  // Expected: ~10ns (sin() call)
```

**Optimization strategies:**
- Inline evaluate() methods (compiler hint)
- Cache frequently used waveforms
- Pre-compute sinusoidal coefficients

### Memory Overhead

**Per ValueOrWaveform:**
- std::optional<double>: 16 bytes
- std::optional<Waveform>: ~80 bytes (variant + largest waveform type)
- std::optional<std::string>: ~32 bytes

**Total:** ~128 bytes (vs 8 bytes for double)

**Impact:** Negligible for config (only a few dozen fields total)

### JSON Parsing Overhead

**Waveform library parsing:** One-time cost at config load (~1ms for 100 waveforms)

**Acceptable:** Config loading is not performance-critical

---

## Future Extensions

### Phase 1.2 Enhancements (Future)

1. **Waveform Composition:**
   ```json
   {
     "type": "sum",
     "waveforms": ["@waveform_1", "@waveform_2"]
   }
   ```

2. **Conditional Waveforms:**
   ```json
   {
     "type": "conditional",
     "condition": "ion_count > 1000",
     "if_true": "@waveform_high_N",
     "if_false": "@waveform_low_N"
   }
   ```

3. **Multi-Dimensional Waveforms:**
   ```json
   {
     "type": "vector",
     "x": "@waveform_x",
     "y": "@waveform_y",
     "z": "@waveform_z"
   }
   ```

4. **GPU Waveform Kernels:**
   - Evaluate waveforms on GPU for large ensembles
   - Pre-compute waveform table for entire simulation

### Phase 2.0 Enhancements (Far Future)

5. **Waveform DSL:**
   ```json
   {
     "type": "expression",
     "formula": "100 * sin(2*pi*freq*t) + 50"
   }
   ```

6. **Real-Time Waveform Updates:**
   - Modify waveforms during simulation
   - Interactive parameter sweeps

---

## Decision Log

### Decision 1: std::variant vs Inheritance

**Options:**
- A: Use std::variant<ConstantWaveform, LinearWaveform, ...>
- B: Use abstract base class + virtual functions

**Decision:** Option A (std::variant)

**Rationale:**
- Zero overhead (no virtual dispatch)
- Type-safe (compile-time checks)
- Pattern matching with std::visit
- Easier to extend (just add to variant)

**Downside:** Cannot add new waveform types at runtime (acceptable for v1.1)

---

### Decision 2: Inline vs Named Waveforms

**Options:**
- A: Only inline waveforms (simpler)
- B: Support both inline and named (reusable)

**Decision:** Option B (support both)

**Rationale:**
- Reusability is critical for multi-domain configs
- Named waveforms self-documenting
- Low implementation cost (~30 LOC)

**Trade-off:** Slightly more complex loader

---

### Decision 3: Waveform Reference Syntax

**Options:**
- A: String reference: `"voltage_V": "my_waveform"`
- B: Object reference: `"voltage_V": {"ref": "my_waveform"}`
- C: @ prefix: `"voltage_V": "@my_waveform"`

**Decision:** Option C (@ prefix)

**Rationale:**
- Visually distinct from static values
- JSON schema can validate with regex: `"^@[A-Za-z0-9_]+$"`
- Familiar pattern (CSS variables, shell, etc.)

**Example:**
```json
{"voltage_V": "@my_waveform"}  // Clear it's a reference
```

---

### Decision 4: Backward Compatibility Strategy

**Options:**
- A: Breaking change (remove sweep flags in v1.1)
- B: Deprecate sweep flags (remove in v2.0)
- C: Dual support forever

**Decision:** Option B (deprecate in v1.1, remove in v2.0)

**Rationale:**
- Gives users time to migrate (1-2 release cycles)
- No surprise breakage
- Clean slate in v2.0

**Migration period:** v1.1 → v2.0 (~6-12 months)

---

## Open Questions

### Question 1: GPU Waveform Evaluation

**Context:** Phase 9 (GPU acceleration) needs waveforms on GPU

**Options:**
- A: Copy waveform parameters to GPU, evaluate in kernels
- B: Pre-compute waveform table on CPU, transfer to GPU
- C: Hybrid: Simple waveforms on GPU, complex on CPU

**Decision:** DEFERRED to Phase 9 design

**Note:** Option B likely best (amortize evaluation cost)

---

### Question 2: Waveform Validation Strictness

**Context:** Should we validate waveform parameters strictly?

**Example:**
```json
{
  "type": "linear",
  "start": 0,
  "end": -100,  // Negative end value - error or warning?
  "duration_s": 0.001
}
```

**Options:**
- A: Strict validation (end < start is error)
- B: Lenient validation (negative values allowed)

**Decision:** Option B (lenient)

**Rationale:**
- Physical parameters can be negative (voltages, frequencies in some contexts)
- User knows their domain better than validator
- Validation should catch typos, not physics

---

## Approval Checklist

Before implementation:
- [ ] Design reviewed by architecture team
- [ ] JSON schema validated with JSON Schema validator
- [ ] Example configs tested manually
- [ ] Performance targets agreed upon
- [ ] Migration strategy approved
- [ ] Documentation plan reviewed
- [ ] Test coverage plan approved

---

## References

- **Issue:** #TODO (create GitHub issue for waveform system)
- **Related Design Docs:**
  - `tmp/SSOT_MIGRATION_SUMMARY.md` (SSOT principles)
  - `docs/ARCHITECTURE.md` (field config architecture)
  - `docs/CONFIG_GUIDE.md` (user-facing config docs)
- **Implementation Branch:** `feature/waveform-system`
- **Target Milestone:** v1.1 (post Phase 6)

---

**END OF DESIGN DOCUMENT**

---

## Appendix A: Full Type Definitions

(See [Proposed Architecture](#proposed-architecture) for complete C++ code)

## Appendix B: JSON Schema (Complete)

(See [JSON Schema & Examples](#json-schema--examples) for complete schema)

## Appendix C: Migration Script (Complete)

(See [Phase 5: Migration Tools](#phase-5-migration-tools-2-hours) for complete script)

## Appendix D: Performance Benchmarks

| Waveform Type | Evaluation Time (ns) | Memory (bytes) |
|---------------|----------------------|----------------|
| Constant | 2 | 8 |
| Linear | 5 | 40 |
| Quadratic | 8 | 48 |
| Sinusoidal | 10 | 40 |
| Pulsed | 6 | 32 |
| Arbitrary (linear) | 15 | 64 + N×16 |

**Notes:**
- Benchmarked on Intel i7-12700K @ 3.6 GHz
- Compiled with `-O3 -march=native`
- Cache-hot measurements (realistic simulation scenario)

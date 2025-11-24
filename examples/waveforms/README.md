# ICARION Waveform Examples

This directory contains example configurations demonstrating the waveform system for time-varying field parameters.

## Available Examples

### 1. Linear Voltage Ramp (`linear_voltage_ramp.json`)

IMS drift tube with linear voltage increase from 0V to 500V over 1ms.

**Features:**
- Linear waveform type
- Clamped after duration (holds 500V)
- Use case: Voltage scanning in IMS

**Run:**
```bash
./build/icarion_cli examples/waveforms/linear_voltage_ramp.json
```

### 2. Frequency Chirp (`frequency_chirp.json`)

RF quadrupole with linear frequency sweep from 1MHz to 2MHz over 10ms.

**Features:**
- Linear frequency chirp
- RF field with time-varying frequency
- Use case: Mass scanning, resonant excitation

**Run:**
```bash
./build/icarion_cli examples/waveforms/frequency_chirp.json
```

### 3. Amplitude Modulation (`amplitude_modulation.json`)

AC field with sinusoidal amplitude modulation at 100Hz.

**Features:**
- Sinusoidal waveform
- Offset + amplitude oscillation
- Use case: Parametric excitation, resonance studies

**Run:**
```bash
./build/icarion_cli examples/waveforms/amplitude_modulation.json
```

### 4. Reusable Named Waveforms (`reusable_waveforms.json`)

Multi-domain simulation with shared waveforms defined in library.

**Features:**
- Top-level `waveforms` library
- Multiple domains referencing same waveforms
- Use case: Complex instruments with synchronized waveforms

**Run:**
```bash
./build/icarion_cli examples/waveforms/reusable_waveforms.json
```

### 5. Arbitrary Waveform (`arbitrary_waveform.json`)

Custom voltage profile from time-value pairs with linear interpolation.

**Features:**
- Arbitrary waveform type
- User-defined time points and values
- Use case: Experimental waveform reproduction

**Run:**
```bash
./build/icarion_cli examples/waveforms/arbitrary_waveform.json
```

## Waveform Types Summary

| Type | Parameters | Use Case |
|------|------------|----------|
| `linear` | start, end, duration_s | Voltage/frequency sweeps |
| `quadratic` | a, b, c | Acceleration profiles |
| `sinusoidal` | amplitude, frequency_Hz | Modulation, oscillation |
| `pulsed` | low, high, pulse_start_s, pulse_width_s | Gating, injection |
| `arbitrary` | times[], values[] | Custom experimental profiles |
| `constant` | value | Static value (backward compatible) |

## Creating Custom Waveforms

### Inline Waveform

```json
"fields": {
  "DC": {
    "axial_V": {
      "type": "linear",
      "start": 0,
      "end": 500,
      "duration_s": 0.001
    }
  }
}
```

### Named Waveform Library

```json
{
  "waveforms": {
    "my_ramp": {
      "type": "linear",
      "start": 0,
      "end": 500,
      "duration_s": 0.001
    }
  },
  "domains": [{
    "fields": {
      "DC": {
        "axial_V": "@my_ramp"
      }
    }
  }]
}
```

## Validation

Validate waveform configurations:

```bash
# JSON schema validation
python3 schema/validate_config.py examples/waveforms/linear_voltage_ramp.json

# Waveform unit tests
cd build && ctest -R Waveform
```

## Documentation

See `docs/CONFIG_GUIDE.md` section "Waveforms (Time-Varying Parameters)" for complete documentation.

## Migration from v1.0

Convert old sweep flags to new waveform format:

```bash
python3 scripts/migrate_to_waveforms.py old_config.json new_config.json
```

**Deprecated (v1.0):**
```json
{
  "AC": {
    "voltage_V": 0,
    "enable_voltage_sweep": true,
    "amplitude_slope_V_s": 500000,
    "rise_time_s": 0.001
  }
}
```

**New (v1.1):**
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

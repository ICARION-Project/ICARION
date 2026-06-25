# Fields and waveforms

Fields are configured per domain under `domains[].fields`. ICARION v1.0.1
supports analytical DC/RF/AC fields, magnetic fields, and imported field array
terms.

```json
"fields": {
  "DC": {
    "EN_Td": 10.0
  },
  "RF": {
    "voltage_V": 620.0,
    "frequency_Hz": 1000000.0,
    "phase_rad": 0.0
  }
}
```

The input keys are uppercase: `DC`, `RF`, and `AC`. In the HDF5 output these are
stored as lowercase groups under `/domains/<name>/fields/dc`, `rf`, and `ac`.

## DC fields

`DC` accepts:

- `EN_Td`: reduced electric field in Townsend.
- `axial_V`: axial voltage.
- `radial_V`: radial voltage (for Orbitraps).
- `quad_V`: quadrupolar DC voltage (between adjacent electrodes).

If both `EN_Td` and `axial_V` are given, `axial_V` takes precedence. `EN_Td` is
converted using the domain gas pressure, temperature, and length. This is
convenient for IMS examples, but it means that changing the domain environment
also changes the converted axial field.

## RF and AC fields

`RF` accepts:

- `voltage_V`: RF amplitude.
- `frequency_Hz`: RF frequency.
- `phase_rad`: phase offset in radians.

`AC` accepts:

- `voltage_V`: AC excitation amplitude.
- `frequency_Hz`: AC excitation frequency.

The exact force interpretation also depends on the instrument geometry. For
example, an RF voltage in an LQIT or quadrupole is not equivalent to an axial IMS
drift voltage.

## Magnetic fields

Magnetic fields are configured under `B`:

```json
"fields": {
  "B": {
    "enabled": true,
    "field_strength_T": [0.0, 0.0, 3.0]
  }
}
```

This is mainly relevant for FT-ICR examples. If `enabled` is false or
omitted, the magnetic field contribution is inactive.

## Waveforms

Field values can be static numbers, inline waveform objects, or references to a
named waveform using `@waveform_id`.

```json
"waveforms": {
  "voltage_ramp": {
    "type": "linear",
    "start": 0.0,
    "end": 500.0,
    "end_time_s": 0.001
  }
},
"domains": [
  {
    "name": "drift",
    "fields": {
      "DC": {
        "axial_V": "@voltage_ramp"
      }
    }
  }
]
```

Supported waveform families include:

- `constant`
- `linear`
- `quadratic`
- `exponential`
- `sinusoidal`
- `pwm`
- `pulsed`
- `arbitrary`

Named waveforms can be declared globally or domainl-ocally. Domain-local
waveforms are useful when two domains use the same waveform name for different
physical meanings.

## Field arrays

Precomputed field arrays are configured with `field_array_terms`:

```json
"fields": {
  "field_array_terms": [
    {
      "file": "examples/field_arrays/ims_drift_200Vcm.h5",
      "scale_type": "Constant",
      "constant_V": 300.0
    }
  ]
}
```

Supported scale types are:

- `Constant`
- `DC_Axial`
- `DC_Quad`
- `DC_Radial`
- `RF`

Field array files are treated as external input data. Keep them under version
control or archive them with the configuration and output file if the simulation
must be reproducible.

## Practical caveats

- Fields are domain-local; moving an ion into another domain changes the active
  domain field model.
- `physics.collision_model` is global, not domain-local.
- Time-dependent fields require a timestep small enough to resolve the waveform.
- HDF5 field metadata stores static or `t=0` values plus waveform metadata; it is
  not a replacement for the original JSON configuration.
- Deprecated sweep fields still exist in the schema, but new examples should use
  waveforms instead.

Start with `examples/ims/ims_basic.json`, then compare with
`examples/waveforms/linear_voltage_ramp.json` and
`examples/ims/ims_field_array_multi_domain.json`.

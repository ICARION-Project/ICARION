# Fields and waveforms

Fields are configured per domain under `domains[].fields`. ICARION v1.1.0
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

Here, `start` and `end` decode the values of the applied axial voltage at the start time or end time, specified by `end_time_s`.

Supported waveform families include:

- `constant`
- `linear`
- `quadratic`
- `exponential`
- `sinusoidal`
- `pwm`
- `pulsed`
- `arbitrary`

Named waveforms can be declared globally or inside individual domains. Domain-local
waveforms are useful when two domains use the same waveform name for different
physical meanings or settings.

## Instrument field quick reference

The following table summarizes the analytical field parameters most commonly
used by the supported instrument types. Field arrays can replace or supplement
these analytical fields when imported data are available.

| Instrument | Typical field parameters | Meaning |
| --- | --- | --- |
| `IMS` | `DC.EN_Td` or `DC.axial_V` | Drift field along the domain axis; `EN_Td` is converted using pressure, temperature, and domain length. |
| `TIMS` | `TIMS.axial_field_initial_*`, `TIMS.axial_field_final_*`, optional gas flow | Trapped IMS axial ramp field. The runtime interpolates from an initial axial profile to a final axial profile over the configured ramp. |
| `TOF` | `DC.axial_V` | Acceleration voltage over `geometry.acc_length_m`; field is zero outside the acceleration region. |
| `Quadrupole` / `QuadrupoleRF` | `RF.voltage_V`, `RF.frequency_Hz`, optional `DC.quad_V`, optional `DC.axial_V` | Quadrupolar RF/DC confinement plus optional axial potential gradient. |
| `LQIT` | `RF.voltage_V`, `RF.frequency_Hz`, optional `DC.axial_V`, optional `AC.voltage_V`/`AC.frequency_Hz` | Linear ion trap RF confinement, (approximated) harmonic axial DC, and optional AC excitation (in Cartesian `x` coordinate). |
| `Orbitrap` | `DC.radial_V` | Orbitrap radial potential; requires Orbitrap-specific radii (see below). |
| `FTICR` | `B.field_strength_T`, optional `DC.radial_V` | Uniform magnetic field for cyclotron motion plus optional trapping DC. |

`NoFixedInstrument` is also accepted for custom domains. In that case,
the active fields depend on the configured field model rather than a named
instrument preset.

## Geometry quick reference

Every domain requires `geometry.origin_m` and `geometry.length_m`. Most
cylindrical examples also use `geometry.radius_m`. Some instruments need
additional geometry parameters:

| Instrument | Geometry parameters to check | Notes |
| --- | --- | --- |
| `IMS` | `origin_m`, `length_m`, `radius_m` | `length_m` is used when converting `EN_Td` to `axial_V`. |
| `TIMS` | `origin_m`, `length_m`, `radius_m` | `length_m` defines the axial ramp coordinate; TIMS examples commonly use axial gas flow. |
| `TOF` | `origin_m`, `length_m`, `radius_m`, `acc_length_m`, optional `end_aperture_m` | `acc_length_m` defines the acceleration section. |
| `Quadrupole` / `QuadrupoleRF` | `origin_m`, `length_m`, `radius_m` | `radius_m` is the effective quadrupole radius used in the analytical RF/DC field. |
| `LQIT` | `origin_m`, `length_m`, `radius_m` | `radius_m` sets the radial field scale; `length_m` is used for harmonic axial DC terms. |
| `Orbitrap` | `origin_m`, `length_m`, `radius_in_m`, `radius_out_m`, `radius_char_m` | `radius_in_m` and `radius_out_m` are the inner/outer electrode radii at `z = 0`, not constant radii along the full electrode length; `radius_char_m` is the characteristic Orbitrap radius controlling the hyperlogarithmic field curvature. |
| `FTICR` | `origin_m`, `length_m`, `radius_m` | `length_m` and `radius_m` set the trapping DC scale. |

These are practical starting points, not a calibration guide. Always check the
example closest to your instrument and run timestep and output sanity checks before
interpreting physics results.

## TIMS axial ramp fields

TIMS domains can use `instrument: "TIMS"` with a `fields.TIMS` block. The
analytical model evaluates an axial field profile of the form:

```text
E_z(z,t) = (1 - f(t)) * E_initial(z) + f(t) * E_final(z)
```

where `f(t)` is the ramp fraction. Profiles may be configured as uniform values
or tabulated axial points, depending on the example and schema. TIMS examples
usually combine this field with an axial gas flow profile in the domain
environment.

Start from `examples/ims/ims_tims_basic.json` and inspect the validation case
before using TIMS fields for production studies.

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

Multiple field array terms are treated as a linear superposition:

```text
E_total(r, t) = sum_i scale_i(t) * E_i(r)
```

This makes it possible to load separate basis field maps, for example one DC
map and one RF map from a BEM/FEM solver, and let ICARION scale and sum them at
runtime. A common convention is that each imported field array represents a
unit voltage solution and is scaled with one of the factors below.

Supported scale types are:

| Scale type | Meaning |
| --- | --- |
| `Constant` | Multiplies the imported field term by `constant_V`, independent of analytical DC/RF settings. |
| `DC_Axial` | Scales the imported term with the configured axial DC voltage. |
| `DC_Quad` | Scales the imported term with the configured quadrupolar DC voltage. |
| `DC_Radial` | Scales the imported term with the configured radial DC voltage. |
| `RF` | Scales the imported term as `RF.voltage_V(t) * cos(2π f t + phase_rad)`. If `field_array_terms[].frequency_Hz` is greater than zero, it is used for `f`; otherwise ICARION uses `fields.RF.frequency_Hz(t)`. |

For `RF` terms, the amplitude comes from `fields.RF.voltage_V`, not from
`constant_V`. Use `constant_V` only with `scale_type: "Constant"`.

Field array files are treated as external input data. We suggest to keep them under version
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

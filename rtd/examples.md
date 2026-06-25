# Examples

The `examples/` directory contains ready to run configurations. Use them as
known good starting points before writing a configuration from scratch.

From the repository root:

```bash
./build/src/icarion_main examples/ims/ims_basic.json
```

If ICARION is installed:

```bash
icarion examples/ims/ims_basic.json
```

## Recommended first path

1. Run `examples/ims/ims_basic.json` unchanged.
2. Inspect the generated HDF5 file with an analysis script.
3. Run the same file with `--dry-run` and `--validate-config`.
4. Copy the JSON file and change one physical parameter.
5. Compare output against the unchanged baseline.

## IMS examples

| File | Use case |
| --- | --- |
| `examples/ims/ims_basic.json` | Basic two domain (drift + detector region) IMS drift setup. |
| `examples/ims/ims_multi_gas_air.json` | Multi-gas IMS environment. |
| `examples/ims/ims_basic_ion_cloud_file.json` | Ion cloud loaded from file. |
| `examples/ims/ims_with_field_array.json` | IMS with imported field array data. |
| `examples/ims/ims_field_array_time_varying.json` | Time-varying field array scaling. |
| `examples/ims/ims_field_array_multi_domain.json` | Multi-domain field array setup. |

Start here if you care about mobility, gas composition, field arrays, or
arrival time analysis.

## Instrument examples

| File | Use case |
| --- | --- |
| `examples/tof/tof_basic.json` | Basic time-of-flight configuration. |
| `examples/tof/tof_reflectron.json` | Reflectron TOF setup. |
| `examples/orbitrap/orbitrap_basic.json` | Orbitrap electrostatic trapping. |
| `examples/lqit/lqit_basic.json` | Linear quadrupole ion trap. |
| `examples/quadrupole/quadrupole_basic.json` | Quadrupole RF/DC fields. |
| `examples/fticr/fticr_basic.json` | FT-ICR magnetic field example. |

These examples are best used to understand configuration structure and output
format. Validate physics assumptions before treating them as calibrated
instrument models.

## Reactions

| File | Use case |
| --- | --- |
| `examples/reactions/reaction_demo.json` | Minimal stochastic reaction workflow. |

Use this with [Reactions](reactions.md) and [Species database](species-database.md).
Product species must exist in the species database.

## Waveforms

| File | Use case |
| --- | --- |
| `examples/waveforms/linear_voltage_ramp.json` | Linear voltage ramp. |
| `examples/waveforms/frequency_chirp.json` | Time-dependent RF frequency. |
| `examples/waveforms/amplitude_modulation.json` | Sinusoidal amplitude modulation. |
| `examples/waveforms/arbitrary_waveform.json` | Table-driven arbitrary waveform. |
| `examples/waveforms/exponential_pressure_control.json` | Exponential control waveform. |
| `examples/waveforms/reusable_waveforms.json` | Named waveform references. |

See [Fields and waveforms](fields-waveforms.md) before changing waveform
parameters. Time-dependent fields require timestep convergence checks.

## Field arrays and ion clouds

Supporting data/examples:

- `examples/field_arrays/create_example_field_array.py`: creates sample HDF5
  field arrays.
- `examples/ion_clouds/default_cloud.json`: reusable ion cloud input.

Field arrays are external inputs. Archive them with the configuration and output
when a run must be reproducible.

## Parameter studies

Use CLI overrides for small controlled sweeps:

```bash
for model in NoCollisions HSS EHSS; do
  icarion \
    --set physics.collision_model="$model" \
    --output-dir "./results/$model" \
    examples/ims/ims_basic.json
done
```

For larger studies, copy the config into a run directory, record the exact
command line, and keep a fixed `simulation.rng_seed` unless seed variation is
part of the experiment.

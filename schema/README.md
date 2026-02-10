# ICARION Configuration Schema (v1.0)

JSON Schema definitions for validating ICARION configuration files.

## Structure

```
schema/
|-- icarion-config.schema.json    # Master schema (entry point)
|-- simulation.schema.json        # Simulation parameters
|-- physics.schema.json           # Collision models, reactions
|-- ions.schema.json              # Ion initialization
|-- domain.schema.json            # Domain definition
|-- geometry.schema.json          # Trap geometry
|-- environment.schema.json       # Gas conditions
|-- fields.schema.json            # Electric/magnetic fields
|-- waveform.schema.json          # Time-varying fields
|-- output.schema.json            # Output settings
|-- species.schema.json           # Species database format
|-- reactions.schema.json         # Reactions database format
|-- common-types.schema.json      # Reusable definitions
|-- boundary.schema.json          # Boundary conditions
`-- validator.py                  # Validation script
```

## Validation

Requires `jsonschema`:
```bash
pip install jsonschema
```

**Command-line:**
```bash
./schema/validator.py path/to/config.json    # Single file (from repo root)
./schema/validator.py examples/*.json        # Multiple files
./schema/validator.py --verbose config.json  # Detailed output
```

**Python API (from repo root):**
```python
from pathlib import Path
import sys

sys.path.append("schema")
from validator import IcarionConfigValidator

validator = IcarionConfigValidator(Path("schema"))
is_valid, errors = validator.validate_file(Path("config.json"))
```

## Waveform Support

All voltage/frequency parameters support static values or time-varying waveforms.
`environment.pressure_Pa` supports the same syntax:

**Static (backward compatible):**
```json
"voltage_V": 100.0
```

**Time-varying:**
```json
"voltage_V": {"type": "linear", "start": 0, "end": 500, "end_time_s": 0.001}
```

**Waveform types:** `constant`, `linear`, `quadratic`, `exponential`, `sinusoidal`, `pwm`, `pulsed`, `arbitrary`

See `docs/CONFIG_GUIDE.md` and `examples/waveforms/` for details.

## Key Features

- **Type safety**: Numeric types (`number`, `integer`), enums, physical units in descriptions
- **Reusable types**: `vec3`, `voltage`, `frequency`, `instrument_type`, `collision_model`, `gas_species`
- **Validation**: Required fields, numeric/array constraints, enum validation
- **Strict mode**: `additionalProperties: false` catches typos

Most schemas use JSON Schema Draft 2020-12; `ions.schema.json`,
`species.schema.json`, and `reactions.schema.json` still use Draft-07.

## References

- [JSON Schema Specification](https://json-schema.org/specification.html)
- [JSON Schema Validator (Python)](https://python-jsonschema.readthedocs.io/)

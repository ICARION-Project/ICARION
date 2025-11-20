# ICARION Configuration Schema v1.0

This directory contains the formal JSON Schema definitions for ICARION configuration files.

## Structure

```
v1.0/
├── icarion-config.schema.json    # Master schema (entry point)
├── simulation.schema.json        # Simulation parameters
├── physics.schema.json           # Physics configuration
├── output.schema.json            # Output settings
├── domain.schema.json            # Domain definition
├── geometry.schema.json          # Domain geometry
├── environment.schema.json       # Gas environment
├── fields.schema.json            # Electric/magnetic fields
└── common-types.schema.json      # Reusable type definitions
```

## Validation

### Prerequisites

```bash
pip install jsonschema
```

### Command-line Usage

```bash
# Validate single file
./validator.py path/to/config.json

# Validate multiple files
./validator.py examples/*.json

# Use custom schema directory
./validator.py --schema-dir /path/to/schemas config.json

# Verbose output
./validator.py --verbose config.json
```

### Python API

```python
from pathlib import Path
from validator import IcarionConfigValidator

# Initialize validator
validator = IcarionConfigValidator(Path("v1.0"))

# Validate file
is_valid, errors = validator.validate_file(Path("config.json"))

if not is_valid:
    for error in errors:
        print(f"Error: {error}")

# Validate dictionary
config = {...}
is_valid, errors = validator.validate_config(config)
```

## Schema Features

### Type Safety
- All numeric types use `number` or `integer` (JSON Schema compliant)
- Physical units documented in descriptions
- Enums for instrument types, collision models, solvers, gas species

### Reusable Types
Common types defined in `common-types.schema.json`:
- `vec3`: 3D vector [x, y, z]
- `voltage`, `frequency`, `time`, `length`, `pressure`, `temperature`
- `instrument_type`, `collision_model`, `solver_type`, `gas_species`

### Validation Rules
- Required fields enforced
- Numeric constraints (min/max/exclusive)
- Array length constraints
- Enum validation
- Cross-referenced types via `$ref`

### Strict Mode
All schemas use `"additionalProperties": false` to catch typos and unsupported fields.

## Migration from v0.x

Key changes:
- `double` → `number`
- `int` → `integer`
- `uint64` → `integer` with `minimum: 0`
- `Vec3` → `array` with `minItems/maxItems: 3`
- Removed trailing commas (JSON compliance)
- All instrument name aliases now in enum
- SIFDT_MS merged into IMS

## Schema Versioning

- **v1.0**: Current version (this directory)
- **Future versions**: Will be in separate `v1.1/`, `v1.2/` directories
- Master schema `$id` includes version for tracking

## Contributing

When modifying schemas:
1. Keep all schemas JSON Schema Draft 2020-12 compliant
2. Add descriptions for all properties
3. Include physical units in descriptions
4. Test with `validator.py` before committing
5. Update this README if adding new schemas

## References

- [JSON Schema Specification](https://json-schema.org/specification.html)
- [Understanding JSON Schema](https://json-schema.org/understanding-json-schema/)
- [JSON Schema Validator (Python)](https://python-jsonschema.readthedocs.io/)

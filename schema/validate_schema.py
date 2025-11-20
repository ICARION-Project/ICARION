#!/usr/bin/env python3
"""
Validate JSON files against JSON Schema
Usage: python3 validate_schema.py <schema_file> <json_file>
"""

import sys
import json
import argparse
from pathlib import Path

try:
    import jsonschema
    from jsonschema import validate, ValidationError
except ImportError:
    print("Error: jsonschema module not installed")
    print("Install with: pip3 install jsonschema")
    sys.exit(1)


def validate_json(schema_path: Path, json_path: Path) -> bool:
    """Validate a JSON file against a schema."""
    
    # Load schema
    try:
        with open(schema_path, 'r') as f:
            schema = json.load(f)
    except FileNotFoundError:
        print(f"❌ Error: Schema file not found: {schema_path}")
        return False
    except json.JSONDecodeError as e:
        print(f"❌ Error: Invalid JSON in schema file: {e}")
        return False
    
    # Load JSON data
    try:
        with open(json_path, 'r') as f:
            data = json.load(f)
    except FileNotFoundError:
        print(f"❌ Error: JSON file not found: {json_path}")
        return False
    except json.JSONDecodeError as e:
        print(f"❌ Error: Invalid JSON in data file: {e}")
        return False
    
    # Validate
    try:
        validate(instance=data, schema=schema)
        print(f"✅ Valid: {json_path.name} conforms to {schema_path.name}")
        return True
    except ValidationError as e:
        print(f"❌ Validation Error in {json_path.name}:")
        print(f"   Path: {' -> '.join(str(p) for p in e.path)}")
        print(f"   Message: {e.message}")
        if e.context:
            print(f"   Context:")
            for ctx_err in e.context:
                print(f"     - {ctx_err.message}")
        return False


def main():
    parser = argparse.ArgumentParser(
        description="Validate JSON files against JSON Schema"
    )
    parser.add_argument(
        "schema",
        type=Path,
        help="Path to JSON Schema file"
    )
    parser.add_argument(
        "json_file",
        type=Path,
        help="Path to JSON file to validate"
    )
    
    args = parser.parse_args()
    
    success = validate_json(args.schema, args.json_file)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()

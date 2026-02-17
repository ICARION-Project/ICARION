#!/usr/bin/env python3
"""
ICARION Configuration Validator
Validates JSON configuration files against the ICARION v1.0.0 schema.
"""

import json
import sys
from pathlib import Path
from typing import Dict, List, Optional

try:
    from jsonschema import Draft202012Validator, RefResolver, ValidationError
    from jsonschema.exceptions import SchemaError
except ImportError:
    print("ERROR: jsonschema library not found. Install with:")
    print("  pip install jsonschema")
    sys.exit(1)


class IcarionConfigValidator:
    """Validator for ICARION configuration files."""
    
    def __init__(self, schema_dir: Path):
        """
        Initialize validator with schema directory.
        
        Args:
            schema_dir: Path to directory containing schema files
        """
        self.schema_dir = schema_dir
        self.master_schema_path = schema_dir / "icarion-config.schema.json"
        
        if not self.master_schema_path.exists():
            raise FileNotFoundError(f"Master schema not found: {self.master_schema_path}")
        
        # Load all schemas into a store for local resolution
        schema_store = {}
        for schema_file in schema_dir.glob("*.schema.json"):
            with open(schema_file) as f:
                schema = json.load(f)
                if "$id" in schema:
                    schema_store[schema["$id"]] = schema
        
        # Load master schema
        with open(self.master_schema_path) as f:
            self.master_schema = json.load(f)
        
        # Create resolver with pre-loaded schemas
        schema_uri = self.master_schema_path.as_uri()
        self.resolver = RefResolver(
            base_uri=schema_uri,
            referrer=self.master_schema,
            store=schema_store
        )
        
        # Create validator
        self.validator = Draft202012Validator(self.master_schema, resolver=self.resolver)
    
    def validate_file(self, config_path: Path) -> tuple[bool, List[str]]:
        """
        Validate a configuration file.
        
        Args:
            config_path: Path to configuration JSON file
        
        Returns:
            Tuple of (is_valid, error_messages)
        """
        try:
            with open(config_path) as f:
                config = json.load(f)
        except json.JSONDecodeError as e:
            return False, [f"JSON parse error: {e}"]
        except FileNotFoundError:
            return False, [f"File not found: {config_path}"]
        
        return self.validate_config(config)
    
    def validate_config(self, config: Dict) -> tuple[bool, List[str]]:
        """
        Validate a configuration dictionary.
        
        Args:
            config: Configuration dictionary
        
        Returns:
            Tuple of (is_valid, error_messages)
        """
        errors = []
        config = self._strip_comments(config)

        if self._is_ion_cloud(config):
            return self._validate_ion_cloud(config)
        
        try:
            # Validate against schema
            for error in self.validator.iter_errors(config):
                # Format error message with path
                path = " → ".join(str(p) for p in error.absolute_path)
                if path:
                    errors.append(f"[{path}] {error.message}")
                else:
                    errors.append(error.message)
        
        except SchemaError as e:
            return False, [f"Schema error: {e}"]
        
        return len(errors) == 0, errors

    def _strip_comments(self, value):
        """
        Recursively remove comment-style fields (keys starting with "_")
        before validation.
        """
        if isinstance(value, dict):
            return {
                key: self._strip_comments(val)
                for key, val in value.items()
                if not key.startswith("_")
            }
        if isinstance(value, list):
            return [self._strip_comments(item) for item in value]
        return value

    def _is_ion_cloud(self, config: Dict) -> bool:
        """
        Detect standalone ion cloud files (not full simulation configs).
        """
        if not isinstance(config, dict):
            return False
        if "ions" not in config or not isinstance(config.get("ions"), list):
            return False
        has_full_sections = any(
            key in config for key in ("simulation", "physics", "output", "domains")
        )
        return not has_full_sections

    def _validate_ion_cloud(self, config: Dict) -> tuple[bool, List[str]]:
        """
        Validate ion cloud file structure.
        """
        errors = []
        ions = config.get("ions", [])
        if not isinstance(ions, list):
            return False, ["[ions] must be an array"]

        for idx, ion in enumerate(ions):
            path_prefix = f"ions[{idx}]"
            if not isinstance(ion, dict):
                errors.append(f"[{path_prefix}] must be an object")
                continue
            species = ion.get("species")
            if not isinstance(species, str) or not species:
                errors.append(f"[{path_prefix} → species] must be a non-empty string")
            for vec_key in ("pos", "vel"):
                if vec_key in ion:
                    vec = ion.get(vec_key)
                    if (
                        not isinstance(vec, list)
                        or len(vec) < 3
                        or not all(isinstance(v, (int, float)) for v in vec[:3])
                    ):
                        errors.append(f"[{path_prefix} → {vec_key}] must be a 3-element numeric array")
            if "birth_time" in ion:
                if not isinstance(ion.get("birth_time"), (int, float)):
                    errors.append(f"[{path_prefix} → birth_time] must be a number")

        return len(errors) == 0, errors


def main():
    """Command-line interface for validation."""
    import argparse
    
    parser = argparse.ArgumentParser(
        description="Validate ICARION configuration files",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate single file
  %(prog)s config.json
  
  # Validate multiple files
  %(prog)s config1.json config2.json
  
  # Validate all examples
  %(prog)s examples/*.json
  
  # Use custom schema directory
  %(prog)s --schema-dir /path/to/schemas config.json
"""
    )
    
    parser.add_argument(
        "config_files",
        nargs="+",
        type=Path,
        help="Configuration file(s) to validate"
    )
    
    parser.add_argument(
        "--schema-dir",
        type=Path,
        default=Path(__file__).parent,
        help="Directory containing schema files (default: schema/)"
    )
    
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed error messages"
    )
    
    args = parser.parse_args()
    
    # Initialize validator
    try:
        validator = IcarionConfigValidator(args.schema_dir)
    except Exception as e:
        print(f"ERROR: Failed to initialize validator: {e}", file=sys.stderr)
        return 1
    
    # Validate each file
    total_files = len(args.config_files)
    valid_files = 0
    
    for config_file in args.config_files:
        print(f"\n{'='*60}")
        print(f"Validating: {config_file}")
        print('='*60)
        
        is_valid, errors = validator.validate_file(config_file)
        
        if is_valid:
            print("✓ VALID")
            valid_files += 1
        else:
            print("✗ INVALID")
            print(f"\nFound {len(errors)} error(s):")
            for i, error in enumerate(errors, 1):
                print(f"  {i}. {error}")
    
    # Summary
    print(f"\n{'='*60}")
    print(f"Summary: {valid_files}/{total_files} files valid")
    print('='*60)
    
    return 0 if valid_files == total_files else 1


if __name__ == "__main__":
    sys.exit(main())

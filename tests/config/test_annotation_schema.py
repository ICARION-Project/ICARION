#!/usr/bin/env python3
"""Exercise annotation fields through ICARION's public schema validator."""

import copy
import json
import sys
from pathlib import Path


def main() -> int:
    repo_root = Path(sys.argv[1]).resolve()
    sys.path.insert(0, str(repo_root / "schema"))
    from validator import IcarionConfigValidator

    validator = IcarionConfigValidator(repo_root / "schema")
    with (repo_root / "examples/ims/ims_output_diagnostics_minimal.json").open(encoding="utf-8") as stream:
        base = json.load(stream)

    cases = [
        ({"note": "schema annotation"}, True, "metadata.note"),
        ({"note_file": "note.md"}, True, "metadata.note_file"),
        ({"note": "one", "note_file": "note.md"}, False, "both annotation fields"),
        ({"unknown": True}, False, "unknown metadata property"),
        ({"note": ""}, False, "empty inline note"),
    ]
    for metadata, expected, label in cases:
        config = copy.deepcopy(base)
        config["metadata"] = metadata
        valid, errors = validator.validate_config(config)
        if valid != expected:
            raise AssertionError(f"{label}: expected valid={expected}, got {valid}; errors={errors}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

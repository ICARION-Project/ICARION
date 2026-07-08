#!/usr/bin/env python3
"""Summarize TIMS elution validation artifacts."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def find_log(results_dir: Path) -> Path:
    candidates = [
        results_dir / "TIMS_ELUTION_VALIDATION.txt",
        results_dir.parent / "logs" / "TIMS_ELUTION_VALIDATION.txt",
        results_dir.parents[1] / "logs" / "TIMS_ELUTION_VALIDATION.txt" if len(results_dir.parents) > 1 else results_dir,
        Path("validation/logs/TIMS_ELUTION_VALIDATION.txt"),
    ]
    for path in candidates:
        if path.exists():
            return path
    raise FileNotFoundError("TIMS_ELUTION_VALIDATION.txt")


def analyze_tims_results(results_dir: Path) -> dict:
    results_dir = Path(results_dir)
    log_path = find_log(results_dir)
    text = log_path.read_text(encoding="utf-8", errors="replace")
    status = "PASS" if "RESULT:" in text and "PASS" in text.split("RESULT:", 1)[1].splitlines()[0] else "FAIL"
    h5_files = sorted(str(path) for path in results_dir.glob("*.h5"))
    return {
        "validation": "TIMS elution",
        "status": status,
        "artifacts": {
            "results_dir": str(results_dir),
            "log": str(log_path),
            "h5_files": h5_files,
        },
    }


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: analyze_tims_elution.py <results_dir>", file=sys.stderr)
        return 2
    try:
        result = analyze_tims_results(Path(sys.argv[1]))
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    print(json.dumps(result, indent=2))
    return 0 if result["status"] == "PASS" else 1


if __name__ == "__main__":
    raise SystemExit(main())

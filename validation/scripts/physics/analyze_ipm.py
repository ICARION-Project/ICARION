#!/usr/bin/env python3
"""Summarize InteractionPotentialModel validation gate results."""

from __future__ import annotations

import json
import sys
from pathlib import Path


def analyze_ipm_results(results_dir: Path) -> dict:
    results_dir = Path(results_dir)
    summary_path = results_dir / "ipm_gate_summary.json"
    if not summary_path.exists():
        raise FileNotFoundError(summary_path)

    with open(summary_path, "r", encoding="utf-8") as handle:
        summary = json.load(handle)

    gates = summary.get("gates", [])
    statuses = [gate.get("status", "UNKNOWN") for gate in gates]
    if "FAIL" in statuses:
        status = "FAIL"
    elif "SKIPPED" in statuses:
        status = "PARTIAL"
    elif statuses and all(item == "PASS" for item in statuses):
        status = "PASS"
    else:
        status = "UNKNOWN"

    return {
        "validation": "InteractionPotentialModel",
        "status": status,
        "gates": gates,
        "artifacts": summary.get("artifacts", {}),
        "summary": {
            "total_gates": len(gates),
            "passed": sum(1 for gate in gates if gate.get("status") == "PASS"),
            "failed": sum(1 for gate in gates if gate.get("status") == "FAIL"),
            "skipped": sum(1 for gate in gates if gate.get("status") == "SKIPPED"),
        },
    }


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: analyze_ipm.py <results_dir>", file=sys.stderr)
        return 2

    try:
        result = analyze_ipm_results(Path(sys.argv[1]))
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1

    print(json.dumps(result, indent=2))
    return 0 if result["status"] in ("PASS", "PARTIAL") else 1


if __name__ == "__main__":
    raise SystemExit(main())

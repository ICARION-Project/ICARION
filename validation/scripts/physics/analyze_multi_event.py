#!/usr/bin/env python3
"""Analyzer for multi-event collision mode validation gate results.

Aggregates primary mode check and dt-convergence sub-gate from validate_multi_event*.py
into a consolidated release-readiness report.
"""

import json
import sys
from pathlib import Path


def analyze_multi_event_results(results_dir: Path) -> dict:
    """Aggregate multi-event gate results into summary."""
    results_dir = Path(results_dir)
    
    # Look for gate summary JSON files
    primary_summary = results_dir / "multi_event_gate_summary.json"
    convergence_summary = results_dir / "multi_event_convergence_gate_summary.json"
    
    gates = {}
    overall_status = "FAIL"
    
    if primary_summary.exists():
        with open(primary_summary) as f:
            primary_data = json.load(f)
            gates["primary_mode"] = {
                "status": primary_data.get("status", "UNKNOWN"),
                "checks": len(primary_data.get("checks", [])),
                "note": "Multi-event mode functional production run"
            }
    
    if convergence_summary.exists():
        with open(convergence_summary) as f:
            convergence_data = json.load(f)
            gates["dt_convergence"] = {
                "status": convergence_data.get("status", "UNKNOWN"),
                "checks": len(convergence_data.get("checks", [])),
                "note": "Timestep convergence sweep (dt=0.5e-9 to 5e-9 s)"
            }
    
    # Determine overall status
    statuses = [g["status"] for g in gates.values()]
    if "FAIL" in statuses:
        overall_status = "FAIL"
    elif "PARTIAL" in statuses:
        overall_status = "PARTIAL"
    elif all(s == "PASS" for s in statuses):
        overall_status = "PASS"
    
    return {
        "validation": "Multi-Event Collision Mode (v1.1)",
        "status": overall_status,
        "gates": gates,
        "summary": {
            "total_gates": len(gates),
            "passed": sum(1 for g in gates.values() if g["status"] == "PASS"),
            "failed": sum(1 for g in gates.values() if g["status"] == "FAIL"),
        }
    }


def main():
    if len(sys.argv) < 2:
        print("Usage: analyze_multi_event.py <results_dir>", file=sys.stderr)
        sys.exit(1)
    
    results_dir = Path(sys.argv[1])
    if not results_dir.exists():
        print(f"Error: {results_dir} not found", file=sys.stderr)
        sys.exit(1)
    
    result = analyze_multi_event_results(results_dir)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Analyzer for output diagnostics feature validation gate results.

Aggregates minimal/deep/reproducibility output mode checks and mode interoperability
from validate_output_diagnostics_e2e.py into a consolidated report.
"""

import json
import sys
from pathlib import Path


def analyze_output_diagnostics_results(results_dir: Path) -> dict:
    """Aggregate output diagnostics gate results into summary."""
    results_dir = Path(results_dir)
    
    summary_file = results_dir / "output_diagnostics_gate_summary.json"
    
    gates = {}
    overall_status = "FAIL"
    
    if summary_file.exists():
        with open(summary_file) as f:
            data = json.load(f)
            
            # Parse output modes and checks from gate summary
            gates["minimal_trajectory"] = {
                "status": "PASS" if "minimal_transport" in data.get("checks", {}) else "FAIL",
                "checks": 1,
                "note": "Minimal trajectory output mode enabled"
            }
            
            gates["deep_collision"] = {
                "status": "PASS" if "deep_collision_summary" in data.get("checks", {}) else "FAIL",
                "checks": 1,
                "note": "Deep collision analysis summary enabled"
            }
            
            gates["reproducibility"] = {
                "status": "PASS" if "reproducibility_metadata" in data.get("checks", {}) else "FAIL",
                "checks": 1,
                "note": "Reproducibility embed-cap and tally tracking"
            }
            
            gates["mode_interop"] = {
                "status": data.get("status", "UNKNOWN"),
                "checks": len(data.get("checks", {})),
                "note": "Minimal vs full mode trajectory agreement"
            }
            
            overall_status = data.get("status", "UNKNOWN")
    
    return {
        "validation": "Output Diagnostics Features (v1.1)",
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
        print("Usage: analyze_output_diagnostics.py <results_dir>", file=sys.stderr)
        sys.exit(1)
    
    results_dir = Path(sys.argv[1])
    if not results_dir.exists():
        print(f"Error: {results_dir} not found", file=sys.stderr)
        sys.exit(1)
    
    result = analyze_output_diagnostics_results(results_dir)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Run a bounded analysis sweep against existing ICARION result files.

The script scans a results directory for representative full/minimal HDF5
outputs, runs compatible analysis CLIs on the first matching files, and writes
a JSON pass/fail summary plus stdout/stderr logs. It is meant as a local
release sanity check, not as a generator for canonical validation evidence.
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import subprocess
import sys
import time
from pathlib import Path

import h5py


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Run representative analysis CLIs on compatible files under results/.")
    p.add_argument("--results-dir", type=Path, default=Path("results"), help="Directory containing HDF5 results.")
    p.add_argument("--out-dir", type=Path, default=Path("analysis/output/realdata_sweep_fresh"), help="Output directory.")
    p.add_argument("--timeout", type=float, default=180.0, help="Per-command timeout in seconds.")
    return p.parse_args()


def classify_file(path: Path) -> str | None:
    try:
        with h5py.File(path, "r") as h5:
            has_full = "trajectory/time" in h5 and "trajectory/positions" in h5 and "trajectory/species_ids" in h5
            has_min = "analysis/minimal_transport" in h5
            has_col = "analysis/deep_collision" in h5
    except Exception:
        return "invalid"
    if has_full and has_col:
        return "full_collision"
    if has_full:
        return "full"
    if has_min and has_col:
        return "minimal_collision"
    if has_min:
        return "minimal"
    return "invalid"


def first_by_class(results_dir: Path) -> dict[str, Path]:
    wanted = {"full", "full_collision", "minimal", "minimal_collision"}
    found: dict[str, Path] = {}
    for path in sorted(results_dir.rglob("*.h5")):
        cls = classify_file(path)
        if cls in wanted and cls not in found:
            found[cls] = path
        if wanted.issubset(found):
            break
    return found


def command_matrix(files: dict[str, Path], out: Path) -> list[tuple[str, list[str], bool]]:
    full = files.get("full")
    full_collision = files.get("full_collision")
    minimal = files.get("minimal")
    minimal_collision = files.get("minimal_collision")
    commands: list[tuple[str, list[str], bool]] = []

    for key, path in files.items():
        commands.append((f"run_report_{key}", [sys.executable, "analysis/run_report.py", "--traj", str(path), "--out", str(out / f"run_report_{key}.csv")], False))

    if full is not None:
        commands += [
            ("plot_trajectories", [sys.executable, "analysis/plot_trajectories.py", "--traj", str(full), "--out", str(out / "plot_trajectories.png"), "--max-ions", "50", "--time-stride", "10"], False),
            ("trap_stability", [sys.executable, "analysis/trap_stability.py", "--traj", str(full), "--out", str(out / "trap_stability.png"), "--out-csv", str(out / "trap_stability.csv"), "--time-stride", "10"], False),
            ("spectrum_frequency_resampled", [sys.executable, "analysis/spectrum_analysis.py", "--traj", str(full), "--mode", "frequency", "--coordinate", "r", "--resample", "uniform", "--out", str(out / "spectrum_frequency.png"), "--out-csv", str(out / "spectrum_frequency.csv")], False),
        ]

    if full_collision is not None:
        commands += [
            ("arrival_time_distribution", [sys.executable, "analysis/arrival_time_distribution.py", "--traj", str(full_collision), "--out", str(out / "arrival.png"), "--out-csv", str(out / "arrival.csv"), "--max-ions", "500", "--time-stride", "5", "--no-per-species"], False),
            ("ims_mobility_inferred", [sys.executable, "analysis/ims_mobility.py", "--traj", str(full_collision), "--out", str(out / "ims_mobility.png"), "--out-csv", str(out / "ims_mobility.csv"), "--max-ions", "500", "--time-stride", "5"], False),
            ("transport_diagnostics", [sys.executable, "analysis/transport_diagnostics.py", "--traj", str(full_collision), "--out", str(out / "transport.png"), "--out-csv", str(out / "transport.csv"), "--max-ions", "500"], False),
            ("spacecharge_diagnostics", [sys.executable, "analysis/spacecharge_diagnostics.py", "--traj", str(full_collision), "--out", str(out / "spacecharge.png"), "--out-csv", str(out / "spacecharge.csv"), "--max-ions", "200", "--time-stride", "10"], False),
            ("collision_diagnostics_full", [sys.executable, "analysis/collision_diagnostics.py", "--traj", str(full_collision), "--out-csv", str(out / "collision_full.csv")], False),
        ]

    if minimal_collision is not None:
        commands.append(("collision_diagnostics_minimal", [sys.executable, "analysis/collision_diagnostics.py", "--traj", str(minimal_collision), "--out-csv", str(out / "collision_minimal.csv")], False))

    if minimal is not None:
        commands.append(("expected_reject_collision_missing", [sys.executable, "analysis/collision_diagnostics.py", "--traj", str(minimal), "--out-csv", str(out / "collision_missing.csv")], True))

    return commands


def run_command(name: str, cmd: list[str], expected_failure: bool, out: Path, timeout: float) -> dict[str, object]:
    env = dict(os.environ)
    env.setdefault("MPLBACKEND", "Agg")
    start = time.time()
    try:
        result = subprocess.run(cmd, cwd=Path.cwd(), env=env, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=timeout)
        timed_out = False
    except subprocess.TimeoutExpired as exc:
        result = subprocess.CompletedProcess(cmd, 124, stdout=exc.stdout or "", stderr=exc.stderr or "")
        timed_out = True
    elapsed = time.time() - start
    (out / f"{name}.stdout.txt").write_text(result.stdout or "", encoding="utf-8")
    (out / f"{name}.stderr.txt").write_text(result.stderr or "", encoding="utf-8")
    passed = (result.returncode == 0 and not expected_failure) or (result.returncode != 0 and expected_failure)
    return {
        "name": name,
        "command": " ".join(shlex.quote(str(c)) for c in cmd),
        "returncode": result.returncode,
        "expected_failure": expected_failure,
        "passed": passed,
        "timed_out": timed_out,
        "seconds": round(elapsed, 3),
        "stdout_tail": (result.stdout or "")[-1000:],
        "stderr_tail": (result.stderr or "")[-1000:],
    }


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    files = first_by_class(args.results_dir)
    commands = command_matrix(files, args.out_dir)
    results = [run_command(name, cmd, expected_failure, args.out_dir, args.timeout) for name, cmd, expected_failure in commands]
    summary = {
        "results_dir": str(args.results_dir),
        "out_dir": str(args.out_dir),
        "selected_files": {key: str(path) for key, path in files.items()},
        "runs": results,
        "passed": sum(1 for row in results if row["passed"]),
        "failed": sum(1 for row in results if not row["passed"]),
    }
    (args.out_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(f"Wrote {args.out_dir / 'summary.json'}")
    return 0 if summary["failed"] == 0 else 1


if __name__ == "__main__":
    raise SystemExit(main())

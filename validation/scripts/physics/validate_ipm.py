#!/usr/bin/env python3
"""InteractionPotentialModel validation gate runner.

This runner validates the current IPM path with freshly generated
`ipm_offline_samples` data:

1. Gate A: InteractionPotentialCollisionHandler unit/loader safety tests.
2. Gate B: build a mini HDF5 sample table with interaction_potential_precompute.
3. Gate C: run a short IMS-style simulation using the generated sample table.
4. Gate D: check finite trajectory and momentum diagnostics.

The generated sample is intentionally small. It is release validation evidence
for format/wiring/runtime consistency, not a high-fidelity mobility benchmark.
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
VALIDATION_DIR = REPO_ROOT / "validation"
DATA_SPECIES_DB = REPO_ROOT / "data" / "species_database_v1.json"
EXAMPLE_CONFIG = REPO_ROOT / "examples" / "ims" / "ims_ipm_basic.json"
DEFAULT_BIN = REPO_ROOT / "build" / "src" / "icarion_main"
DEFAULT_PRECOMPUTE_BIN = REPO_ROOT / "build" / "src" / "interaction_potential_precompute"
DEFAULT_BUILD_DIR = REPO_ROOT / "build"

RUN_DIR_ENV = os.environ.get("ICARION_VALIDATION_RUN_DIR")
if RUN_DIR_ENV:
    RUN_DIR = Path(RUN_DIR_ENV)
    OUT_DIR = RUN_DIR / "results" / "physics" / "ipm"
    LOG_DIR = RUN_DIR / "logs" / "ipm"
else:
    OUT_DIR = VALIDATION_DIR / "results" / "physics" / "ipm"
    LOG_DIR = VALIDATION_DIR / "logs" / "ipm"


@dataclass
class GateResult:
    gate: str
    status: str
    note: str


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run InteractionPotentialModel validation gates")
    parser.add_argument("--icarion-bin", type=Path, default=DEFAULT_BIN)
    parser.add_argument("--precompute-bin", type=Path, default=DEFAULT_PRECOMPUTE_BIN)
    parser.add_argument("--build-dir", type=Path, default=DEFAULT_BUILD_DIR)
    parser.add_argument("--species", default="H3O+")
    parser.add_argument("--gas", default="He")
    parser.add_argument("--pressure-pa", type=float, default=80.0)
    parser.add_argument("--temperature-k", type=float, default=300.0)
    parser.add_argument("--ion-count", type=int, default=120)
    parser.add_argument("--dt-s", type=float, default=2.0e-9)
    parser.add_argument("--total-time-s", type=float, default=5.0e-6)
    parser.add_argument("--write-interval", type=int, default=20)
    parser.add_argument("--n-orientations", type=int, default=2)
    parser.add_argument("--n-trials", type=int, default=16)
    parser.add_argument("--v-bins", type=int, default=5)
    parser.add_argument("--v-min", type=float, default=200.0)
    parser.add_argument("--v-max", type=float, default=20000.0)
    parser.add_argument("--threads", type=int, default=1)
    parser.add_argument("--repeat-ctest", type=int, default=1)
    parser.add_argument("--timeout", type=int, default=900)
    parser.add_argument("--skip-ctest", action="store_true")
    parser.add_argument("--skip-sim", action="store_true")
    return parser.parse_args()


def run_cmd(cmd: list[str], cwd: Path, log_path: Path, timeout: int) -> tuple[int, float, str]:
    start = time.perf_counter()
    proc = subprocess.run(
        cmd,
        cwd=str(cwd),
        capture_output=True,
        text=True,
        timeout=timeout,
    )
    elapsed = time.perf_counter() - start
    output = (proc.stdout or "") + (("\n" + proc.stderr) if proc.stderr else "")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(output, encoding="utf-8")
    return proc.returncode, elapsed, output


def copy_data_tree_for_precompute(work_dir: Path) -> Path:
    """Create a small self-contained input tree preserving geometry_file paths."""

    input_dir = work_dir / "precompute_input"
    molecules_dir = input_dir / "molecules"
    molecules_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(DATA_SPECIES_DB, input_dir / "species_database_v1.json")
    shutil.copy2(REPO_ROOT / "data" / "molecules" / "H3O+.json", molecules_dir / "H3O+.json")
    return input_dir / "species_database_v1.json"


def write_validation_forcefield(path: Path) -> None:
    """Write the minimal LJ parameter set needed for H3O+/He validation."""

    params = {
        "gases": {
            "He": {
                "sigma_A": 2.551,
                "epsilon_eV": 0.0008798,
                "alpha_A3": 0.204956,
            }
        },
        "elements": {
            "H": {
                "sigma_A": 2.261,
                "epsilon_eV": 0.059579817756640736,
            },
            "O": {
                "sigma_A": 2.4344,
                "epsilon_eV": 0.1034324933362249,
            },
        },
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(params, handle, indent=2)


def inspect_sample_file(path: Path, expected_species: str, expected_gas: str) -> tuple[bool, str]:
    if not path.exists():
        return False, "sample file was not created"

    try:
        with h5py.File(path, "r") as handle:
            fmt = handle.attrs.get("format", "")
            species = handle.attrs.get("species_id", "")
            gas = handle.attrs.get("gas", "")
            if isinstance(fmt, bytes):
                fmt = fmt.decode()
            if isinstance(species, bytes):
                species = species.decode()
            if isinstance(gas, bytes):
                gas = gas.decode()

            if fmt != "ipm_offline_samples":
                return False, f"unexpected format attribute: {fmt!r}"
            if species != expected_species or gas != expected_gas:
                return False, f"metadata mismatch species={species!r} gas={gas!r}"

            required = ["logv_bins", "orientations_quat", "sigma_mt_m2", "b_max_m", "dp_stats"]
            missing = [name for name in required if name not in handle]
            if missing:
                return False, f"missing datasets: {', '.join(missing)}"

            logv = np.asarray(handle["logv_bins"][:], dtype=float)
            sigma = np.asarray(handle["sigma_mt_m2"][:], dtype=float)
            b_max = np.asarray(handle["b_max_m"][:], dtype=float)
            dp_stats = np.asarray(handle["dp_stats"][:], dtype=float)

            if logv.size == 0 or sigma.size == 0 or b_max.size == 0 or dp_stats.size == 0:
                return False, "one or more required datasets are empty"
            if not all(np.isfinite(arr).all() for arr in (logv, sigma, b_max, dp_stats)):
                return False, "sample datasets contain NaN/Inf"
            if not np.any(sigma > 0.0):
                return False, "sigma_mt_m2 has no positive entries"
            if not np.all(b_max > 0.0):
                return False, "b_max_m contains non-positive entries"

            return True, (
                f"format ok; orientations={handle.attrs.get('n_orientations')}; "
                f"v_bins={logv.size}; sigma_range=[{sigma.min():.3e}, {sigma.max():.3e}] m^2"
            )
    except Exception as exc:
        return False, f"failed to inspect sample file: {exc}"


def write_runtime_species_db(path: Path, sample_path: Path, species: str) -> None:
    # Use the full species database so the global reaction-database fallback can
    # validate species names even though reactions are disabled for this run.
    with open(DATA_SPECIES_DB, "r", encoding="utf-8") as handle:
        species_db = json.load(handle)
    species_db["species"][species]["ipm_samples_file"] = str(sample_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as handle:
        json.dump(species_db, handle, indent=2)


def write_runtime_config(args: argparse.Namespace, path: Path, out_dir: Path, sample_path: Path) -> tuple[Path, Path]:
    with open(EXAMPLE_CONFIG, "r", encoding="utf-8") as handle:
        config = json.load(handle)

    runtime_species_db = out_dir / "species_database_ipm_validation.json"
    write_runtime_species_db(runtime_species_db, sample_path, args.species)

    vrel_csv = out_dir / "ipm_vrel.csv"
    momentum_csv = out_dir / "ipm_momentum.csv"

    config["simulation"]["total_time_s"] = args.total_time_s
    config["simulation"]["dt_s"] = args.dt_s
    config["simulation"]["write_interval"] = args.write_interval
    config["simulation"]["rng_seed"] = 123
    config["physics"]["ipm_vrel_log_prefix"] = str(vrel_csv)
    config["physics"]["ipm_momentum_log_prefix"] = str(momentum_csv)
    config["species_database"] = str(runtime_species_db)
    config["output"]["folder"] = str(out_dir)
    config["output"]["trajectory_file"] = "ipm_validation.h5"
    config["output"]["trajectory_mode"] = "minimal"
    config["output"]["print_progress"] = False
    config["ions"]["species"][0]["id"] = args.species
    config["ions"]["species"][0]["count"] = args.ion_count
    config["domains"][0]["env"]["pressure_Pa"] = args.pressure_pa
    config["domains"][0]["env"]["temperature_K"] = args.temperature_k
    config["domains"][0]["env"]["gas_species"] = args.gas

    with open(path, "w", encoding="utf-8") as handle:
        json.dump(config, handle, indent=2)
    return vrel_csv, momentum_csv


def analyze_trajectory(h5_path: Path) -> tuple[bool, str]:
    if not h5_path.exists():
        return False, "trajectory file missing"

    with h5py.File(h5_path, "r") as handle:
        minimal_root = "/analysis/minimal_transport"
        required = [
            f"{minimal_root}/final_pos_x",
            f"{minimal_root}/final_pos_y",
            f"{minimal_root}/final_pos_z",
            f"{minimal_root}/final_vel_x",
            f"{minimal_root}/final_vel_y",
            f"{minimal_root}/final_vel_z",
        ]
        if all(path in handle for path in required):
            arrays = [np.asarray(handle[path][:], dtype=float) for path in required]
            if any(arr.size == 0 for arr in arrays):
                return False, "minimal_transport contains empty final arrays"
            if not all(np.isfinite(arr).all() for arr in arrays):
                return False, "minimal_transport final arrays contain NaN/Inf"
            return True, "finite minimal_transport final position/velocity arrays"

        if "/trajectory/positions" in handle:
            pos = np.asarray(handle["/trajectory/positions"][:], dtype=float)
            if pos.size == 0:
                return False, "empty position trajectory dataset"
            if not np.isfinite(pos).all():
                return False, "position dataset contains NaN/Inf"
            return True, "finite position trajectory dataset"

    return False, "no supported trajectory datasets found"


def analyze_momentum_csv(csv_path: Path, sample_path: Path) -> tuple[bool, str]:
    if not csv_path.exists():
        return False, "momentum CSV not created"

    rows = [line.strip() for line in csv_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    if len(rows) < 2:
        return False, "momentum CSV has no data row"

    headers = rows[0].split(",")
    values = rows[1].split(",")
    if len(headers) != len(values):
        return False, "momentum CSV row/header mismatch"
    data = dict(zip(headers, values))

    required = ["attempts", "collisions", "mean_k_s-1", "mean_vrel_m_s", "mean_dp_par_SI"]
    parsed: dict[str, float] = {}
    for key in required:
        try:
            parsed[key] = float(data[key])
        except Exception:
            return False, f"momentum CSV parse failed for {key}"

    if not all(np.isfinite(value) for value in parsed.values()):
        return False, "momentum CSV contains NaN/Inf"
    if parsed["attempts"] <= 0.0 or parsed["collisions"] <= 0.0:
        return False, "momentum CSV shows no attempts/collisions"

    try:
        with h5py.File(sample_path, "r") as handle:
            sample_v = np.exp(np.asarray(handle["logv_bins"][:], dtype=float))
            sample_v_min = float(np.min(sample_v))
            sample_v_max = float(np.max(sample_v))
    except Exception as exc:
        return False, f"failed to read sample velocity coverage: {exc}"

    mean_vrel = parsed["mean_vrel_m_s"]
    if mean_vrel < sample_v_min or mean_vrel > sample_v_max:
        return False, (
            f"mean_vrel={mean_vrel:.1f} m/s outside sample coverage "
            f"[{sample_v_min:.1f}, {sample_v_max:.1f}] m/s"
        )

    return True, (
        f"finite momentum observables; attempts={parsed['attempts']:.0f}; "
        f"collisions={parsed['collisions']:.0f}; mean_vrel={mean_vrel:.1f} m/s; "
        f"sample_v=[{sample_v_min:.1f}, {sample_v_max:.1f}] m/s"
    )


def write_summary(summary_path: Path, md_path: Path, payload: dict[str, Any]) -> None:
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    with open(summary_path, "w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)

    lines = [
        "# InteractionPotentialModel Validation Summary",
        "",
        f"Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "## Gates",
    ]
    for gate in payload["gates"]:
        lines.append(f"- {gate['gate']}: {gate['status']} - {gate['note']}")

    lines += [
        "",
        "## Artifacts",
    ]
    for key, value in payload["artifacts"].items():
        lines.append(f"- {key}: {value}")

    lines += [
        "",
        "## Scope",
        "- Fresh `ipm_offline_samples` HDF5 data is generated during the run.",
        "- This validates IPM sample format, config wiring, runtime loading, and finite observables.",
        "- It is not a high-fidelity mobility benchmark; use larger precompute settings for release evidence.",
    ]
    md_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    gates: list[GateResult] = []
    artifacts: dict[str, str] = {}

    if not args.precompute_bin.exists():
        print(f"ERROR: precompute binary not found: {args.precompute_bin}", file=sys.stderr)
        return 2
    if not args.icarion_bin.exists() and not args.skip_sim:
        print(f"ERROR: icarion binary not found: {args.icarion_bin}", file=sys.stderr)
        return 2

    if args.skip_ctest:
        gates.append(GateResult("Gate A", "SKIPPED", "skipped via --skip-ctest"))
    else:
        repeat_ok = True
        repeat_notes: list[str] = []
        for idx in range(args.repeat_ctest):
            log_path = LOG_DIR / f"gateA_ctest_ipm_run{idx + 1}.log"
            rc, elapsed, _ = run_cmd(
                ["ctest", "--output-on-failure", "-R", "InteractionPotentialCollisionHandler"],
                cwd=args.build_dir,
                log_path=log_path,
                timeout=args.timeout,
            )
            artifacts[f"gateA_ctest_run{idx + 1}_log"] = str(log_path)
            repeat_notes.append(f"run {idx + 1}: rc={rc} elapsed_s={elapsed:.2f}")
            if rc != 0:
                repeat_ok = False
        status = "PASS" if repeat_ok else "FAIL"
        gates.append(GateResult("Gate A", status, "; ".join(repeat_notes)))

    input_db = copy_data_tree_for_precompute(OUT_DIR)
    forcefield_path = OUT_DIR / "precompute_input" / "ipm_validation_forcefield.json"
    write_validation_forcefield(forcefield_path)
    sample_path = OUT_DIR / "precomputed_ipm" / f"{args.species}_ipm_samples_{args.gas}.h5"
    sample_path.parent.mkdir(parents=True, exist_ok=True)
    precompute_log = LOG_DIR / "gateB_precompute_ipm.log"
    precompute_cmd = [
        str(args.precompute_bin),
        "--input", str(input_db),
        "--species", args.species,
        "--gas", args.gas,
        "--output", str(sample_path),
        "--gas-params", str(forcefield_path),
        "--element-params", str(forcefield_path),
        "--orient-grid", "random",
        "--n-orientations", str(args.n_orientations),
        "--n-trials", str(args.n_trials),
        "--v-bins", str(args.v_bins),
        "--v-min", str(args.v_min),
        "--v-max", str(args.v_max),
        "--threads", str(args.threads),
        "--seed", "12345",
    ]
    rc_pre, elapsed_pre, _ = run_cmd(precompute_cmd, REPO_ROOT, precompute_log, args.timeout)
    artifacts["precompute_log"] = str(precompute_log)
    artifacts["sample_file"] = str(sample_path)

    sample_ok = False
    sample_note = f"precompute failed rc={rc_pre}"
    if rc_pre == 0:
        sample_ok, sample_note = inspect_sample_file(sample_path, args.species, args.gas)
    gates.append(GateResult("Gate B", "PASS" if sample_ok else "FAIL", f"{sample_note}; precompute_s={elapsed_pre:.2f}"))

    sim_ok = True
    sim_note = "skipped via --skip-sim"
    traj_path = OUT_DIR / "ipm_validation.h5"
    runtime_config = OUT_DIR / "ipm_validation_config.json"
    vrel_csv = OUT_DIR / "ipm_vrel.csv"
    momentum_csv = OUT_DIR / "ipm_momentum.csv"
    artifacts["runtime_config"] = str(runtime_config)
    artifacts["trajectory"] = str(traj_path)
    artifacts["vrel_csv"] = str(vrel_csv)
    artifacts["momentum_csv"] = str(momentum_csv)

    if args.skip_sim:
        gates.append(GateResult("Gate C", "SKIPPED", sim_note))
        gates.append(GateResult("Gate D", "SKIPPED", sim_note))
    elif not sample_ok:
        gates.append(GateResult("Gate C", "FAIL", "simulation skipped because sample gate failed"))
        gates.append(GateResult("Gate D", "FAIL", "observable checks skipped because sample gate failed"))
        sim_ok = False
    else:
        vrel_csv, momentum_csv = write_runtime_config(args, runtime_config, OUT_DIR, sample_path)
        sim_log = LOG_DIR / "gateC_simulation.log"
        rc_sim, elapsed_sim, _ = run_cmd(
            [str(args.icarion_bin), str(runtime_config), "--threads", "1"],
            cwd=REPO_ROOT,
            log_path=sim_log,
            timeout=args.timeout,
        )
        artifacts["simulation_log"] = str(sim_log)
        if rc_sim != 0:
            sim_ok = False
            sim_note = f"simulation failed rc={rc_sim}"
        else:
            sim_note = f"simulation completed in {elapsed_sim:.2f}s"
        gates.append(GateResult("Gate C", "PASS" if sim_ok else "FAIL", sim_note))

        if sim_ok:
            traj_ok, traj_note = analyze_trajectory(traj_path)
            mom_ok, mom_note = analyze_momentum_csv(momentum_csv, sample_path)
            obs_ok = traj_ok and mom_ok
            gates.append(GateResult("Gate D", "PASS" if obs_ok else "FAIL", f"traj={traj_note}; momentum={mom_note}"))
            sim_ok = obs_ok
        else:
            gates.append(GateResult("Gate D", "FAIL", "observable checks skipped because simulation failed"))

    payload = {
        "script": "validate_ipm.py",
        "inputs": {
            "icarion_bin": str(args.icarion_bin),
            "precompute_bin": str(args.precompute_bin),
            "build_dir": str(args.build_dir),
            "species": args.species,
            "gas": args.gas,
            "pressure_pa": args.pressure_pa,
            "temperature_k": args.temperature_k,
            "ion_count": args.ion_count,
            "dt_s": args.dt_s,
            "total_time_s": args.total_time_s,
            "n_orientations": args.n_orientations,
            "n_trials": args.n_trials,
            "v_bins": args.v_bins,
            "v_min": args.v_min,
            "v_max": args.v_max,
        },
        "gates": [gate.__dict__ for gate in gates],
        "artifacts": artifacts,
    }

    summary_json = OUT_DIR / "ipm_gate_summary.json"
    summary_md = OUT_DIR / "ipm_gate_summary.md"
    artifacts["summary_json"] = str(summary_json)
    artifacts["summary_md"] = str(summary_md)
    payload["artifacts"] = artifacts
    write_summary(summary_json, summary_md, payload)

    print(f"Wrote summary: {summary_json}")
    for gate in gates:
        print(f"{gate.gate}: {gate.status} - {gate.note}")

    failed = any(gate.status == "FAIL" for gate in gates)
    return 1 if failed else 0


if __name__ == "__main__":
    raise SystemExit(main())

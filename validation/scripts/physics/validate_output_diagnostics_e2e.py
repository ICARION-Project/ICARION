#!/usr/bin/env python3
"""E2E validation gate for output diagnostics features.

Covers v1.1 output-path additions with a runtime check:
- output.trajectory_mode=minimal -> /analysis/minimal_transport
- output.deep_analysis.mode=sampled_events -> /analysis/deep_collision
- reproducibility embed cap via ICARION_HDF5_MAX_EMBED_BYTES
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import h5py
import numpy as np


REPO_ROOT = Path(__file__).resolve().parents[3]
VALIDATION_DIR = REPO_ROOT / "validation"
SPECIES_DB = REPO_ROOT / "data" / "species_database_v1.json"
DEFAULT_BIN = REPO_ROOT / "build" / "src" / "icarion_main"

RUN_DIR_ENV = os.environ.get("ICARION_VALIDATION_RUN_DIR")
if RUN_DIR_ENV:
    RUN_DIR = Path(RUN_DIR_ENV)
    OUT_DIR = RUN_DIR / "results" / "physics" / "output_diagnostics"
    LOG_DIR = RUN_DIR / "logs" / "output_diagnostics"
else:
    OUT_DIR = VALIDATION_DIR / "results" / "physics" / "output_diagnostics"
    LOG_DIR = VALIDATION_DIR / "logs" / "output_diagnostics"


@dataclass
class CheckResult:
    name: str
    status: str
    note: str


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Validate output diagnostics E2E")
    ap.add_argument("--icarion-bin", type=Path, default=DEFAULT_BIN)
    ap.add_argument("--species", default="H3O+")
    ap.add_argument("--gas", default="He")
    ap.add_argument("--pressure-pa", type=float, default=1200.0)
    ap.add_argument("--temperature-k", type=float, default=300.0)
    ap.add_argument("--ion-count", type=int, default=400)
    ap.add_argument("--dt-s", type=float, default=4.0e-9)
    ap.add_argument("--total-time-s", type=float, default=1.0e-5)
    ap.add_argument("--write-interval", type=int, default=20)
    ap.add_argument("--embed-cap-bytes", type=int, default=1024)
    ap.add_argument("--timeout", type=int, default=900)
    return ap.parse_args()


def build_config(
    args: argparse.Namespace, out_folder: Path, traj_name: str, traj_mode: str = "minimal"
) -> dict[str, Any]:
    return {
        "simulation": {
            "total_time_s": args.total_time_s,
            "dt_s": args.dt_s,
            "write_interval": args.write_interval,
            "integrator": "RK4",
            "enable_openmp": True,
            "enable_gpu": False,
            "rng_seed": 42,
        },
        "physics": {
            "collision_model": "HSS",
            "enable_space_charge": False,
            "enable_reactions": False,
        },
        "species_database": str(SPECIES_DB),
        "output": {
            "folder": str(out_folder),
            "trajectory_file": traj_name,
            "trajectory_mode": traj_mode,
            "print_progress": False,
            "deep_analysis": {
                "mode": "sampled_events",
                "domain_filter_index": -1,
                "sample_every_n": 5,
                "max_events_per_ion": 200,
            },
        },
        "ions": {
            "species": [
                {
                    "id": args.species,
                    "count": args.ion_count,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.0],
                        "std": [4e-5, 4e-5, 4e-5],
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": args.temperature_k,
                    },
                }
            ]
        },
        "domains": [
            {
                "name": "output_diag_gate",
                "instrument": "IMS",
                "geometry": {
                    "origin_m": [0.0, 0.0, -0.02],
                    "length_m": 0.04,
                    "radius_m": 0.01,
                },
                "env": {
                    "temperature_K": args.temperature_k,
                    "pressure_Pa": args.pressure_pa,
                    "gas_species": args.gas,
                    "gas_velocity_m_s": [0.0, 0.0, 0.0],
                },
                "fields": {"DC": {"axial_V": 40.0}},
            }
        ],
    }


def run_case(
    args: argparse.Namespace,
    cfg_path: Path,
    log_path: Path,
    extra_env: dict[str, str] | None = None,
) -> tuple[int, float, str]:
    env = os.environ.copy()
    if extra_env:
        env.update(extra_env)

    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(args.icarion_bin), str(cfg_path), "--threads", "1"],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=args.timeout,
        env=env,
    )
    elapsed = time.perf_counter() - t0
    out = (proc.stdout or "") + ("\n" + proc.stderr if proc.stderr else "")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(out, encoding="utf-8")
    return proc.returncode, elapsed, out


def check_minimal_and_deep(h5_path: Path) -> list[CheckResult]:
    checks: list[CheckResult] = []

    with h5py.File(h5_path, "r") as f:
        # 1) metadata/config trajectory mode
        traj_mode = ""
        if "/metadata/config/trajectory_mode" in f:
            raw = f["/metadata/config/trajectory_mode"][()]
            if isinstance(raw, bytes):
                traj_mode = raw.decode("utf-8")
            else:
                traj_mode = str(raw)

        if traj_mode == "minimal":
            checks.append(CheckResult("trajectory_mode", "PASS", "metadata trajectory_mode=minimal"))
        else:
            checks.append(CheckResult("trajectory_mode", "FAIL", f"unexpected trajectory_mode={traj_mode!r}"))

        # 2) minimal transport group + finite arrays
        required_minimal = [
            "/analysis/minimal_transport/final_pos_x",
            "/analysis/minimal_transport/final_pos_y",
            "/analysis/minimal_transport/final_pos_z",
            "/analysis/minimal_transport/final_vel_x",
            "/analysis/minimal_transport/final_vel_y",
            "/analysis/minimal_transport/final_vel_z",
            "/analysis/minimal_transport/active",
        ]
        missing = [p for p in required_minimal if p not in f]
        if missing:
            checks.append(CheckResult("minimal_transport_presence", "FAIL", f"missing datasets: {missing}"))
        else:
            finites_ok = True
            n = None
            for p in required_minimal[:-1]:  # numeric arrays
                arr = np.asarray(f[p][:], dtype=float)
                if arr.size == 0 or not np.isfinite(arr).all():
                    finites_ok = False
                if n is None:
                    n = arr.size
                elif n != arr.size:
                    finites_ok = False
            checks.append(
                CheckResult(
                    "minimal_transport_finite",
                    "PASS" if finites_ok else "FAIL",
                    "minimal transport arrays finite and consistent" if finites_ok else "invalid or inconsistent minimal arrays",
                )
            )

        # 3) trajectory snapshots should be absent in minimal mode
        has_positions = "/trajectory/positions" in f
        checks.append(
            CheckResult(
                "trajectory_snapshots_disabled",
                "PASS" if not has_positions else "FAIL",
                "trajectory snapshots disabled in minimal mode" if not has_positions else "/trajectory/positions exists unexpectedly",
            )
        )

        # 4) deep collision summary presence
        deep_required = [
            "/analysis/deep_collision/mode",
            "/analysis/deep_collision/collisions_total",
            "/analysis/deep_collision/mean_delta_pz_kgms",
        ]
        deep_missing = [p for p in deep_required if p not in f]
        if deep_missing:
            checks.append(CheckResult("deep_collision_presence", "FAIL", f"missing deep datasets: {deep_missing}"))
        else:
            mode_raw = f["/analysis/deep_collision/mode"][()]
            mode = mode_raw.decode("utf-8") if isinstance(mode_raw, bytes) else str(mode_raw)
            collisions = np.asarray(f["/analysis/deep_collision/collisions_total"][:], dtype=np.int64)
            deep_ok = mode == "sampled_events" and collisions.size > 0
            checks.append(
                CheckResult(
                    "deep_collision_summary",
                    "PASS" if deep_ok else "FAIL",
                    f"mode={mode}, ions={collisions.size}, collisions_sum={int(collisions.sum())}",
                )
            )
            # events group is optional when no sampled events are retained
            has_events = "/analysis/deep_collision/events/time_s" in f
            checks.append(
                CheckResult(
                    "deep_collision_events",
                    "PASS" if has_events else "PARTIAL",
                    "sampled event stream present" if has_events else "events stream absent (allowed when no events retained)",
                )
            )

    return checks


def check_embed_cap(h5_path: Path, expected_cap: int) -> list[CheckResult]:
    checks: list[CheckResult] = []
    with h5py.File(h5_path, "r") as f:
        cap_path = "/metadata/reproducibility/embed_max_bytes"
        if cap_path not in f:
            checks.append(CheckResult("embed_cap_metadata", "FAIL", "embed_max_bytes missing"))
            return checks

        cap = int(np.asarray(f[cap_path][()]).item())
        checks.append(
            CheckResult(
                "embed_cap_value",
                "PASS" if cap == expected_cap else "FAIL",
                f"embed_max_bytes={cap}, expected={expected_cap}",
            )
        )

        # Heuristic evidence that cap had effect: at least one blob status exists.
        blob_root = "/metadata/reproducibility/input_blobs"
        if blob_root in f:
            statuses = []
            grp = f[blob_root]
            for key in grp.keys():
                if key.endswith("_status"):
                    raw = grp[key][()]
                    statuses.append(raw.decode("utf-8") if isinstance(raw, bytes) else str(raw))
            if statuses:
                note = ",".join(sorted(set(statuses))[:5])
                checks.append(CheckResult("embed_blob_statuses", "PASS", f"observed statuses: {note}"))
            else:
                checks.append(CheckResult("embed_blob_statuses", "PARTIAL", "no *_status datasets found under input_blobs"))
        else:
            checks.append(CheckResult("embed_blob_statuses", "FAIL", "input_blobs group missing"))

    return checks


def compare_minimal_vs_default(
    minimal_h5: Path, default_h5: Path, tol_rel: float = 1e-10, tol_abs: float = 1e-20
) -> list[CheckResult]:
    """Compare minimal-mode and full-mode outputs with identical seed.
    
    Verifies that minimal mode stores the same physical state as full mode
    for core trajectory endpoints and statistics.
    """
    checks: list[CheckResult] = []

    with h5py.File(minimal_h5, "r") as minimal_f, h5py.File(default_h5, "r") as default_f:
        # Check if both files have the required datasets
        required_minimal = [
            "/analysis/minimal_transport/final_pos_x",
            "/analysis/minimal_transport/final_pos_y",
            "/analysis/minimal_transport/final_pos_z",
            "/analysis/minimal_transport/final_vel_x",
            "/analysis/minimal_transport/final_vel_y",
            "/analysis/minimal_transport/final_vel_z",
            "/analysis/minimal_transport/active",
        ]

        # Full mode should have trajectory snapshots
        if "/trajectory/positions" not in default_f:
            checks.append(
                CheckResult("full_mode_trajectories", "FAIL", "/trajectory/positions missing in full mode")
            )
            return checks

        # Extract final states from full mode trajectories
        default_pos = default_f["/trajectory/positions"][-1, :, :]  # last snapshot, all ions
        default_vel = default_f["/trajectory/velocities"][-1, :, :]
        
        # Determine active status: ions with death_time_s == -1 are still alive
        death_times = np.asarray(default_f["/ions/death_time_s"][:], dtype=float)
        default_active = death_times < 0  # -1 means still alive
        
        if default_pos.shape[0] == 0:
            checks.append(CheckResult("full_mode_data", "FAIL", "full mode has no trajectory data"))
            return checks

        minimal_pos_x = np.asarray(minimal_f["/analysis/minimal_transport/final_pos_x"][:], dtype=float)
        minimal_pos_y = np.asarray(minimal_f["/analysis/minimal_transport/final_pos_y"][:], dtype=float)
        minimal_pos_z = np.asarray(minimal_f["/analysis/minimal_transport/final_pos_z"][:], dtype=float)
        minimal_vel_x = np.asarray(minimal_f["/analysis/minimal_transport/final_vel_x"][:], dtype=float)
        minimal_vel_y = np.asarray(minimal_f["/analysis/minimal_transport/final_vel_y"][:], dtype=float)
        minimal_vel_z = np.asarray(minimal_f["/analysis/minimal_transport/final_vel_z"][:], dtype=float)
        minimal_active = np.asarray(minimal_f["/analysis/minimal_transport/active"][:], dtype=bool)

        # Stack minimal positions
        minimal_pos = np.column_stack([minimal_pos_x, minimal_pos_y, minimal_pos_z])

        # Check shape consistency
        if minimal_pos.shape[0] != default_pos.shape[0]:
            checks.append(
                CheckResult(
                    "data_shape_consistency",
                    "FAIL",
                    f"ion count mismatch: minimal={minimal_pos.shape[0]}, full={default_pos.shape[0]}",
                )
            )
            return checks

        checks.append(
            CheckResult(
                "data_shape_consistency",
                "PASS",
                f"both modes have {minimal_pos.shape[0]} ions",
            )
        )

        # Compare positions
        pos_diff = np.abs(minimal_pos - default_pos[:, :3])
        pos_max_rel = np.max(np.absolute(pos_diff / (np.abs(default_pos[:, :3]) + 1e-20)))
        pos_max_abs = np.max(pos_diff)

        pos_ok = pos_max_rel < tol_rel or pos_max_abs < tol_abs
        checks.append(
            CheckResult(
                "final_position_agreement",
                "PASS" if pos_ok else "FAIL",
                f"max_rel={pos_max_rel:.2e}, max_abs={pos_max_abs:.2e}m",
            )
        )

        # Compare velocities
        minimal_vel = np.column_stack([minimal_vel_x, minimal_vel_y, minimal_vel_z])
        vel_diff = np.abs(minimal_vel - default_vel[:, :3])
        vel_max_rel = np.max(np.absolute(vel_diff / (np.abs(default_vel[:, :3]) + 1e-20)))
        vel_max_abs = np.max(vel_diff)

        vel_ok = vel_max_rel < tol_rel or vel_max_abs < tol_abs
        checks.append(
            CheckResult(
                "final_velocity_agreement",
                "PASS" if vel_ok else "FAIL",
                f"max_rel={vel_max_rel:.2e}, max_abs={vel_max_abs:.2e}m/s",
            )
        )

        # Compare active status (should be identical)
        active_match = np.array_equal(minimal_active, default_active)
        checks.append(
            CheckResult(
                "active_status_agreement",
                "PASS" if active_match else "FAIL",
                f"minimal active={int(minimal_active.sum())}/{len(minimal_active)}, full active={int(default_active.sum())}/{len(default_active)}",
            )
        )

        # Summary statistics comparison
        mean_pos_diff = np.mean(np.linalg.norm(pos_diff, axis=1))
        mean_vel_diff = np.mean(np.linalg.norm(vel_diff, axis=1))
        checks.append(
            CheckResult(
                "mean_state_difference",
                "PASS" if mean_pos_diff < 1e-15 and mean_vel_diff < 1e-15 else "PARTIAL",
                f"mean position delta={mean_pos_diff:.2e}m, mean velocity delta={mean_vel_diff:.2e}m/s",
            )
        )

    return checks


def summarize_status(checks: list[CheckResult]) -> str:
    statuses = {c.status for c in checks}
    if "FAIL" in statuses:
        return "FAIL"
    if "PARTIAL" in statuses:
        return "PARTIAL"
    return "PASS"


def main() -> int:
    args = parse_args()
    if not args.icarion_bin.exists():
        raise SystemExit(f"icarion binary not found: {args.icarion_bin}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    base_cfg = OUT_DIR / "output_diag_primary_config.json"
    base_h5 = OUT_DIR / "output_diag_primary.h5"
    base_log = LOG_DIR / "output_diag_primary.log"

    default_cfg = OUT_DIR / "output_diag_default_config.json"
    default_h5 = OUT_DIR / "output_diag_default.h5"
    default_log = LOG_DIR / "output_diag_default.log"

    cap_cfg = OUT_DIR / "output_diag_cap_config.json"
    cap_h5 = OUT_DIR / "output_diag_cap.h5"
    cap_log = LOG_DIR / "output_diag_cap.log"

    with open(base_cfg, "w", encoding="utf-8") as f:
        json.dump(build_config(args, OUT_DIR, base_h5.name, traj_mode="minimal"), f, indent=2)
    with open(default_cfg, "w", encoding="utf-8") as f:
        json.dump(build_config(args, OUT_DIR, default_h5.name, traj_mode="full"), f, indent=2)
    with open(cap_cfg, "w", encoding="utf-8") as f:
        json.dump(build_config(args, OUT_DIR, cap_h5.name, traj_mode="minimal"), f, indent=2)

    checks: list[CheckResult] = []

    rc, elapsed, _ = run_case(args, base_cfg, base_log)
    if rc != 0:
        checks.append(CheckResult("primary_run", "FAIL", f"simulation rc={rc}"))
    elif not base_h5.exists():
        checks.append(CheckResult("primary_run", "FAIL", "HDF5 output missing"))
    else:
        checks.append(CheckResult("primary_run", "PASS", f"completed in {elapsed:.2f}s"))
        checks.extend(check_minimal_and_deep(base_h5))

    rc_default, elapsed_default, _ = run_case(args, default_cfg, default_log)
    if rc_default != 0:
        checks.append(CheckResult("full_mode_run", "FAIL", f"simulation rc={rc_default}"))
    elif not default_h5.exists():
        checks.append(CheckResult("full_mode_run", "FAIL", "HDF5 output missing"))
    else:
        checks.append(CheckResult("full_mode_run", "PASS", f"completed in {elapsed_default:.2f}s"))
        # Compare minimal vs full if both succeeded
        if base_h5.exists() and default_h5.exists():
            checks.extend(compare_minimal_vs_default(base_h5, default_h5))

    rc2, elapsed2, _ = run_case(
        args,
        cap_cfg,
        cap_log,
        extra_env={"ICARION_HDF5_MAX_EMBED_BYTES": str(args.embed_cap_bytes)},
    )
    if rc2 != 0:
        checks.append(CheckResult("embed_cap_run", "FAIL", f"simulation rc={rc2}"))
    elif not cap_h5.exists():
        checks.append(CheckResult("embed_cap_run", "FAIL", "cap-run HDF5 output missing"))
    else:
        checks.append(CheckResult("embed_cap_run", "PASS", f"completed in {elapsed2:.2f}s"))
        checks.extend(check_embed_cap(cap_h5, args.embed_cap_bytes))

    status = summarize_status(checks)

    payload = {
        "script": "validate_output_diagnostics_e2e.py",
        "status": status,
        "inputs": {
            "species": args.species,
            "gas": args.gas,
            "pressure_pa": args.pressure_pa,
            "temperature_k": args.temperature_k,
            "ion_count": args.ion_count,
            "dt_s": args.dt_s,
            "total_time_s": args.total_time_s,
            "embed_cap_bytes": args.embed_cap_bytes,
        },
        "checks": [c.__dict__ for c in checks],
        "artifacts": {
            "primary_h5": str(base_h5),
            "primary_log": str(base_log),
            "default_h5": str(default_h5),
            "default_log": str(default_log),
            "cap_h5": str(cap_h5),
            "cap_log": str(cap_log),
        },
    }

    out_json = OUT_DIR / "output_diagnostics_gate_summary.json"
    out_md = OUT_DIR / "output_diagnostics_gate_summary.md"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)

    lines = [
        "# Output Diagnostics E2E Gate Summary",
        "",
        f"Status: {status}",
        "",
        "## Checks",
    ]
    for c in checks:
        lines.append(f"- {c.name}: {c.status} - {c.note}")
    lines += [
        "",
        "## Artifacts",
        f"- JSON: {out_json}",
        f"- Primary (Minimal) H5: {base_h5}",
        f"- Full Mode H5: {default_h5}",
        f"- Cap Test (Minimal) H5: {cap_h5}",
    ]
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote: {out_json}")
    print(f"Wrote: {out_md}")
    print(f"Gate: {status}")

    return 1 if status == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())

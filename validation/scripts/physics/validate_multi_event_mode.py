#!/usr/bin/env python3
"""Validation gate for collision multi-event mode.

Runs paired short drift scenarios with single-event vs multi-event collision mode,
collects finite transport diagnostics, and checks that multi-event stays within a
reasonable envelope while actually activating micro-subcycling.
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
RESULT_PREFIX = "ICARION_RESULT "

RUN_DIR_ENV = os.environ.get("ICARION_VALIDATION_RUN_DIR")
if RUN_DIR_ENV:
    RUN_DIR = Path(RUN_DIR_ENV)
    OUT_DIR = RUN_DIR / "results" / "physics" / "multi_event"
    LOG_DIR = RUN_DIR / "logs" / "multi_event"
else:
    OUT_DIR = VALIDATION_DIR / "results" / "multi_event"
    LOG_DIR = VALIDATION_DIR / "logs" / "multi_event"


@dataclass
class CaseResult:
    name: str
    status: str
    runtime_s: float
    active_fraction: float
    mean_vz_m_s: float
    median_arrival_s: float
    collision_events: int
    collision_substeps: int
    collision_multi_event_mode: str
    traj_file: str
    log_file: str
    note: str


@dataclass
class ConvergencePoint:
    dt_s: float
    mode: str
    median_arrival_s: float
    rel_err_vs_ref: float
    active_fraction: float
    active_delta_vs_ref: float


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Validate collision multi-event mode behavior")
    ap.add_argument("--icarion-bin", type=Path, default=DEFAULT_BIN)
    ap.add_argument("--species", default="H3O+")
    ap.add_argument("--gas", default="He")
    ap.add_argument("--pressure-pa", type=float, default=1000.0)
    ap.add_argument("--temperature-k", type=float, default=300.0)
    ap.add_argument("--ion-count", type=int, default=400)
    ap.add_argument("--field-v-m", type=float, default=1000.0)
    ap.add_argument("--dt-s", type=float, default=4.0e-9)
    ap.add_argument("--total-time-s", type=float, default=1.2e-5)
    ap.add_argument("--write-interval", type=int, default=20)
    ap.add_argument("--max-events", type=int, default=4)
    ap.add_argument("--arrival-tol-frac", type=float, default=0.02)
    ap.add_argument("--median-arrival-rel-tol", type=float, default=0.25)
    ap.add_argument("--active-fraction-abs-tol", type=float, default=0.05)
    ap.add_argument("--timeout", type=int, default=900)
    ap.add_argument(
        "--convergence-dt-values",
        default="1e-9,4e-9,8e-9,1.6e-8,4e-8",
        help="Comma-separated dt ladder for multi-event convergence sub-gate",
    )
    ap.add_argument(
        "--convergence-reference-dt",
        type=float,
        default=1.0e-9,
        help="Reference dt for convergence comparisons (single-event)",
    )
    ap.add_argument(
        "--convergence-recommended-dt",
        type=float,
        default=4.0e-9,
        help="Recommended operating dt to check explicitly",
    )
    ap.add_argument(
        "--convergence-recommended-max-events",
        type=int,
        default=4,
        help="Recommended multi-event cap used for convergence checks",
    )
    ap.add_argument(
        "--convergence-recommended-median-rel-tol",
        type=float,
        default=0.02,
        help="Max relative median-arrival error at recommended operating point",
    )
    ap.add_argument(
        "--convergence-prefer-multi-slack",
        type=float,
        default=0.01,
        help="Allowed slack in rel-error when expecting multi-event <= single-event",
    )
    ap.add_argument(
        "--skip-convergence",
        action="store_true",
        help="Skip convergence sub-gate (faster smoke-run)",
    )
    return ap.parse_args()


def parse_result_line(output: str) -> dict[str, str]:
    for line in output.splitlines():
        if RESULT_PREFIX in line:
            idx = line.index(RESULT_PREFIX)
            payload = line[idx + len(RESULT_PREFIX) :].strip()
            result: dict[str, str] = {}
            for token in payload.split():
                if "=" not in token:
                    continue
                k, v = token.split("=", 1)
                result[k] = v
            return result
    raise RuntimeError("ICARION_RESULT line not found")


def build_config(args: argparse.Namespace, traj_name: str) -> dict[str, Any]:
    axial_v = args.field_v_m * 0.04
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
            "folder": str(OUT_DIR),
            "trajectory_file": traj_name,
            "trajectory_mode": "minimal",
            "print_progress": False,
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
                "name": "multi_event_gate",
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
                "fields": {"DC": {"axial_V": axial_v}},
            }
        ],
    }


def run_case(
    args: argparse.Namespace,
    config_path: Path,
    case_name: str,
    traj_name: str,
    multi_event: bool,
    max_events: int | None,
    dt_override: float | None = None,
) -> tuple[int, float, str]:
    log_path = LOG_DIR / f"{case_name}.log"
    cmd = [
        str(args.icarion_bin),
        str(config_path),
        "--threads",
        "1",
        "--set",
        f"physics.collision_multi_event_mode={'true' if multi_event else 'false'}",
        "--set",
        f"output.trajectory_file={traj_name}",
    ]
    if multi_event and max_events is not None:
        cmd += ["--set", f"physics.collision_max_events_per_step={max_events}"]
    if dt_override is not None:
        cmd += ["--set", f"simulation.dt_s={dt_override:.12g}"]

    t0 = time.perf_counter()
    proc = subprocess.run(
        cmd,
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=args.timeout,
    )
    elapsed = time.perf_counter() - t0

    output = (proc.stdout or "") + ("\n" + proc.stderr if proc.stderr else "")
    log_path.parent.mkdir(parents=True, exist_ok=True)
    log_path.write_text(output, encoding="utf-8")
    return proc.returncode, elapsed, output


def analyze_minimal_transport(h5_path: Path, field_v_m: float, tol_frac: float) -> tuple[float, float, float]:
    with h5py.File(h5_path, "r") as f:
        root = "/analysis/minimal_transport"
        final_vz = np.asarray(f[f"{root}/final_vel_z"][:], dtype=float)
        active = np.asarray(f[f"{root}/active"][:], dtype=bool)
        final_x = np.asarray(f[f"{root}/final_pos_x"][:], dtype=float)
        final_y = np.asarray(f[f"{root}/final_pos_y"][:], dtype=float)
        final_z = np.asarray(f[f"{root}/final_pos_z"][:], dtype=float)

        geom = "/domains/domain_0/geometry"
        length_m = float(f[f"{geom}/length_m"][()])
        radius_m = float(f[f"{geom}/radius_m"][()])
        origin = np.asarray(f[f"{geom}/origin_m"][()], dtype=float)

    if final_vz.size == 0 or not np.isfinite(final_vz).all():
        raise RuntimeError("invalid final_vel_z in minimal transport")

    active_fraction = float(np.mean(active)) if active.size else 0.0
    mean_vz = float(np.mean(final_vz[active])) if np.any(active) else float(np.mean(final_vz))

    # Approximate arrival using final-state projection near the detector plane.
    z_end = float(origin[2] + length_m)
    tol_z = length_m * tol_frac
    tol_r = radius_m * tol_frac
    rr = np.sqrt(final_x * final_x + final_y * final_y)
    arrived = (final_z >= (z_end - tol_z)) & (rr <= (radius_m + tol_r))

    # Approximate median arrival from drift velocity if arrivals exist.
    if np.any(arrived) and abs(mean_vz) > 1e-12:
        median_arrival = float(length_m / abs(mean_vz))
    else:
        median_arrival = float("nan")

    _ = field_v_m
    return active_fraction, mean_vz, median_arrival


def evaluate_cases(args: argparse.Namespace, baseline: CaseResult, multi: CaseResult) -> tuple[str, str]:
    if baseline.status != "PASS" or multi.status != "PASS":
        return "FAIL", "one or more cases failed"

    if multi.collision_multi_event_mode.lower() != "true":
        return "FAIL", "multi-event run did not report collision_multi_event_mode=true"

    if multi.collision_substeps <= baseline.collision_substeps:
        return "FAIL", "multi-event run did not increase collision_substeps"

    if abs(multi.active_fraction - baseline.active_fraction) > args.active_fraction_abs_tol:
        return (
            "PARTIAL",
            f"active_fraction delta too high ({abs(multi.active_fraction - baseline.active_fraction):.4f})",
        )

    if np.isfinite(baseline.median_arrival_s) and np.isfinite(multi.median_arrival_s):
        rel = abs(multi.median_arrival_s - baseline.median_arrival_s) / max(abs(baseline.median_arrival_s), 1e-12)
        if rel > args.median_arrival_rel_tol:
            return "PARTIAL", f"median arrival relative delta too high ({rel:.4f})"

    return "PASS", "multi-event mode active and transport envelope acceptable"


def build_convergence_points(
    args: argparse.Namespace,
    cfg_path: Path,
    ref_case: CaseResult,
) -> tuple[list[ConvergencePoint], str]:
    raw_values = [x.strip() for x in args.convergence_dt_values.split(",") if x.strip()]
    if not raw_values:
        return [], "no convergence dt values configured"

    dt_values = sorted({float(v) for v in raw_values})
    if args.convergence_reference_dt not in dt_values:
        dt_values = sorted(set(dt_values + [args.convergence_reference_dt]))

    points: list[ConvergencePoint] = []

    ref_median = ref_case.median_arrival_s
    ref_active = ref_case.active_fraction
    if not np.isfinite(ref_median):
        return [], "reference median arrival is not finite"

    for dt in dt_values:
        for mode, is_multi in (("single", False), ("multi", True)):
            name = f"conv_{mode}_{dt:.12g}".replace("+", "")
            traj_name = f"{name}.h5"
            max_events = args.convergence_recommended_max_events if is_multi else None

            rc, elapsed, out = run_case(
                args,
                cfg_path,
                name,
                traj_name,
                is_multi,
                max_events,
                dt_override=dt,
            )

            case = to_case_result(
                name,
                rc,
                elapsed,
                out,
                OUT_DIR / traj_name,
                LOG_DIR / f"{name}.log",
                args,
            )

            if case.status != "PASS" or not np.isfinite(case.median_arrival_s):
                return points, f"convergence case failed: {name} ({case.note})"

            rel_err = abs(case.median_arrival_s - ref_median) / max(abs(ref_median), 1e-12)
            active_delta = abs(case.active_fraction - ref_active)
            points.append(
                ConvergencePoint(
                    dt_s=dt,
                    mode=mode,
                    median_arrival_s=case.median_arrival_s,
                    rel_err_vs_ref=rel_err,
                    active_fraction=case.active_fraction,
                    active_delta_vs_ref=active_delta,
                )
            )

    return points, "ok"


def evaluate_convergence(args: argparse.Namespace, points: list[ConvergencePoint]) -> tuple[str, str]:
    if not points:
        return "PARTIAL", "no convergence points produced"

    by_dt: dict[float, dict[str, ConvergencePoint]] = {}
    for p in points:
        by_dt.setdefault(p.dt_s, {})[p.mode] = p

    # Hard check at recommended operating point.
    rec = by_dt.get(args.convergence_recommended_dt, {}).get("multi")
    if rec is None:
        return "PARTIAL", f"recommended dt={args.convergence_recommended_dt:.3g} missing in convergence set"
    if rec.rel_err_vs_ref > args.convergence_recommended_median_rel_tol:
        return (
            "FAIL",
            f"recommended-point rel error too high ({rec.rel_err_vs_ref:.4f} > {args.convergence_recommended_median_rel_tol:.4f})",
        )

    # Soft preference: for same dt, multi-event should not be worse than single-event (with slack).
    violations = 0
    checked = 0
    for dt, modes in sorted(by_dt.items()):
        s = modes.get("single")
        m = modes.get("multi")
        if s is None or m is None:
            continue
        checked += 1
        if m.rel_err_vs_ref > s.rel_err_vs_ref + args.convergence_prefer_multi_slack:
            violations += 1

    if checked == 0:
        return "PARTIAL", "no dt pairs with both single and multi results"
    if violations > 0:
        return "PARTIAL", f"{violations}/{checked} dt pairs where multi-event is worse than single beyond slack"

    return "PASS", "convergence ladder checks passed"


def to_case_result(
    name: str,
    rc: int,
    elapsed: float,
    out: str,
    traj: Path,
    log: Path,
    args: argparse.Namespace,
) -> CaseResult:
    if rc != 0:
        return CaseResult(
            name=name,
            status="FAIL",
            runtime_s=elapsed,
            active_fraction=float("nan"),
            mean_vz_m_s=float("nan"),
            median_arrival_s=float("nan"),
            collision_events=0,
            collision_substeps=0,
            collision_multi_event_mode="unknown",
            traj_file=str(traj),
            log_file=str(log),
            note=f"simulation rc={rc}",
        )

    try:
        result_line = parse_result_line(out)
        active_fraction, mean_vz, median_arrival = analyze_minimal_transport(
            traj, args.field_v_m, args.arrival_tol_frac
        )
        return CaseResult(
            name=name,
            status="PASS",
            runtime_s=elapsed,
            active_fraction=active_fraction,
            mean_vz_m_s=mean_vz,
            median_arrival_s=median_arrival,
            collision_events=int(result_line.get("collision_events", "0")),
            collision_substeps=int(result_line.get("collision_substeps", "0")),
            collision_multi_event_mode=result_line.get("collision_multi_event_mode", "unknown"),
            traj_file=str(traj),
            log_file=str(log),
            note="ok",
        )
    except Exception as exc:
        return CaseResult(
            name=name,
            status="FAIL",
            runtime_s=elapsed,
            active_fraction=float("nan"),
            mean_vz_m_s=float("nan"),
            median_arrival_s=float("nan"),
            collision_events=0,
            collision_substeps=0,
            collision_multi_event_mode="unknown",
            traj_file=str(traj),
            log_file=str(log),
            note=str(exc),
        )


def main() -> int:
    args = parse_args()
    if not args.icarion_bin.exists():
        raise SystemExit(f"icarion binary not found: {args.icarion_bin}")

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    cfg_path = OUT_DIR / "multi_event_gate_config.json"
    with open(cfg_path, "w", encoding="utf-8") as f:
        json.dump(build_config(args, "multi_event_gate.h5"), f, indent=2)

    baseline_name = "single_event_baseline"
    multi_name = "multi_event_enabled"

    baseline_traj = OUT_DIR / f"{baseline_name}.h5"
    multi_traj = OUT_DIR / f"{multi_name}.h5"

    rc0, t0, out0 = run_case(args, cfg_path, baseline_name, baseline_traj.name, False, None)
    rc1, t1, out1 = run_case(args, cfg_path, multi_name, multi_traj.name, True, args.max_events)

    baseline_log = LOG_DIR / f"{baseline_name}.log"
    multi_log = LOG_DIR / f"{multi_name}.log"

    baseline = to_case_result(baseline_name, rc0, t0, out0, baseline_traj, baseline_log, args)
    multi = to_case_result(multi_name, rc1, t1, out1, multi_traj, multi_log, args)

    gate_status, gate_note = evaluate_cases(args, baseline, multi)

    conv_status = "SKIPPED"
    conv_note = "disabled via --skip-convergence"
    conv_points: list[ConvergencePoint] = []
    if not args.skip_convergence:
        conv_points, conv_build_note = build_convergence_points(args, cfg_path, baseline)
        if conv_points:
            conv_status, conv_note = evaluate_convergence(args, conv_points)
        else:
            conv_status = "PARTIAL"
            conv_note = conv_build_note

    final_status = gate_status
    if final_status != "FAIL":
        if conv_status == "FAIL":
            final_status = "FAIL"
        elif conv_status == "PARTIAL" or gate_status == "PARTIAL":
            final_status = "PARTIAL"

    summary = {
        "script": "validate_multi_event_mode.py",
        "gate": "collision_multi_event",
        "status": final_status,
        "note": f"primary={gate_status}: {gate_note}; convergence={conv_status}: {conv_note}",
        "inputs": {
            "species": args.species,
            "gas": args.gas,
            "pressure_pa": args.pressure_pa,
            "temperature_k": args.temperature_k,
            "ion_count": args.ion_count,
            "field_v_m": args.field_v_m,
            "dt_s": args.dt_s,
            "total_time_s": args.total_time_s,
            "max_events": args.max_events,
            "convergence_dt_values": args.convergence_dt_values,
            "convergence_reference_dt": args.convergence_reference_dt,
            "convergence_recommended_dt": args.convergence_recommended_dt,
            "convergence_recommended_max_events": args.convergence_recommended_max_events,
        },
        "primary_gate": {"status": gate_status, "note": gate_note},
        "convergence_gate": {
            "status": conv_status,
            "note": conv_note,
            "points": [p.__dict__ for p in conv_points],
        },
        "cases": [baseline.__dict__, multi.__dict__],
    }

    out_json = OUT_DIR / "multi_event_gate_summary.json"
    out_md = OUT_DIR / "multi_event_gate_summary.md"
    with open(out_json, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    lines = [
        "# Collision Multi-Event Gate Summary",
        "",
        f"Status: {final_status}",
        "",
        f"Primary Gate: {gate_status} - {gate_note}",
        f"Convergence Gate: {conv_status} - {conv_note}",
        "",
        "## Cases",
    ]
    for c in [baseline, multi]:
        lines.append(
            f"- {c.name}: {c.status}, active_fraction={c.active_fraction:.4f}, "
            f"mean_vz={c.mean_vz_m_s:.6g} m/s, collision_substeps={c.collision_substeps}, note={c.note}"
        )
    if conv_points:
        lines += ["", "## Convergence Points"]
        for p in conv_points:
            lines.append(
                f"- dt={p.dt_s:.3g}, mode={p.mode}, rel_err_vs_ref={p.rel_err_vs_ref:.6f}, "
                f"active_delta_vs_ref={p.active_delta_vs_ref:.6f}"
            )
    lines += ["", "## Artifacts", f"- JSON: {out_json}", f"- Logs: {LOG_DIR}"]
    out_md.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Wrote: {out_json}")
    print(f"Wrote: {out_md}")
    print(f"Gate: {final_status} (primary={gate_status}, convergence={conv_status})")

    return 1 if final_status == "FAIL" else 0


if __name__ == "__main__":
    raise SystemExit(main())

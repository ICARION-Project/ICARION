#!/usr/bin/env python3
"""Validate zero-field diffusion against the Einstein relation.

Runs ICARION transport simulations with no imposed electric field and compares
the measured diffusion coefficient from the mean-squared displacement (MSD)
against the Einstein estimate

    D = k_B * T * K / q,   K = K0 * (P0 / P) * (T / T0)

for selected species and collision models.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import subprocess
import time
from dataclasses import dataclass
from pathlib import Path

import h5py
import matplotlib.pyplot as plt
import numpy as np


K_B = 1.380649e-23
Q_E = 1.602176634e-19
P0_PA = 101325.0
T0_K = 273.15

REPO_ROOT = Path(__file__).resolve().parents[3]
VALIDATION_DIR = REPO_ROOT / "validation"
SPECIES_DB_PATH = REPO_ROOT / "data" / "species_database_v1.json"
ICARION_BIN_DEFAULT = REPO_ROOT / "build" / "src" / "icarion_main"

RUN_DIR = os.environ.get("ICARION_VALIDATION_RUN_DIR")
if RUN_DIR:
    RUN_DIR = Path(RUN_DIR)
    OUT_ROOT = RUN_DIR / "results" / "physics" / "diffusion"
    LOG_DIR = RUN_DIR / "logs"
    FIG_DIR = RUN_DIR / "figures" / "physics" / "diffusion"
else:
    OUT_ROOT = VALIDATION_DIR / "results" / "diffusion_einstein"
    LOG_DIR = VALIDATION_DIR / "logs"
    FIG_DIR = VALIDATION_DIR / "figures" / "physics" / "diffusion_einstein"


@dataclass
class CaseResult:
    model: str
    species: str
    gas_species: str
    pressure_pa: float
    temperature_k: float
    ion_count: int
    total_time_s: float
    runtime_s: float
    k0_m2_vs: float
    k_m2_vs: float
    d_sim_m2_s: float
    d_ein_m2_s: float
    d_ratio: float
    msd_fit_r2: float
    fit_window: str
    status: str
    note: str
    trajectory_file: str
    plot_file: str


def parse_args() -> argparse.Namespace:
    ap = argparse.ArgumentParser(description="Validate diffusion against the Einstein relation")
    ap.add_argument("--icarion-bin", type=Path, default=ICARION_BIN_DEFAULT)
    ap.add_argument("--species", nargs="+", default=["H3O+", "PentanalH+"])
    ap.add_argument("--models", nargs="+", default=["HSS", "EHSS"])
    ap.add_argument("--gas-species", default="He")
    ap.add_argument("--pressure-pa", type=float, default=0.2)
    ap.add_argument("--temperature-k", type=float, default=300.0)
    ap.add_argument("--total-time-s", type=float, default=2.8e-3)
    ap.add_argument("--dt-s", type=float, default=1.0e-8)
    ap.add_argument("--write-interval", type=int, default=1000)
    ap.add_argument("--ion-count", type=int, default=400)
    ap.add_argument("--fit-start-frac", type=float, default=0.2)
    ap.add_argument("--fit-end-frac", type=float, default=0.8)
    ap.add_argument("--timeout", type=int, default=1800)
    return ap.parse_args()


def _load_species_db() -> dict:
    with open(SPECIES_DB_PATH, "r", encoding="utf-8") as f:
        return json.load(f)["species"]


def _load_k0_si(species_db: dict, species_id: str) -> float:
    mobility_cm2_vs = species_db[species_id]["mobility_cm2Vs"]
    return float(mobility_cm2_vs) * 1e-4


def _mobility_at_conditions(k0_m2_vs: float, pressure_pa: float, temperature_k: float) -> float:
    return k0_m2_vs * (P0_PA / pressure_pa) * (temperature_k / T0_K)


def _einstein_diffusion(k_m2_vs: float, temperature_k: float, charge_c: float = Q_E) -> float:
    return K_B * temperature_k * k_m2_vs / charge_c


def _case_slug(model: str, species: str, pressure_pa: float, temperature_k: float) -> str:
    species_slug = species.replace("+", "plus").replace(" ", "_")
    return f"diffusion_{model}_{species_slug}_{int(round(temperature_k))}K_{pressure_pa:g}Pa"


def _create_config(
    species: str,
    model: str,
    gas_species: str,
    pressure_pa: float,
    temperature_k: float,
    total_time_s: float,
    dt_s: float,
    write_interval: int,
    ion_count: int,
    output_name: str,
) -> dict:
    return {
        "simulation": {
            "total_time_s": total_time_s,
            "dt_s": dt_s,
            "write_interval": write_interval,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": 42,
        },
        "physics": {
            "collision_model": model,
            "enable_space_charge": False,
            "enable_reactions": False,
        },
        "species_database_path": str(SPECIES_DB_PATH),
        "output": {
            "folder": str(OUT_ROOT),
            "trajectory_file": output_name,
            "print_progress": False,
        },
        "ions": {
            "species": [
                {
                    "species_id": species,
                    "count": ion_count,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.0],
                        "std": [2e-4, 2e-4, 2e-4],
                    },
                    "velocity": {
                        "type": "thermal",
                        "temperature_K": temperature_k,
                    },
                }
            ]
        },
        "domains": [
            {
                "domain_index": 0,
                "name": f"Diffusion {model} {species}",
                "instrument": "IMS",
                "geometry": {
                    "type": "cylinder",
                    "radius_m": 0.05,
                    "length_m": 0.20,
                    "origin_m": [0.0, 0.0, -0.10],
                },
                "env": {
                    "gas_species": gas_species,
                    "temperature_K": temperature_k,
                    "pressure_Pa": pressure_pa,
                    "gas_velocity_m_s": [0.0, 0.0, 0.0],
                },
                "fields": {
                    "DC": {
                        "uniform_field_V_m": [0.0, 0.0, 0.0],
                    }
                },
            }
        ],
    }


def _run_case(icarion_bin: Path, config_path: Path, timeout_s: int) -> tuple[int, float, str]:
    t0 = time.perf_counter()
    proc = subprocess.run(
        [str(icarion_bin), str(config_path)],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
        timeout=timeout_s,
    )
    runtime_s = time.perf_counter() - t0
    log = (proc.stdout or "") + ("\n" + proc.stderr if proc.stderr else "")
    return proc.returncode, runtime_s, log


def _analyze_diffusion(h5_path: Path, fit_start_frac: float, fit_end_frac: float) -> tuple[float, float, str, np.ndarray, np.ndarray, np.ndarray]:
    with h5py.File(h5_path, "r") as f:
        time_s = np.asarray(f["/trajectory/time"][:], dtype=float)
        positions = np.asarray(f["/trajectory/positions"][:], dtype=float)

    disp = positions - positions[0:1, :, :]
    disp_centered = disp - np.mean(disp, axis=1, keepdims=True)
    msd = np.mean(np.sum(disp_centered * disp_centered, axis=2), axis=1)

    if len(time_s) < 10:
        raise RuntimeError(f"Not enough trajectory frames for diffusion fit: {len(time_s)}")

    start_idx = max(1, int(len(time_s) * fit_start_frac))
    end_idx = max(start_idx + 3, int(len(time_s) * fit_end_frac))
    end_idx = min(end_idx, len(time_s))

    fit_t = time_s[start_idx:end_idx]
    fit_msd = msd[start_idx:end_idx]
    slope, intercept = np.polyfit(fit_t, fit_msd, 1)
    fit_pred = slope * fit_t + intercept
    ss_res = float(np.sum((fit_msd - fit_pred) ** 2))
    ss_tot = float(np.sum((fit_msd - np.mean(fit_msd)) ** 2))
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0.0 else 1.0
    d_sim = slope / 6.0
    return d_sim, r2, f"{start_idx}:{end_idx}", time_s, msd, slope * time_s + intercept


def _plot_case(plot_path: Path, time_s: np.ndarray, msd: np.ndarray, fit_line: np.ndarray, d_sim: float, d_ein: float, r2: float, label: str) -> None:
    plot_path.parent.mkdir(parents=True, exist_ok=True)
    plt.figure(figsize=(7, 4.5))
    plt.plot(time_s * 1e3, msd * 1e6, label="MSD", lw=2)
    plt.plot(time_s * 1e3, fit_line * 1e6, "--", label="Linear fit", lw=2)
    plt.xlabel("Time [ms]")
    plt.ylabel("MSD [mm$^2$]")
    plt.title(label)
    plt.legend(title=f"D_sim={d_sim:.3e} m$^2$/s\nD_ein={d_ein:.3e} m$^2$/s\nR$^2$={r2:.4f}")
    plt.tight_layout()
    plt.savefig(plot_path, dpi=160)
    plt.close()


def _write_csv(results: list[CaseResult], csv_path: Path) -> None:
    csv_path.parent.mkdir(parents=True, exist_ok=True)
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([
            "model", "species", "gas_species", "P_Pa", "T_K", "ion_count", "total_time_s",
            "runtime_s", "K0_m2Vs", "K_m2Vs", "D_sim_m2s", "D_ein_m2s", "D_ratio",
            "msd_fit_r2", "fit_window", "status", "note", "trajectory_file", "plot_file",
        ])
        for row in results:
            writer.writerow([
                row.model, row.species, row.gas_species, row.pressure_pa, row.temperature_k,
                row.ion_count, row.total_time_s, row.runtime_s, row.k0_m2_vs, row.k_m2_vs,
                row.d_sim_m2_s, row.d_ein_m2_s, row.d_ratio, row.msd_fit_r2, row.fit_window,
                row.status, row.note, row.trajectory_file, row.plot_file,
            ])


def _write_report(results: list[CaseResult], report_path: Path, args: argparse.Namespace) -> None:
    report_path.parent.mkdir(parents=True, exist_ok=True)
    ok = sum(1 for r in results if r.status == "ok")
    lines = [
        "# Diffusion vs Einstein Validation",
        "",
        f"- gas_species: `{args.gas_species}`",
        f"- pressure_Pa: `{args.pressure_pa}`",
        f"- temperature_K: `{args.temperature_k}`",
        f"- dt_s: `{args.dt_s}`",
        f"- total_time_s: `{args.total_time_s}`",
        f"- ion_count: `{args.ion_count}`",
        f"- fit_window_frac: `{args.fit_start_frac}`:`{args.fit_end_frac}`",
        f"- cases: `{len(results)}`",
        f"- ok: `{ok}`",
        "",
        "| Model | Species | Status | D_sim [m^2/s] | D_ein [m^2/s] | Ratio | R^2 | Fit Window |",
        "|---|---|---|---:|---:|---:|---:|---|",
    ]
    for r in results:
        lines.append(
            f"| {r.model} | {r.species} | {r.status} | {r.d_sim_m2_s:.6e} | {r.d_ein_m2_s:.6e} | {r.d_ratio:.3f} | {r.msd_fit_r2:.4f} | {r.fit_window} |"
        )
    lines.append("")
    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    if not args.icarion_bin.exists():
        raise SystemExit(f"Missing ICARION binary: {args.icarion_bin}")

    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    LOG_DIR.mkdir(parents=True, exist_ok=True)
    FIG_DIR.mkdir(parents=True, exist_ok=True)

    species_db = _load_species_db()
    results: list[CaseResult] = []

    for model in args.models:
        for species in args.species:
            if species not in species_db:
                results.append(
                    CaseResult(model, species, args.gas_species, args.pressure_pa, args.temperature_k,
                               args.ion_count, args.total_time_s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, "n/a", "missing_species", "species not in species database", "", "")
                )
                continue

            slug = _case_slug(model, species, args.pressure_pa, args.temperature_k)
            cfg_path = OUT_ROOT / f"{slug}.json"
            h5_path = OUT_ROOT / f"{slug}.h5"
            plot_path = FIG_DIR / f"{slug}.png"
            cfg = _create_config(
                species=species,
                model=model,
                gas_species=args.gas_species,
                pressure_pa=args.pressure_pa,
                temperature_k=args.temperature_k,
                total_time_s=args.total_time_s,
                dt_s=args.dt_s,
                write_interval=args.write_interval,
                ion_count=args.ion_count,
                output_name=h5_path.name,
            )
            cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")

            try:
                rc, runtime_s, log = _run_case(args.icarion_bin, cfg_path, args.timeout)
                if rc != 0:
                    results.append(
                        CaseResult(model, species, args.gas_species, args.pressure_pa, args.temperature_k,
                                   args.ion_count, args.total_time_s, runtime_s, 0.0, 0.0, 0.0, 0.0, 0.0,
                                   0.0, "n/a", "failed", log.strip()[-300:], str(h5_path), str(plot_path))
                    )
                    continue

                k0 = _load_k0_si(species_db, species)
                k = _mobility_at_conditions(k0, args.pressure_pa, args.temperature_k)
                d_ein = _einstein_diffusion(k, args.temperature_k)
                d_sim, r2, fit_window, time_s, msd, fit_line = _analyze_diffusion(
                    h5_path, args.fit_start_frac, args.fit_end_frac
                )
                ratio = d_sim / d_ein if d_ein > 0.0 else 0.0
                _plot_case(plot_path, time_s, msd, fit_line, d_sim, d_ein, r2, f"{model} {species}")
                results.append(
                    CaseResult(model, species, args.gas_species, args.pressure_pa, args.temperature_k,
                               args.ion_count, args.total_time_s, runtime_s, k0, k, d_sim, d_ein, ratio,
                               r2, fit_window, "ok", "", str(h5_path), str(plot_path))
                )
            except subprocess.TimeoutExpired:
                results.append(
                    CaseResult(model, species, args.gas_species, args.pressure_pa, args.temperature_k,
                               args.ion_count, args.total_time_s, float(args.timeout), 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, "n/a", "timeout", f"timeout after {args.timeout}s", str(h5_path), str(plot_path))
                )
            except Exception as exc:
                results.append(
                    CaseResult(model, species, args.gas_species, args.pressure_pa, args.temperature_k,
                               args.ion_count, args.total_time_s, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               0.0, "n/a", "failed", str(exc), str(h5_path), str(plot_path))
                )

    csv_path = OUT_ROOT / "diffusion_results.csv"
    report_path = OUT_ROOT / "diffusion_report.md"
    _write_csv(results, csv_path)
    _write_report(results, report_path, args)
    print(f"Wrote: {csv_path}")
    print(f"Wrote: {report_path}")
    failures = [r for r in results if r.status != "ok"]
    return 1 if failures else 0


if __name__ == "__main__":
    raise SystemExit(main())
#!/usr/bin/env python3
"""ICARION Physics Validation: Friction+OU vs HSD+OU Thermalization Comparison.

Both models use deterministic damping + Ornstein-Uhlenbeck thermal noise.
Starting from cold (0.1 K), ions should thermalize to the gas temperature (300 K).

Tolerance: final ion temperature within ±20% of T_gas.
Fixed seed: 42 (reproducibility).

Test species: H3O+ in He
Conditions: P=100 Pa, T=300 K, E=0 V/m
Duration: 5 µs (~38 collision times for H3O+/He at 100 Pa)
"""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import h5py
import numpy as np

COMMON_DIR = Path(__file__).resolve().parents[1] / "common"
if str(COMMON_DIR) not in sys.path:
    sys.path.append(str(COMMON_DIR))
from hdf5_utils import load_species_ids  # noqa: E402

# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------
T_GAS_K = 300.0
P_PA = 100.0
ION_SPECIES = "H3O+"
ION_MASS_AMU = 19.02
N_IONS = 500
DT_S = 1e-8
TOTAL_TIME_S = 5e-6
WRITE_INTERVAL = 50
RNG_SEED = 42
TOLERANCE_PCT = 20.0  # documented tolerance

MODELS = [
    {"label": "Friction+OU", "collision_model": "Friction", "enable_ou_thermalization": True},
    {"label": "HSD+OU",      "collision_model": "HSD",      "enable_ou_thermalization": True},
]

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_ROOT = Path(__file__).resolve().parents[3]
ICARION_BIN = Path(os.environ.get("ICARION_BIN", PROJECT_ROOT / "build" / "src" / "icarion_main"))
SPECIES_DB = PROJECT_ROOT / "data" / "species_database_v1.json"

import os
_RUN_DIR = os.environ.get("ICARION_VALIDATION_RUN_DIR")
if _RUN_DIR:
    RESULTS_DIR = Path(_RUN_DIR) / "results" / "ou_comparison"
else:
    RESULTS_DIR = PROJECT_ROOT / "validation" / "results" / "ou_thermalization_comparison"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

# ---------------------------------------------------------------------------
# Physical constants
# ---------------------------------------------------------------------------
K_B = 1.380649e-23
AMU_TO_KG = 1.66053906660e-27


def _build_config(model_cfg: dict, h5_name: str) -> dict:
    return {
        "simulation": {
            "total_time_s": TOTAL_TIME_S,
            "dt_s": DT_S,
            "write_interval": WRITE_INTERVAL,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": RNG_SEED,
        },
        "physics": {
            "collision_model": model_cfg["collision_model"],
            "enable_ou_thermalization": model_cfg["enable_ou_thermalization"],
            "enable_space_charge": False,
            "enable_reactions": False,
        },
        "species_database": str(SPECIES_DB),
        "output": {
            "folder": str(RESULTS_DIR),
            "trajectory_file": h5_name,
            "print_progress": False,
        },
        "ions": {
            "species": [
                {
                    "species_id": ION_SPECIES,
                    "count": N_IONS,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.0],
                        "std": [0.001, 0.001, 0.001],
                    },
                    "velocity": {"type": "thermal", "temperature_K": 0.1},
                }
            ]
        },
        "domains": [
            {
                "name": "thermalization_chamber",
                "instrument": "IMS",
                "geometry": {
                    "origin_m": [0.0, 0.0, -5.0],
                    "radius_m": 10,
                    "length_m": 10,
                },
                "environment": {
                    "pressure_Pa": P_PA,
                    "temperature_K": T_GAS_K,
                    "gas_species": ["He"],
                },
                "fields": {
                    "electric": {"type": "none"},
                    "magnetic": {"type": "none"},
                },
            }
        ],
    }


def _analyze(h5_path: Path) -> dict:
    with h5py.File(h5_path, "r") as f:
        velocities = f["trajectory/velocities"][:]  # (frames, ions, 3)
        times = f["trajectory/time"][:]
        species_ids = load_species_ids(f)
        if species_ids.ndim == 1:
            species_ids = species_ids[np.newaxis, :]

    n_frames, n_ions, _ = velocities.shape
    T_ions = np.zeros(n_frames)

    for i in range(n_frames):
        active = np.array([s != "" for s in species_ids[i]])
        if active.sum() > 0:
            v2 = np.sum(velocities[i, active] ** 2, axis=1)
            T_ions[i] = ION_MASS_AMU * AMU_TO_KG * v2.mean() / (3.0 * K_B)

    final_start = int(0.9 * n_frames)
    T_final = float(np.mean(T_ions[final_start:]))
    T_std = float(np.std(T_ions[final_start:]))
    error_pct = abs(T_final - T_GAS_K) / T_GAS_K * 100.0
    return {"T_final": T_final, "T_std": T_std, "error_pct": error_pct, "T_ions": T_ions, "times": times}


def main() -> int:
    if not ICARION_BIN.exists():
        print(f"ERROR: ICARION binary not found: {ICARION_BIN}", file=sys.stderr)
        return 1

    print("=" * 60)
    print("Friction+OU vs HSD+OU Thermalization Comparison")
    print("=" * 60)
    print(f"  Ion:       {ION_SPECIES}  (mass {ION_MASS_AMU} amu)")
    print(f"  Gas:       He,  P={P_PA} Pa,  T={T_GAS_K} K")
    print(f"  Ions:      {N_IONS},  dt={DT_S:.0e} s,  total={TOTAL_TIME_S:.0e} s")
    print(f"  Seed:      {RNG_SEED} (fixed)")
    print(f"  Tolerance: final T within ±{TOLERANCE_PCT:.0f}% of T_gas")
    print()

    results: list[dict] = []

    for model_cfg in MODELS:
        label = model_cfg["label"]
        slug = label.replace("+", "_plus_")
        h5_name = f"ou_comparison_{slug}.h5"
        h5_path = RESULTS_DIR / h5_name
        cfg_path = RESULTS_DIR / f"ou_comparison_{slug}.json"

        cfg = _build_config(model_cfg, h5_name)
        cfg_path.write_text(json.dumps(cfg, indent=2), encoding="utf-8")

        print(f"Running {label} …")
        proc = subprocess.run(
            [str(ICARION_BIN), str(cfg_path)],
            capture_output=True,
            text=True,
            cwd=PROJECT_ROOT,
        )
        if proc.returncode != 0:
            print(f"  FAILED (exit {proc.returncode})")
            print(proc.stderr[-500:])
            results.append({"label": label, "status": "failed", "T_final": float("nan"),
                            "T_std": float("nan"), "error_pct": float("nan")})
            continue

        analysis = _analyze(h5_path)
        ok = analysis["error_pct"] <= TOLERANCE_PCT
        status = "ok" if ok else "fail"
        symbol = "✅" if ok else "❌"
        print(f"  {symbol} T_final = {analysis['T_final']:.1f} ± {analysis['T_std']:.1f} K  "
              f"(error {analysis['error_pct']:.1f}%,  tolerance ±{TOLERANCE_PCT:.0f}%)")
        results.append({
            "label": label,
            "collision_model": model_cfg["collision_model"],
            "enable_ou": model_cfg["enable_ou_thermalization"],
            "status": status,
            "T_final_K": analysis["T_final"],
            "T_std_K": analysis["T_std"],
            "T_expected_K": T_GAS_K,
            "error_pct": analysis["error_pct"],
            "tolerance_pct": TOLERANCE_PCT,
            "h5": str(h5_path),
        })

    # Write report
    report_lines = [
        "# Friction+OU vs HSD+OU Thermalization Comparison\n",
        f"- Ion: `{ION_SPECIES}` ({ION_MASS_AMU} amu)",
        f"- Gas: `He`,  P=`{P_PA} Pa`,  T=`{T_GAS_K} K`",
        f"- Ions: `{N_IONS}`,  dt=`{DT_S:.0e} s`,  total=`{TOTAL_TIME_S:.0e} s`",
        f"- Seed: `{RNG_SEED}` (fixed for reproducibility)",
        f"- Tolerance: final T within ±`{TOLERANCE_PCT:.0f}%` of T_gas\n",
        "| Model | T_final [K] | ±σ [K] | Error [%] | Tolerance [%] | Status |",
        "|---|---:|---:|---:|---:|---|",
    ]
    for r in results:
        if r["status"] == "failed":
            report_lines.append(f"| {r['label']} | — | — | — | ±{TOLERANCE_PCT:.0f}% | ❌ FAILED |")
        else:
            sym = "✅ ok" if r["status"] == "ok" else "❌ fail"
            report_lines.append(
                f"| {r['label']} | {r['T_final_K']:.1f} | {r['T_std_K']:.1f} | "
                f"{r['error_pct']:.1f} | ±{r['tolerance_pct']:.0f}% | {sym} |"
            )

    report_lines += [
        "",
        "## Interpretation",
        "",
        "Both Friction+OU and HSD+OU should thermalize H3O+ ions from a 0.1 K cold start to the gas temperature (300 K).",
        "The two models use different damping mechanisms:",
        "- **Friction+OU**: damping gamma derived from ion reduced mobility K0, OU adds thermal noise.",
        "  Thermalizes via the fluctuation-dissipation relation anchored to K0.",
        "- **HSD+OU**: damping gamma from hard-sphere collision cross-section, OU adds thermal noise.",
        "  Thermalizes via momentum-transfer rate across CCS.",
        "",
        "Discrepancies between the two final temperatures are expected at short simulation times because",
        "the thermalization time constant differs (Friction: governed by τ = m/(qγ); HSD: governed by",
        "τ = 1/(n·CCS·v_th)). Both converge to T_gas given sufficient run length.",
        "",
        f"Documented tolerance: ±{TOLERANCE_PCT:.0f}% of T_gas at t={TOTAL_TIME_S:.0e} s, P={P_PA} Pa, seed={RNG_SEED}.",
    ]

    report_path = RESULTS_DIR / "ou_comparison_report.md"
    report_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")

    summary_path = RESULTS_DIR / "ou_comparison_summary.json"
    summary_path.write_text(json.dumps(results, indent=2), encoding="utf-8")

    print()
    print(f"Report: {report_path}")
    print(f"JSON:   {summary_path}")

    failures = [r for r in results if r["status"] != "ok"]
    if failures:
        print(f"\nFailed: {[r['label'] for r in failures]}", file=sys.stderr)
        return 1
    print("\nAll models passed.")
    return 0


if __name__ == "__main__":
    sys.exit(main())

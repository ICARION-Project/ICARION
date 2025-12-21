#!/usr/bin/env python3
"""
Generate LQIT AC excitation sweep configs (pressure × frequency) and optional run script.

Default grid:
  pressures: [1e-4, 1e-3, 1e-2, 1e-1, 1] Pa
  f_AC: 50 kHz .. 80 kHz, step 250 Hz (121 points)

Example:
  python tmp/generate_lqit_ac_sweep.py \
      --output-dir tmp/lqit_ac_sweep/configs \
      --run-script tmp/run_lqit_ac_sweep.sh
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, List


DEFAULT_PRESSURES = [1e-3, 1e-2, 1e-1, 1.0, 10.0]
DEFAULT_F_START = 50_000.0
DEFAULT_F_STOP = 80_000.0
DEFAULT_F_STEP = 250.0

# AC amplitude as fraction of RF amplitude (lower-bound values from plan)
AC_SCALE_DEFAULT = {
    1e-3: 0.003,
    1e-2: 0.003,
    1e-1: 0.003,
    1.0: 0.003,
    10.0: 0.003,
}


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate LQIT AC sweep configs.")
    p.add_argument(
        "--output-dir",
        type=Path,
        default=Path("tmp/lqit_ac_sweep/configs"),
        help="Where to write generated JSON configs",
    )
    p.add_argument(
        "--run-script",
        type=Path,
        default=None,
        help="Optional path to write a helper bash script with all runs",
    )
    p.add_argument(
        "--binary",
        type=str,
        default="./build/src/icarion_main",
        help="Binary to call in the run script",
    )
    p.add_argument(
        "--pressures",
        type=float,
        nargs="+",
        default=DEFAULT_PRESSURES,
        help="Pressure grid in Pa",
    )
    p.add_argument(
        "--f-start",
        type=float,
        default=DEFAULT_F_START,
        help="Start frequency for AC sweep [Hz]",
    )
    p.add_argument(
        "--f-stop",
        type=float,
        default=DEFAULT_F_STOP,
        help="Stop frequency for AC sweep [Hz]",
    )
    p.add_argument(
        "--f-step",
        type=float,
        default=DEFAULT_F_STEP,
        help="Frequency step [Hz]",
    )
    p.add_argument(
        "--ac-scale",
        type=float,
        nargs=2,
        action="append",
        metavar=("PRESSURE_PA", "SCALE"),
        help="Override AC scale factor (V_AC = scale * V_RF) for a pressure",
    )
    p.add_argument(
        "--rf-voltage",
        type=float,
        default=50.0,
        help="RF amplitude [V]",
    )
    p.add_argument(
        "--rf-frequency",
        type=float,
        default=2e6,
        help="RF frequency [Hz]",
    )
    p.add_argument(
        "--axial-voltage",
        type=float,
        default=10.0,
        help="Axial DC voltage [V]",
    )
    p.add_argument(
        "--ion-count",
        type=int,
        default=100,
        help="Ions per species",
    )
    p.add_argument(
        "--sim-time",
        type=float,
        default=5e-3,
        help="Total simulation time [s]",
    )
    p.add_argument(
        "--dt",
        type=float,
        default=1e-8,
        help="Timestep [s]",
    )
    return p.parse_args()


def build_base_config(
    rf_voltage: float,
    rf_frequency: float,
    axial_voltage: float,
    ion_count: int,
) -> Dict:
    ion_species = []
    for mass in (80, 88, 96):
        for ccs in (25, 35, 50):
            ion_species.append(
                {
                    "species_id": f"Mass{mass}_CCS{ccs}",
                    "count": ion_count,
                    "position": {
                        "type": "gaussian",
                        "center": [0.0, 0.0, 0.0],
                        "std": [5e-4, 5e-4, 1e-3],
                    },
                    "velocity": {"type": "thermal", "temperature_K": 300.0},
                }
            )

    return {
        "title": "LQIT AC excitation sweep (auto-generated)",
        # Use absolute path to avoid nesting issues when configs live in subdirs
        "species_database": str(Path("tmp/species_database_lqit.json").resolve()),
        "simulation": {
            "total_time_s": None,  # filled later
            "dt_s": None,          # filled later
            "write_interval": 1000,
            "integrator": "RK4",
            "rng_seed": 42,
        },
        "physics": {
            "collision_model": "HSS",
            "enable_space_charge": False,
            "enable_reactions": False,
        },
        "output": {
            "folder": None,           # filled later
            "trajectory_file": None,  # filled later
            "print_progress": True,
        },
        "ions": {"species": ion_species},
        "domains": [
            {
                "name": "LQIT_ac_excitation",
                "instrument": "LQIT",
                "geometry": {
                    "origin_m": [0.0, 0.0, -0.02],
                    "length_m": 0.04,
                    "radius_m": 0.004,
                },
                "env": {
                    "temperature_K": 300.0,
                    "pressure_Pa": None,  # filled per run
                    "gas_species": "He",
                },
                "fields": {
                    "RF": {
                        "voltage_V": rf_voltage,
                        "frequency_Hz": rf_frequency,
                        "phase_rad": 0.0,
                    },
                    "DC": {
                        "axial_V": axial_voltage,
                        "quad_V": 0.0,
                        "radial_V": 0.0,
                    },
                    "AC": {
                        "voltage_V": None,      # filled per run
                        "frequency_Hz": None,   # filled per run
                        "phase_rad": 0.0,
                    },
                },
                "boundary": {"type": "Absorption"},
            }
        ],
    }


def generate_configs(
    base_cfg: Dict,
    pressures: List[float],
    freqs: List[float],
    ac_scale: Dict[float, float],
    rf_voltage: float,
    sim_time: float,
    dt: float,
    out_dir: Path,
) -> List[Path]:
    out_paths: List[Path] = []
    for p in pressures:
        scale = ac_scale.get(p, list(ac_scale.values())[-1])
        vac = scale * rf_voltage
        for f in freqs:
            cfg = json.loads(json.dumps(base_cfg))  # deep copy via JSON
            cfg["simulation"]["total_time_s"] = sim_time
            cfg["simulation"]["dt_s"] = dt
            # Absolute output folder so it's independent of config location
            cfg["output"]["folder"] = str(Path("results/lqit_ac_sweep").resolve())
            cfg["output"]["trajectory_file"] = (
                f"lqit_ac_p{p:g}Pa_f{int(f)}Hz.h5"
            )
            cfg["domains"][0]["env"]["pressure_Pa"] = p
            cfg["domains"][0]["fields"]["AC"]["voltage_V"] = vac
            cfg["domains"][0]["fields"]["AC"]["frequency_Hz"] = f

            out_dir_p = out_dir / f"p{p:g}Pa"
            out_dir_p.mkdir(parents=True, exist_ok=True)
            out_path = out_dir_p / f"f{int(f)}Hz.json"
            with out_path.open("w") as fh:
                json.dump(cfg, fh, indent=2)
            out_paths.append(out_path)
    return out_paths


def write_run_script(paths: List[Path], binary: str, script_path: Path) -> None:
    script_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
    for cfg in paths:
        lines.append(f'{binary} "{cfg}"')
    script_path.write_text("\n".join(lines) + "\n")
    script_path.chmod(0o755)
    print(f"Wrote run script with {len(paths)} calls: {script_path}")


def main() -> None:
    args = parse_args()
    ac_scale = {**AC_SCALE_DEFAULT}
    if args.ac_scale:
        for p_val, scale in args.ac_scale:
            ac_scale[float(p_val)] = float(scale)

    freqs = []
    f = args.f_start
    while f <= args.f_stop + 1e-9:
        freqs.append(round(f, 6))
        f += args.f_step

    base_cfg = build_base_config(
        rf_voltage=args.rf_voltage,
        rf_frequency=args.rf_frequency,
        axial_voltage=args.axial_voltage,
        ion_count=args.ion_count,
    )

    out_paths = generate_configs(
        base_cfg=base_cfg,
        pressures=args.pressures,
        freqs=freqs,
        ac_scale=ac_scale,
        rf_voltage=args.rf_voltage,
        sim_time=args.sim_time,
        dt=args.dt,
        out_dir=args.output_dir,
    )
    print(f"Wrote {len(out_paths)} configs to {args.output_dir}")

    if args.run_script:
        write_run_script(out_paths, args.binary, args.run_script)


if __name__ == "__main__":
    main()

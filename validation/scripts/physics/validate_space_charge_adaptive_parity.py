#!/usr/bin/env python3
"""
Adaptive SC vs. fixed-step RK4 parity test (CPU, Direct SC)

Purpose:
    Compare adaptive RK45 (stage-synchronous space charge) against a
    high-accuracy fixed-step RK4 reference for a small ion cloud using the
    Direct space-charge model. This is a quick regression (not a long bench).

Setup:
    - Single cylindrical domain, no external fields
    - Direct space charge, 200 ions
    - Reference: RK4 dt=5e-11 s, t=5e-8 s
    - Test: RK45 adaptive, dt=5e-9 s, default tolerances

Outputs:
    - Log: validation/logs/ADAPTIVE_SC_PARITY.txt
    - Summary: RMS pos/vel error and acceptance stats
"""

import json
import subprocess
import h5py
import numpy as np
from pathlib import Path
from datetime import datetime

REPO_ROOT = Path(__file__).resolve().parents[3]
ICARION_BIN = REPO_ROOT / "build" / "src" / "icarion_main"
VALIDATION_DIR = REPO_ROOT / "validation"
RESULTS_DIR = VALIDATION_DIR / "results" / "adaptive_sc_parity"
LOGS_DIR = VALIDATION_DIR / "logs"
SPECIES_DB_PATH = REPO_ROOT / "data" / "species_database_v1.json"

ION_SPECIES_ID = "Ar+"  # lightweight, deterministic
NUM_IONS = 200

TOTAL_TIME_S = 5e-8
RK4_DT_S = 5e-11
RK45_DT_S = 5e-9

def log(msg, handle):
    print(msg)
    handle.write(msg + "\n")
    handle.flush()

def config_template(integrator: str, dt_s: float, traj_name: str):
    return {
        "simulation": {
            "total_time_s": TOTAL_TIME_S,
            "dt_s": dt_s,
            "write_interval": 0,
            "integrator": integrator,
            "enable_gpu": False,
            "enable_openmp": False,
            "rng_seed": 123
        },
        "physics": {
            "collision_model": "NoCollisions",
            "enable_space_charge": True,
            "enable_space_charge_gpu": False
        },
        "output": {
            "folder": str(RESULTS_DIR),
            "trajectory_file": traj_name,
            "print_progress": False
        },
        "ions": {
            "species": [{
                "id": ION_SPECIES_ID,
                "count": NUM_IONS,
                "position": {
                    "type": "gaussian",
                    "center": [0.0, 0.0, 0.001],
                    "std": [5e-5, 5e-5, 5e-5]
                },
                "velocity": {
                    "type": "gaussian",
                    "mean": [0.0, 0.0, 0.0],
                    "std": [0.0, 0.0, 0.0]
                }
            }]
        },
        "species_database_path": str(SPECIES_DB_PATH),
        "domains": [{
            "name": "sc_parity",
            "instrument": "IMS",
            "geometry": {
                "origin_m": [0.0, 0.0, 0.0],
                "length_m": 0.01,
                "radius_m": 0.005
            },
            "env": {
                "temperature_K": 300.0,
                "pressure_Pa": 100.0,
                "gas_species": "He"
            },
            "fields": {}
        }]
    }

def run(config: dict, name: str, log_handle):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    config_path = RESULTS_DIR / f"{name}.json"
    with open(config_path, "w") as f:
        json.dump(config, f, indent=2)

    log(f"▶ Running {name} ({config['simulation']['integrator']} dt={config['simulation']['dt_s']})", log_handle)
    res = subprocess.run(
        [str(ICARION_BIN), str(config_path)],
        capture_output=True,
        text=True,
        timeout=300
    )
    if res.returncode != 0:
        log(f"  ❌ failed (rc={res.returncode}) stderr tail:\n{res.stderr[-400:]}", log_handle)
        raise RuntimeError(f"{name} failed")
    log(f"  ✓ done", log_handle)
    return RESULTS_DIR / config["output"]["trajectory_file"]

def load_state(h5_path: Path):
    with h5py.File(h5_path, "r") as f:
        pos = f["/trajectory/positions"][-1]  # last step
        vel = f["/trajectory/velocities"][-1]
    return pos, vel

def main():
    LOGS_DIR.mkdir(exist_ok=True)
    log_path = LOGS_DIR / "ADAPTIVE_SC_PARITY.txt"
    with open(log_path, "w") as log_handle:
        log(f"Adaptive SC parity run {datetime.utcnow().isoformat()}Z", log_handle)
        # Reference RK4
        cfg_rk4 = config_template("RK4", RK4_DT_S, "sc_parity_rk4.h5")
        ref_path = run(cfg_rk4, "rk4_ref", log_handle)

        # Adaptive RK45
        cfg_rk45 = config_template("RK45", RK45_DT_S, "sc_parity_rk45.h5")
        # make sure adaptive SC is enabled
        cfg_rk45["simulation"]["enable_openmp"] = False
        cfg_rk45["simulation"]["write_interval"] = 0
        rk45_path = run(cfg_rk45, "rk45_adaptive", log_handle)

        ref_pos, ref_vel = load_state(ref_path)
        rk_pos, rk_vel = load_state(rk45_path)

        pos_err = np.linalg.norm(ref_pos - rk_pos, axis=1)
        vel_err = np.linalg.norm(ref_vel - rk_vel, axis=1)

        rms_pos = float(np.sqrt(np.mean(pos_err ** 2)))
        rms_vel = float(np.sqrt(np.mean(vel_err ** 2)))

        log(f"RMS position error: {rms_pos:.3e} m", log_handle)
        log(f"RMS velocity error: {rms_vel:.3e} m/s", log_handle)

        # Simple acceptance criteria: errors should be small relative to cloud size/thermal speeds
        if rms_pos > 1e-5 or rms_vel > 1e-1:
            log(f"❌ Parity FAILED (errors too large)", log_handle)
            raise SystemExit(1)
        log(f"✅ Parity PASSED", log_handle)

if __name__ == "__main__":
    main()

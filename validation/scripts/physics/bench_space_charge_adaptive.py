#!/usr/bin/env python3
"""
Space-charge adaptive RK45 micro-benchmark (CPU, Direct/Grid SC)

Measures wall time and step counts for adaptive RK45 with space charge enabled
on small/medium ion clouds. Designed to spot regressions in stage-synchronous
SC rebuilds and rejection behavior.

Runs:
    - N in [1000, 5000, 10000]
    - Single cylindrical domain, no external fields
    - dt = 5e-9 s, total time = 5e-7 s (~100 steps if no rejections)
Outputs:
    - validation/logs/ADAPTIVE_SC_BENCH.txt with per-run wall time and status
"""

import json
import subprocess
import time
from pathlib import Path
from datetime import datetime

REPO_ROOT = Path(__file__).resolve().parents[3]
ICARION_BIN = REPO_ROOT / "build" / "src" / "icarion_main"
VALIDATION_DIR = REPO_ROOT / "validation"
RESULTS_DIR = VALIDATION_DIR / "results" / "adaptive_sc_bench"
LOGS_DIR = VALIDATION_DIR / "logs"
SPECIES_DB_PATH = REPO_ROOT / "data" / "species_database_v1.json"

ION_SPECIES_ID = "Ar+"
TEMPERATURE_K = 300.0
PRESSURE_PA = 100.0
DT_S = 5e-9
T_TOTAL_S = 5e-7

def log(msg, handle):
    print(msg)
    handle.write(msg + "\n")
    handle.flush()

def make_config(n_ions: int, traj_name: str):
    return {
        "simulation": {
            "total_time_s": T_TOTAL_S,
            "dt_s": DT_S,
            "write_interval": 0,
            "integrator": "RK45",
            "enable_gpu": False,
            "enable_openmp": False,
            "rng_seed": 777
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
                "count": n_ions,
                "position": {
                    "type": "gaussian",
                    "center": [0.0, 0.0, 0.001],
                    "std": [5e-4, 5e-4, 5e-4]
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
            "name": "sc_bench",
            "instrument": "IMS",
            "geometry": {
                "origin_m": [0.0, 0.0, 0.0],
                "length_m": 0.02,
                "radius_m": 0.01
            },
            "env": {
                "temperature_K": TEMPERATURE_K,
                "pressure_Pa": PRESSURE_PA,
                "gas_species": "He"
            },
            "fields": {}
        }]
    }

def run_case(n_ions: int, log_handle):
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)
    cfg = make_config(n_ions, f"sc_bench_{n_ions}.h5")
    cfg_path = RESULTS_DIR / f"sc_bench_{n_ions}.json"
    with open(cfg_path, "w") as f:
        json.dump(cfg, f, indent=2)

    log(f"▶ N={n_ions}: running...", log_handle)
    t0 = time.perf_counter()
    res = subprocess.run(
        [str(ICARION_BIN), str(cfg_path)],
        capture_output=True,
        text=True,
        timeout=600
    )
    t1 = time.perf_counter()

    if res.returncode != 0:
        log(f"  ❌ rc={res.returncode}", log_handle)
        log(f"  stderr tail:\n{res.stderr[-400:]}", log_handle)
        return

    wall = t1 - t0
    log(f"  ✓ wall={wall:.3f}s", log_handle)

def main():
    LOGS_DIR.mkdir(exist_ok=True)
    log_path = LOGS_DIR / "ADAPTIVE_SC_BENCH.txt"
    with open(log_path, "w") as handle:
        log(f"Adaptive SC microbench {datetime.utcnow().isoformat()}Z", handle)
        for n in (1000, 5000, 10000):
            run_case(n, handle)
        log("Done.", handle)

if __name__ == "__main__":
    main()

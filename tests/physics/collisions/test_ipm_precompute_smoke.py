#!/usr/bin/env python3
# ICARION: Ion Collision And Reaction IntegratiON
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2026 Christoph Schaefer

"""Smoke-test the IPM precompute tool and its HDF5 quality diagnostics."""

from __future__ import annotations

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import h5py
import numpy as np


def write_forcefield(path: Path) -> None:
    params = {
        "gases": {
            "He": {
                "sigma_A": 2.551,
                "epsilon_eV": 0.0008798,
                "alpha_A3": 0.204956,
            }
        },
        "elements": {
            "H": {"sigma_A": 2.261, "epsilon_eV": 0.059579817756640736},
            "O": {"sigma_A": 2.4344, "epsilon_eV": 0.1034324933362249},
        },
    }
    path.write_text(json.dumps(params), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", type=Path, required=True)
    parser.add_argument("--precompute-bin", type=Path, required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory(prefix="icarion_ipm_precompute_") as tmp:
        tmp_dir = Path(tmp)
        forcefield = tmp_dir / "forcefield.json"
        output = tmp_dir / "H3O+_ipm_samples_He.h5"
        write_forcefield(forcefield)

        cmd = [
            str(args.precompute_bin),
            "--input", str(args.repo_root / "data" / "species_database_v1.json"),
            "--species", "H3O+",
            "--gas", "He",
            "--output", str(output),
            "--gas-params", str(forcefield),
            "--element-params", str(forcefield),
            "--potential", "lj1264",
            "--polarization", "partial",
            "--mixing-rule", "lb",
            "--orient-grid", "random",
            "--n-orientations", "1",
            "--n-trials", "8",
            "--v-bins", "3",
            "--v-min", "500",
            "--v-max", "5000",
            "--max-non-asymptotic-frac", "0",
            "--max-steps", "1000000",
            "--threads", "1",
            "--seed", "12345",
            "--store-full-cdf",
        ]
        subprocess.run(cmd, cwd=args.repo_root, check=True)

        with h5py.File(output, "r") as handle:
            for dataset in (
                "sigma_event_m2",
                "sigma_mt_m2",
                "b_max_m",
                "attempted_trajectories",
                "accepted_trajectories",
                "rejected_non_asymptotic",
                "max_energy_step_error",
                "max_energy_cumulative_error",
                "cdf_values",
                "dp_samples",
            ):
                if dataset not in handle:
                    raise AssertionError(f"missing dataset: {dataset}")

            rejected = np.asarray(handle["rejected_non_asymptotic"][:], dtype=np.int64)
            if np.any(rejected != 0):
                raise AssertionError("mini precompute produced non-asymptotic rejects")

            for dataset in ("sigma_event_m2", "sigma_mt_m2", "b_max_m"):
                values = np.asarray(handle[dataset][:], dtype=float)
                if not np.all(np.isfinite(values)) or not np.all(values > 0.0):
                    raise AssertionError(f"{dataset} contains non-finite or non-positive values")

            for dataset in ("max_energy_step_error", "max_energy_cumulative_error"):
                values = np.asarray(handle[dataset][:], dtype=float)
                if not np.all(np.isfinite(values)) or np.any(values < 0.0):
                    raise AssertionError(f"{dataset} contains invalid values")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

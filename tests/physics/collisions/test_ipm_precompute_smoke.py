#!/usr/bin/env python3
# ICARION: Ion Collision And Reaction IntegratiON
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2026 Christoph Schaefer

"""Smoke-test the IPM precompute tool and its HDF5 quality diagnostics."""

from __future__ import annotations

import argparse
import hashlib
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
    parser.add_argument("--loader-bin", type=Path, required=True)
    parser.add_argument("--expected-version", required=True)
    parser.add_argument("--expected-config", required=True)
    parser.add_argument("--expected-cuda", type=int, choices=(0, 1), required=True)
    parser.add_argument("--expected-core-only", type=int, choices=(0, 1), required=True)
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
            legacy_attributes = {
                "format", "version", "species_id", "gas", "gas_model", "mixing_rule",
                "polarization", "potential", "param_model", "units", "n_orientations",
                "v_bins", "seed", "epsilon_deflection", "sigma_scale", "epsilon_scale",
                "pol_damp_A", "mmff_energy_scale", "mmff_distance_scale",
                "max_non_asymptotic_fraction", "max_steps",
            }
            missing_attributes = legacy_attributes.difference(handle.attrs)
            if missing_attributes:
                raise AssertionError(f"missing legacy attributes: {sorted(missing_attributes)}")
            format_value = handle.attrs["format"]
            if isinstance(format_value, bytes):
                format_value = format_value.decode()
            if format_value != "ipm_offline_samples":
                raise AssertionError("legacy format attribute changed")
            if int(handle.attrs["version"]) != 1:
                raise AssertionError("IPM format version changed")

            expected_groups = {"schema", "software", "system", "rng", "species", "neutral", "precompute", "inputs", "completion"}
            if expected_groups.difference(handle["metadata"]):
                raise AssertionError("incomplete /metadata hierarchy")

            def text(path: str) -> str:
                value = handle[path][()]
                return value.decode() if isinstance(value, bytes) else str(value)

            if text("metadata/species/species_id") != "H3O+" or text("metadata/neutral/gas_id") != "He":
                raise AssertionError("incorrect resolved species or gas")
            if int(handle["metadata/rng/base_seed"][()]) != 12345:
                raise AssertionError("incorrect metadata seed")
            if not bool(handle["metadata/rng/seed_supplied_explicitly"][()]):
                raise AssertionError("explicit seed not marked explicit")
            if text("metadata/schema/metadata_schema_version") != "1.0.0":
                raise AssertionError("incorrect metadata schema version")
            if int(handle["metadata/schema/data_format_version"][()]) != 1:
                raise AssertionError("incorrect metadata data-format version")
            if text("metadata/software/icarion_version") != args.expected_version:
                raise AssertionError("incorrect ICARION version")
            if not text("metadata/software/git_hash"):
                raise AssertionError("missing Git hash")
            if text("metadata/software/git_state") not in {"clean", "dirty", "unknown"}:
                raise AssertionError("invalid Git state")
            if text("metadata/software/git_state_capture") != "configure_time":
                raise AssertionError("Git state capture time is unclear")
            if text("metadata/software/build_type") != args.expected_config:
                raise AssertionError("exact CMake build configuration not stored")
            if bool(handle["metadata/software/cuda_enabled"][()]) != bool(args.expected_cuda):
                raise AssertionError("incorrect CUDA build metadata")
            expected_mode = "core-only" if args.expected_core_only else "full"
            if text("metadata/software/build_mode") != expected_mode:
                raise AssertionError("incorrect build-mode metadata")
            if "hostname" in handle["metadata/system"] or "username" in handle["metadata/system"]:
                raise AssertionError("publication-sensitive identity stored in IPM metadata")
            if float(handle["metadata/precompute/temperature_K"][()]) != 298.15:
                raise AssertionError("incorrect resolved temperature")
            if not bool(handle["metadata/completion/success"][()]) or bool(handle["metadata/completion/is_checkpoint"][()]):
                raise AssertionError("final output not marked complete")

            species_path = args.repo_root / "data" / "species_database_v1.json"
            geometry_path = args.repo_root / "data" / "molecules" / "H3O+.json"
            expected_species_hash = hashlib.sha256(species_path.read_bytes()).hexdigest()
            expected_geometry_hash = hashlib.sha256(geometry_path.read_bytes()).hexdigest()
            if text("metadata/inputs/hashes/species_database_sha256") != expected_species_hash:
                raise AssertionError("species database SHA-256 mismatch")
            if text("metadata/inputs/hashes/molecular_geometry_sha256") != expected_geometry_hash:
                raise AssertionError("geometry SHA-256 mismatch")
            if not bool(handle["metadata/inputs/hashes/molecular_geometry_used"][()]):
                raise AssertionError("geometry input not marked used")
            if not text("metadata/inputs/hashes/molecular_geometry_role"):
                raise AssertionError("geometry input role missing")
            for key in handle["metadata/inputs/hashes"]:
                if key.endswith("_filename") and Path(text(f"metadata/inputs/hashes/{key}")).is_absolute():
                    raise AssertionError("absolute input path leaked into metadata")
            if text("metadata/inputs/blobs/species_database") != species_path.read_bytes().decode("utf-8"):
                raise AssertionError("embedded species database is unreadable or changed")
            if text("metadata/inputs/blobs/molecular_geometry") != geometry_path.read_bytes().decode("utf-8"):
                raise AssertionError("embedded geometry is unreadable or changed")
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

        subprocess.run([str(args.loader_bin), str(output)], check=True)

        repeat_output = tmp_dir / "repeat_fixed_seed.h5"
        repeat_cmd = list(cmd)
        repeat_cmd[repeat_cmd.index(str(output))] = str(repeat_output)
        subprocess.run(repeat_cmd, cwd=args.repo_root, check=True)
        with h5py.File(output, "r") as first, h5py.File(repeat_output, "r") as second:
            for dataset in ("logv_bins", "orientations_quat", "sigma_event_m2", "sigma_mt_m2", "b_max_m",
                            "attempted_trajectories", "accepted_trajectories", "rejected_non_asymptotic",
                            "max_energy_step_error", "max_energy_cumulative_error", "cdf_offsets", "cdf_counts",
                            "cdf_values", "dp_samples", "dp_stats"):
                if not np.array_equal(first[dataset][...], second[dataset][...]):
                    raise AssertionError(f"fixed-seed numerical payload changed: {dataset}")

        auto_output = tmp_dir / "auto_seed.h5"
        auto_cmd = [value for value in cmd if value not in ("--seed", "12345")]
        auto_cmd[auto_cmd.index(str(output))] = str(auto_output)
        subprocess.run(auto_cmd, cwd=args.repo_root, check=True)
        with h5py.File(auto_output, "r") as handle:
            if bool(handle["metadata/rng/seed_supplied_explicitly"][()]):
                raise AssertionError("automatically generated seed marked explicit")

        zero_output = tmp_dir / "explicit_zero_seed.h5"
        zero_cmd = list(cmd)
        zero_cmd[zero_cmd.index(str(output))] = str(zero_output)
        zero_cmd[zero_cmd.index("12345")] = "0"
        subprocess.run(zero_cmd, cwd=args.repo_root, check=True)
        with h5py.File(zero_output, "r") as handle:
            if int(handle["metadata/rng/base_seed"][()]) != 0:
                raise AssertionError("explicit zero seed was replaced")
            if not bool(handle["metadata/rng/seed_supplied_explicitly"][()]):
                raise AssertionError("explicit zero seed not marked explicit")

        checkpoint_output = tmp_dir / "real_checkpoint.h5"
        checkpoint_cmd = list(cmd)
        checkpoint_cmd[checkpoint_cmd.index(str(output))] = str(checkpoint_output)
        checkpoint_cmd.remove("--store-full-cdf")
        checkpoint_cmd.extend(["--compact-dp-stats", "--checkpoint-cells", "1"])
        checkpoint_cmd[checkpoint_cmd.index("--n-orientations") + 1] = "2"
        checkpoint_cmd.append("--stop-after-checkpoint")
        checkpoint_result = subprocess.run(checkpoint_cmd, cwd=args.repo_root)
        if checkpoint_result.returncode != 2:
            raise AssertionError(f"intentional checkpoint stop returned {checkpoint_result.returncode}, expected 2")
        with h5py.File(checkpoint_output, "r") as handle:
            if bool(handle["metadata/completion/success"][()]):
                raise AssertionError("checkpoint incorrectly marked successful")
            if not bool(handle["metadata/completion/is_checkpoint"][()]):
                raise AssertionError("checkpoint not marked incomplete")
            original_start = handle["metadata/completion/start_timestamp_utc"][()]
            if isinstance(original_start, bytes):
                original_start = original_start.decode()

        changed_seed_cmd = [arg for arg in checkpoint_cmd if arg != "--stop-after-checkpoint"] + ["--resume"]
        changed_seed_cmd[changed_seed_cmd.index("12345")] = "12346"
        changed_seed = subprocess.run(changed_seed_cmd, cwd=args.repo_root, text=True, capture_output=True)
        if changed_seed.returncode == 0 or "mismatch for seed" not in changed_seed.stderr:
            raise AssertionError(f"resume did not reject changed seed: rc={changed_seed.returncode}, stderr={changed_seed.stderr!r}")

        changed_forcefield = tmp_dir / "changed_forcefield.json"
        changed_forcefield.write_bytes(forcefield.read_bytes() + b"\n")
        changed_input_cmd = [arg for arg in checkpoint_cmd if arg != "--stop-after-checkpoint"] + ["--resume"]
        changed_input_cmd[changed_input_cmd.index(str(forcefield))] = str(changed_forcefield)
        changed_input_cmd[changed_input_cmd.index(str(forcefield))] = str(changed_forcefield)
        changed_input = subprocess.run(changed_input_cmd, cwd=args.repo_root, text=True, capture_output=True)
        if changed_input.returncode == 0 or "input hash" not in changed_input.stderr:
            raise AssertionError(f"resume did not reject changed input: rc={changed_input.returncode}, stderr={changed_input.stderr!r}")

        resume_cmd = [arg for arg in checkpoint_cmd if arg != "--stop-after-checkpoint"] + ["--resume"]
        subprocess.run(resume_cmd, cwd=args.repo_root, check=True)
        with h5py.File(checkpoint_output, "r") as handle:
            resumed_start = handle["metadata/completion/start_timestamp_utc"][()]
            if isinstance(resumed_start, bytes):
                resumed_start = resumed_start.decode()
            if resumed_start != original_start:
                raise AssertionError("resume did not preserve original start timestamp")
            if not bool(handle["metadata/completion/resume_mode_used"][()]):
                raise AssertionError("resume usage not recorded")
            if not bool(handle["metadata/completion/success"][()]):
                raise AssertionError("resumed output not marked successful")
            if bool(handle["metadata/completion/is_checkpoint"][()]):
                raise AssertionError("resumed final output still marked as checkpoint")

        subprocess.run([str(args.loader_bin), str(checkpoint_output)], check=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

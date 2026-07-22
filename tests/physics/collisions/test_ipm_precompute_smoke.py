#!/usr/bin/env python3
# ICARION: Ion Collision And Reaction IntegratiON
# SPDX-License-Identifier: GPL-3.0-only
# Copyright (c) 2026 Christoph Schaefer

"""Smoke-test the IPM precompute tool and its HDF5 quality diagnostics."""

from __future__ import annotations

import argparse
import hashlib
import json
import shutil
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
            },
            "N2": {
                "sigma_A": 3.798,
                "epsilon_eV": 0.00694,
                "alpha_A3": 1.74,
                "alpha_parallel_A3": 2.45,
                "alpha_perp_A3": 1.39,
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

        scan_with_note = list(cmd)
        scan_with_note.extend(["--scan-deflection", "--scan-v-rel", "1000", "--scan-b-steps", "2",
                               "--note", "must not be discarded"])
        scan_result = subprocess.run(scan_with_note, cwd=args.repo_root, text=True, capture_output=True)
        expected_scan_diagnostic = (
            "--note and --note-file are only supported for HDF5 IPM sample or checkpoint output, "
            "not --scan-deflection."
        )
        if scan_result.returncode != 1 or expected_scan_diagnostic not in scan_result.stderr:
            raise AssertionError(
                f"scan annotation rejection failed: rc={scan_result.returncode}, stderr={scan_result.stderr!r}"
            )

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
            if "annotations" in handle["metadata"]:
                raise AssertionError("annotation group written without a note")

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
            if text("metadata/neutral/sigma_source") != "parameter_file":
                raise AssertionError("file-derived He sigma source not recorded")
            if text("metadata/neutral/epsilon_source") != "parameter_file":
                raise AssertionError("file-derived He epsilon source not recorded")
            if text("metadata/neutral/polarizability_source") != "parameter_file":
                raise AssertionError("file-derived He polarizability source not recorded")
            if not bool(handle["metadata/inputs/hashes/gas_parameter_file_used"][()]):
                raise AssertionError("used He gas parameter snapshot marked unused")
            if not bool(handle["metadata/inputs/hashes/element_parameter_file_used"][()]):
                raise AssertionError("shared element parameter input not marked used")
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

        def provenance_run(name: str, extra: list[str], gas: str = "He", params: Path = forcefield) -> Path:
            result = tmp_dir / name
            run_cmd = list(cmd)
            run_cmd[run_cmd.index(str(output))] = str(result)
            run_cmd[run_cmd.index("He")] = gas
            run_cmd[run_cmd.index(str(forcefield))] = str(params)
            run_cmd[run_cmd.index(str(forcefield))] = str(params)
            run_cmd[run_cmd.index("--n-trials") + 1] = "1"
            run_cmd[run_cmd.index("--v-bins") + 1] = "1"
            run_cmd.extend(extra)
            subprocess.run(run_cmd, cwd=args.repo_root, check=True)
            return result

        parameter_output = provenance_run("parameter_fallback.h5", [], gas="N2")
        with h5py.File(parameter_output, "r") as handle:
            def parameter_text(path: str) -> str:
                value = handle[path][()]
                return value.decode() if isinstance(value, bytes) else str(value)
            if parameter_text("metadata/neutral/sigma_source") != "parameter_file" or \
                    parameter_text("metadata/neutral/epsilon_source") != "parameter_file":
                raise AssertionError("N2 parameter-file fallback source not recorded")
            if not bool(handle["metadata/inputs/hashes/gas_parameter_file_used"][()]):
                raise AssertionError("used N2 gas parameter file marked unused")

        override_output = provenance_run(
            "cli_override.h5", ["--gas-sigma-A", "2.7", "--gas-epsilon-eV", "0.001"])
        with h5py.File(override_output, "r") as handle:
            def override_text(path: str) -> str:
                value = handle[path][()]
                return value.decode() if isinstance(value, bytes) else str(value)
            if override_text("metadata/neutral/sigma_source") != "cli_override" or \
                    override_text("metadata/neutral/epsilon_source") != "cli_override":
                raise AssertionError("CLI gas override source not recorded")
            if override_text("metadata/neutral/polarizability_source") != "parameter_file":
                raise AssertionError("file-derived He polarizability lost after LJ override")
            if not bool(handle["metadata/inputs/hashes/gas_parameter_file_used"][()]):
                raise AssertionError("He gas parameter file with used polarizability marked unused")

        he_no_alpha_forcefield = tmp_dir / "he_no_alpha_forcefield.json"
        he_no_alpha_params = json.loads(forcefield.read_text(encoding="utf-8"))
        he_no_alpha_params["gases"]["He"].pop("alpha_A3")
        he_no_alpha_forcefield.write_text(json.dumps(he_no_alpha_params), encoding="utf-8")
        he_no_alpha_output = provenance_run(
            "he_builtin_polarizability.h5", [], params=he_no_alpha_forcefield)
        with h5py.File(he_no_alpha_output, "r") as handle:
            def he_no_alpha_text(path: str) -> str:
                value = handle[path][()]
                return value.decode() if isinstance(value, bytes) else str(value)
            if he_no_alpha_text("metadata/neutral/sigma_source") != "parameter_file" or \
                    he_no_alpha_text("metadata/neutral/epsilon_source") != "parameter_file":
                raise AssertionError("He LJ parameter-file source not retained without alpha_A3")
            if he_no_alpha_text("metadata/neutral/polarizability_source") != "built_in":
                raise AssertionError("built-in He polarizability source not recorded")
            if not bool(handle["metadata/inputs/hashes/gas_parameter_file_used"][()]):
                raise AssertionError("He gas parameter file supplying LJ values marked unused")

        anisotropic_file_output = provenance_run(
            "n2_anisotropic_file.h5",
            ["--gas-model", "diatomic", "--polarization", "pairwise", "--n2-aniso-pol"],
            gas="N2")
        with h5py.File(anisotropic_file_output, "r") as handle:
            def anisotropic_text(path: str) -> str:
                value = handle[path][()]
                return value.decode() if isinstance(value, bytes) else str(value)
            if anisotropic_text("metadata/neutral/polarizability_parallel_source") != "parameter_file" or \
                    anisotropic_text("metadata/neutral/polarizability_perpendicular_source") != "parameter_file":
                raise AssertionError("file-derived anisotropic N2 source not recorded")

        isotropic_forcefield = tmp_dir / "isotropic_forcefield.json"
        isotropic_params = json.loads(forcefield.read_text(encoding="utf-8"))
        isotropic_params["gases"]["N2"].pop("alpha_A3")
        isotropic_params["gases"]["N2"].pop("alpha_parallel_A3")
        isotropic_params["gases"]["N2"].pop("alpha_perp_A3")
        isotropic_forcefield.write_text(json.dumps(isotropic_params), encoding="utf-8")
        anisotropic_builtin_output = provenance_run(
            "n2_anisotropic_builtin_fallback.h5",
            ["--gas-model", "diatomic", "--polarization", "pairwise", "--n2-aniso-pol"],
            gas="N2", params=isotropic_forcefield)
        with h5py.File(anisotropic_builtin_output, "r") as handle:
            def fallback_text(path: str) -> str:
                value = handle[path][()]
                return value.decode() if isinstance(value, bytes) else str(value)
            if fallback_text("metadata/neutral/polarizability_source") != "built_in":
                raise AssertionError("built-in isotropic N2 source not recorded")
            if fallback_text("metadata/neutral/polarizability_parallel_source") != "derived_isotropic_fallback" or \
                    fallback_text("metadata/neutral/polarizability_perpendicular_source") != "derived_isotropic_fallback":
                raise AssertionError("anisotropic N2 isotropic fallback source not recorded")

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

        annotated_output = tmp_dir / "annotated.h5"
        annotated_cmd = list(cmd)
        annotated_cmd[annotated_cmd.index(str(output))] = str(annotated_output)
        annotated_cmd.extend(["--note", "line one\nline two"])
        subprocess.run(annotated_cmd, cwd=args.repo_root, check=True)
        with h5py.File(annotated_output, "r") as annotated, h5py.File(output, "r") as plain:
            annotation = annotated["metadata/annotations"]
            if annotation["note"][()].decode() != "line one\nline two":
                raise AssertionError("inline annotation bytes changed")
            if annotation["source"][()].decode() != "inline" or annotation["source_filename"][()].decode() != "":
                raise AssertionError("inline annotation provenance incorrect")
            expected_note_hash = hashlib.sha256(b"line one\nline two").hexdigest()
            if annotation["note_sha256"][()].decode() != expected_note_hash:
                raise AssertionError("inline annotation hash incorrect")
            for dataset in ("logv_bins", "orientations_quat", "sigma_event_m2", "sigma_mt_m2", "b_max_m",
                            "cdf_offsets", "cdf_counts", "cdf_values", "dp_samples", "dp_stats"):
                if not np.array_equal(annotated[dataset][...], plain[dataset][...]):
                    raise AssertionError(f"annotation changed numerical IPM dataset: {dataset}")

        note_file = tmp_dir / "ipm_note.md"
        note_file.write_bytes(b"file note\r\nexact\n")
        file_note_output = tmp_dir / "file_note.h5"
        file_note_cmd = list(cmd)
        file_note_cmd[file_note_cmd.index(str(output))] = str(file_note_output)
        file_note_cmd.extend(["--note-file", str(note_file)])
        subprocess.run(file_note_cmd, cwd=args.repo_root, check=True)
        with h5py.File(file_note_output, "r") as handle:
            annotation = handle["metadata/annotations"]
            if annotation["note"][()].decode() != "file note\r\nexact\n":
                raise AssertionError("file annotation bytes changed")
            if annotation["source_filename"][()].decode() != "ipm_note.md":
                raise AssertionError("file annotation did not store basename only")

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
        checkpoint_cmd.extend(["--note", "checkpoint annotation"])
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
            if handle["metadata/annotations/note"][()].decode() != "checkpoint annotation":
                raise AssertionError("checkpoint annotation missing")
            original_start = handle["metadata/completion/start_timestamp_utc"][()]
            if isinstance(original_start, bytes):
                original_start = original_start.decode()

        identical_inline_checkpoint = tmp_dir / "identical_inline_checkpoint.h5"
        identical_file_checkpoint = tmp_dir / "identical_file_checkpoint.h5"
        annotation_free_checkpoint = tmp_dir / "annotation_free_checkpoint.h5"
        corrupt_annotation_checkpoint = tmp_dir / "corrupt_annotation_checkpoint.h5"
        for copy_path in (identical_inline_checkpoint, identical_file_checkpoint,
                          annotation_free_checkpoint, corrupt_annotation_checkpoint):
            shutil.copy2(checkpoint_output, copy_path)

        def checkpoint_resume_command(path: Path) -> list[str]:
            result = [arg for arg in checkpoint_cmd if arg != "--stop-after-checkpoint"] + ["--resume"]
            result[result.index(str(checkpoint_output))] = str(path)
            return result

        subprocess.run(checkpoint_resume_command(identical_inline_checkpoint), cwd=args.repo_root, check=True)
        identical_note_file = tmp_dir / "identical_checkpoint_note.txt"
        identical_note_file.write_text("checkpoint annotation", encoding="utf-8")
        identical_file_cmd = checkpoint_resume_command(identical_file_checkpoint)
        note_index = identical_file_cmd.index("--note")
        identical_file_cmd[note_index:note_index + 2] = ["--note-file", str(identical_note_file)]
        subprocess.run(identical_file_cmd, cwd=args.repo_root, check=True)
        with h5py.File(identical_file_checkpoint, "r") as handle:
            if handle["metadata/annotations/source"][()].decode() != "inline":
                raise AssertionError("identical file resume rewrote original annotation provenance")

        with h5py.File(annotation_free_checkpoint, "r+") as handle:
            del handle["metadata/annotations"]
        added_annotation = subprocess.run(checkpoint_resume_command(annotation_free_checkpoint),
                                          cwd=args.repo_root, text=True, capture_output=True)
        if added_annotation.returncode == 0 or "checkpoint has no annotation" not in added_annotation.stderr:
            raise AssertionError("resume added an annotation to an annotation-free checkpoint")

        with h5py.File(corrupt_annotation_checkpoint, "r+") as handle:
            handle["metadata/annotations/note_sha256"][()] = "0" * 64
        corrupt_annotation = subprocess.run(checkpoint_resume_command(corrupt_annotation_checkpoint),
                                            cwd=args.repo_root, text=True, capture_output=True)
        if corrupt_annotation.returncode != 1 or "annotation integrity" not in corrupt_annotation.stderr:
            raise AssertionError("resume accepted corrupted annotation metadata")

        corruption_cases = {
            "missing_dataset": lambda group: group.__delitem__("note"),
            "numeric_note": lambda group: (
                group.__delitem__("note"), group.create_dataset("note", data=np.int32(7))
            ),
            "nonscalar_note": lambda group: (
                group.__delitem__("note"),
                group.create_dataset("note", data=np.array(["checkpoint annotation"],
                                     dtype=h5py.string_dtype(encoding="utf-8")))
            ),
            "fixed_note": lambda group: (
                group.__delitem__("note"), group.create_dataset("note", data=np.bytes_("checkpoint annotation"))
            ),
            "invalid_utf8_note": lambda group: (
                group.__delitem__("note"),
                group.create_dataset("note", data=b"bad\xff", dtype=h5py.string_dtype(encoding="utf-8"))
            ),
        }
        for label, mutate in corruption_cases.items():
            corrupt_path = tmp_dir / f"corrupt_{label}.h5"
            shutil.copy2(checkpoint_output, corrupt_path)
            with h5py.File(corrupt_path, "r+") as handle:
                mutate(handle["metadata/annotations"])
            result = subprocess.run(checkpoint_resume_command(corrupt_path), cwd=args.repo_root,
                                    text=True, capture_output=True)
            if result.returncode != 1 or "annotation integrity" not in result.stderr:
                raise AssertionError(
                    f"resume corruption case {label} was not rejected cleanly: "
                    f"rc={result.returncode}, stderr={result.stderr!r}"
                )

        changed_seed_cmd = [arg for arg in checkpoint_cmd if arg != "--stop-after-checkpoint"] + ["--resume"]
        changed_seed_cmd[changed_seed_cmd.index("12345")] = "12346"
        changed_seed = subprocess.run(changed_seed_cmd, cwd=args.repo_root, text=True, capture_output=True)
        if changed_seed.returncode == 0 or "mismatch for seed" not in changed_seed.stderr:
            raise AssertionError(f"resume did not reject changed seed: rc={changed_seed.returncode}, stderr={changed_seed.stderr!r}")

        changed_annotation_cmd = [arg for arg in checkpoint_cmd if arg != "--stop-after-checkpoint"] + ["--resume"]
        changed_annotation_cmd[changed_annotation_cmd.index("checkpoint annotation")] = "different annotation"
        changed_annotation = subprocess.run(changed_annotation_cmd, cwd=args.repo_root, text=True, capture_output=True)
        if changed_annotation.returncode == 0 or "annotation integrity" not in changed_annotation.stderr:
            raise AssertionError("resume accepted a different annotation")

        changed_forcefield = tmp_dir / "changed_forcefield.json"
        changed_forcefield.write_bytes(forcefield.read_bytes() + b"\n")
        changed_input_cmd = [arg for arg in checkpoint_cmd if arg != "--stop-after-checkpoint"] + ["--resume"]
        changed_input_cmd[changed_input_cmd.index(str(forcefield))] = str(changed_forcefield)
        changed_input_cmd[changed_input_cmd.index(str(forcefield))] = str(changed_forcefield)
        changed_input = subprocess.run(changed_input_cmd, cwd=args.repo_root, text=True, capture_output=True)
        if changed_input.returncode == 0 or "input hash" not in changed_input.stderr:
            raise AssertionError(f"resume did not reject changed input: rc={changed_input.returncode}, stderr={changed_input.stderr!r}")

        resume_cmd = [arg for arg in checkpoint_cmd if arg not in ("--stop-after-checkpoint", "--note", "checkpoint annotation")] + ["--resume"]
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
            if handle["metadata/annotations/note"][()].decode() != "checkpoint annotation":
                raise AssertionError("resume without flags did not preserve annotation")
            if handle["metadata/annotations/source"][()].decode() != "inline":
                raise AssertionError("resume rewrote annotation provenance")

        subprocess.run([str(args.loader_bin), str(checkpoint_output)], check=True)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

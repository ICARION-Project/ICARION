#!/usr/bin/env python3
"""End-to-end tests for user annotations written by the normal ICARION CLI."""

import argparse
import json
import subprocess
import tempfile
from pathlib import Path

import h5py


def decode(dataset) -> str:
    value = dataset[()]
    return value.decode("utf-8") if isinstance(value, bytes) else str(value)


def run(command: list[str], *, expect_success: bool = True) -> subprocess.CompletedProcess[str]:
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    if expect_success and result.returncode != 0:
        raise AssertionError(f"command failed ({result.returncode}): {' '.join(command)}\n{result.stdout}\n{result.stderr}")
    if not expect_success and result.returncode == 0:
        raise AssertionError(f"command unexpectedly succeeded: {' '.join(command)}")
    return result


def assert_annotation(path: Path, note: str, source: str, filename: str) -> None:
    with h5py.File(path, "r") as handle:
        group = handle["metadata/annotations"]
        assert decode(group["note"]) == note
        assert decode(group["source"]) == source
        assert decode(group["source_filename"]) == filename


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo-root", required=True, type=Path)
    parser.add_argument("--icarion-bin", required=True, type=Path)
    args = parser.parse_args()
    repo_root = args.repo_root.resolve()

    with tempfile.TemporaryDirectory(prefix="icarion_annotation_cli_") as directory:
        temporary = Path(directory)
        with (repo_root / "examples/ims/ims_output_diagnostics_minimal.json").open(encoding="utf-8") as stream:
            config = json.load(stream)
        config["simulation"].update({"total_time_s": 4e-9, "dt_s": 4e-9, "write_interval": 1,
                                     "enable_openmp": False})
        config["physics"]["collision_model"] = "NoCollisions"
        config["ions"]["species"][0]["count"] = 1
        config["species_database"] = str((repo_root / "data/species_database_v1.json").resolve())
        config["output"].update({"folder": str(temporary), "print_progress": False})
        config_path = temporary / "config.json"
        config_path.write_text(json.dumps(config), encoding="utf-8")

        def command(output: str) -> list[str]:
            return [str(args.icarion_bin), str(config_path), "--output", output,
                    "--log-level", "ERROR"]

        inline_path = temporary / "inline.h5"
        run(command(inline_path.name) + ["--note", "CLI Grüße 😀"])
        assert_annotation(inline_path, "CLI Grüße 😀", "inline", "")

        note_bytes = " exact file bytes \r\n東京\n".encode("utf-8")
        note_file = temporary / "context.md"
        note_file.write_bytes(note_bytes)
        file_path = temporary / "file.h5"
        run(command(file_path.name) + ["--note-file", str(note_file)])
        assert_annotation(file_path, note_bytes.decode("utf-8"), "file", note_file.name)

        config["metadata"] = {"note": "configured note"}
        config_path.write_text(json.dumps(config), encoding="utf-8")
        override_path = temporary / "override.h5"
        run(command(override_path.name) + ["--note", "CLI override"])
        assert_annotation(override_path, "CLI override", "inline", "")

        conflict = run(command("conflict.h5") + ["--note", "one", "--note-file", str(note_file)],
                       expect_success=False)
        if "mutually exclusive" not in conflict.stderr:
            raise AssertionError(f"missing mutual-exclusion diagnostic: {conflict.stderr}")

        del config["metadata"]
        config_path.write_text(json.dumps(config), encoding="utf-8")
        plain_path = temporary / "plain.h5"
        run(command(plain_path.name))
        with h5py.File(plain_path, "r") as handle:
            if "annotations" in handle["metadata"]:
                raise AssertionError("unannotated simulation unexpectedly contains metadata/annotations")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""Create a compact run summary for ICARION HDF5 output.

This script reads one trajectory or minimal-output HDF5 file and writes a
small table with species counts, time span, trajectory shape, domains, and
selected optional diagnostics. It is intended as a quick inspection tool
before running more specialized analysis scripts.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np

if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.trajectory import decode_hdf5_string, read_domain_environment, read_species_ids


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Summarize species, domains, trajectory shape, and optional diagnostics.")
    p.add_argument("--traj", "--input", dest="traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", "--output", dest="out", type=Path, default=None, help="Optional CSV summary path.")
    return p.parse_args()


def _domain_rows(h5: h5py.File) -> list[dict[str, object]]:
    domains = h5.get("domains")
    if domains is None:
        return []
    rows = []
    for name in sorted(domains.keys()):
        idx = int(name.split("_")[-1]) if name.startswith("domain_") and name.split("_")[-1].isdigit() else -1
        inst = domains[name].attrs.get("instrument", domains[name].attrs.get("type", None))
        if inst is None and "instrument" in domains[name]:
            inst = domains[name]["instrument"][()]
        if inst is None:
            inst = "unknown"
        temp, pressure = read_domain_environment(h5, idx, require_pressure=False, allow_missing=True) if idx >= 0 else (None, None)
        rows.append(
            {
                "domain": name,
                "instrument": decode_hdf5_string(inst),
                "temperature_K": temp,
                "pressure_Pa": pressure,
            }
        )
    return rows


def build_report(path: Path) -> dict[str, object]:
    with h5py.File(path, "r") as h5:
        if "trajectory/time" in h5 and "trajectory/species_ids" in h5:
            traj = h5["trajectory"]
            time = np.asarray(traj["time"], dtype=float)
            species = read_species_ids(traj)
            output_mode = "trajectory"
        elif "analysis/minimal_transport" in h5:
            mt = h5["analysis/minimal_transport"]
            time = np.asarray(mt["ion_time_s"], dtype=float) if "ion_time_s" in mt else np.array([], dtype=float)
            species = _minimal_transport_species(mt)
            output_mode = "minimal_transport"
        else:
            raise KeyError(f"{path} missing '/trajectory' and '/analysis/minimal_transport'")
        unique, counts = np.unique(species, return_counts=True)
        death_times = np.asarray(h5["ions/death_time_s"], dtype=float) if "ions/death_time_s" in h5 else None
        active = int(np.count_nonzero((death_times <= 0.0) | ~np.isfinite(death_times))) if death_times is not None else None
        optional = {
            "minimal_transport": "analysis/minimal_transport" in h5,
            "deep_collision": "analysis/deep_collision" in h5,
        }
        return {
            "file": str(path),
            "output_mode": output_mode,
            "n_frames": int(time.size),
            "n_ions": int(species.size),
            "t_start_s": float(time[0]) if time.size else float("nan"),
            "t_end_s": float(time[-1]) if time.size else float("nan"),
            "species_counts": dict(zip(unique.tolist(), counts.astype(int).tolist())),
            "active_ions_final": active,
            "domains": _domain_rows(h5),
            "optional_diagnostics": optional,
        }


def _minimal_transport_species(mt: h5py.Group) -> np.ndarray:
    if "species_id_indices" not in mt:
        return np.array([], dtype=str)
    indices = np.asarray(mt["species_id_indices"], dtype=int)
    if "species_pool" not in mt:
        return np.array([str(i) for i in indices])
    pool_raw = np.asarray(mt["species_pool"])
    pool = np.array([decode_hdf5_string(v) for v in pool_raw.reshape(-1)])
    clipped = np.clip(indices, 0, max(0, pool.size - 1))
    return pool[clipped]


def print_report(report: dict[str, object]) -> None:
    print(f"file: {report['file']}")
    print(f"output_mode: {report['output_mode']}")
    print(f"frames: {report['n_frames']}  ions: {report['n_ions']}  t_end_s: {report['t_end_s']:.6g}")
    if report["active_ions_final"] is not None:
        print(f"active_ions_final: {report['active_ions_final']}")
    print("species:")
    for species, count in report["species_counts"].items():
        print(f"  {species}: {count}")
    print("domains:")
    for row in report["domains"]:
        print(
            f"  {row['domain']}: instrument={row['instrument']} "
            f"T={row['temperature_K']} K p={row['pressure_Pa']} Pa"
        )
    print("optional_diagnostics:")
    for key, present in report["optional_diagnostics"].items():
        print(f"  {key}: {present}")


def write_csv(report: dict[str, object], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for key in ("file", "output_mode", "n_frames", "n_ions", "t_start_s", "t_end_s", "active_ions_final"):
            writer.writerow([key, report[key]])
        for species, count in report["species_counts"].items():
            writer.writerow([f"species_count.{species}", count])
        for row in report["domains"]:
            writer.writerow([f"domain.{row['domain']}.instrument", row["instrument"]])
            writer.writerow([f"domain.{row['domain']}.temperature_K", row["temperature_K"]])
            writer.writerow([f"domain.{row['domain']}.pressure_Pa", row["pressure_Pa"]])
        for key, present in report["optional_diagnostics"].items():
            writer.writerow([f"optional.{key}", present])


def main() -> int:
    args = parse_args()
    report = build_report(args.traj)
    print_report(report)
    if args.out is not None:
        write_csv(report, args.out)
        print(f"Wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

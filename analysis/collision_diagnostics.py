#!/usr/bin/env python3
"""Summarize optional /analysis/deep_collision HDF5 diagnostics."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import h5py
import numpy as np

if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Summarize deep-collision diagnostics if present in an ICARION HDF5 file.")
    p.add_argument("--traj", "--input", dest="traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out-csv", "--output", dest="out_csv", type=Path, default=Path("analysis/output/collision_diagnostics.csv"), help="Output CSV path.")
    return p.parse_args()


def _visit_datasets(group: h5py.Group, prefix: str = "") -> list[tuple[str, np.ndarray]]:
    rows = []
    for key, item in group.items():
        path = f"{prefix}/{key}" if prefix else key
        if isinstance(item, h5py.Dataset):
            rows.append((path, np.asarray(item)))
        elif isinstance(item, h5py.Group):
            rows.extend(_visit_datasets(item, path))
    return rows


def main() -> int:
    args = parse_args()
    with h5py.File(args.traj, "r") as h5:
        if "analysis/deep_collision" not in h5:
            raise KeyError("Missing /analysis/deep_collision diagnostics in input file.")
        datasets = _visit_datasets(h5["analysis/deep_collision"])

    rows = []
    for name, arr in datasets:
        numeric = np.asarray(arr, dtype=float) if np.issubdtype(arr.dtype, np.number) else None
        rows.append(
            {
                "dataset": name,
                "shape": "x".join(str(s) for s in arr.shape),
                "count": int(arr.size),
                "mean": "" if numeric is None or numeric.size == 0 else float(np.nanmean(numeric)),
                "min": "" if numeric is None or numeric.size == 0 else float(np.nanmin(numeric)),
                "max": "" if numeric is None or numeric.size == 0 else float(np.nanmax(numeric)),
            }
        )

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=["dataset", "shape", "count", "mean", "min", "max"])
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

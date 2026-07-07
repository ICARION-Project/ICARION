#!/usr/bin/env python3
"""Summarize species population kinetics from trajectory species IDs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.trajectory import decode_hdf5_string, open_trajectory


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot species population versus time for reaction-enabled trajectories.")
    p.add_argument("--traj", "--input", dest="traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", "--output", dest="out", type=Path, default=Path("analysis/output/reaction_kinetics.png"), help="Output plot path.")
    p.add_argument("--out-csv", type=Path, default=Path("analysis/output/reaction_kinetics.csv"), help="Output CSV path.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    with open_trajectory(args.traj) as h5:
        traj = h5["trajectory"]
        time = np.asarray(traj["time"], dtype=float)
        species_ds = traj["species_ids"]
        if species_ds.ndim < 2:
            raise ValueError("reaction_kinetics requires time-resolved /trajectory/species_ids with shape (time, ion)")
        species = np.asarray([[decode_hdf5_string(v) for v in row] for row in species_ds[:]])

    labels = sorted(set(species.reshape(-1).tolist()))
    counts = {label: np.count_nonzero(species == label, axis=1) for label in labels}

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for label in labels:
        ax.plot(time, counts[label], label=label)
    ax.set_xlabel("time [s]")
    ax.set_ylabel("ion count")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(args.out, dpi=180)
    plt.close(fig)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", *labels])
        for i, t in enumerate(time):
            writer.writerow([t, *[counts[label][i] for label in labels]])
    print(f"Wrote {args.out}")
    print(f"Wrote {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

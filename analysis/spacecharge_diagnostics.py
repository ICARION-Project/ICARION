#!/usr/bin/env python3
"""Cloud size diagnostics useful for space charge broadening studies."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.trajectory import (
    ensure_nonempty_selection,
    normalize_species_filter,
    open_trajectory,
    read_positions_subset,
    select_ions_from_trajectory,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute trajectory cloud centroid and RMS width versus time.")
    p.add_argument("--traj", "--input", dest="traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", "--output", dest="out", type=Path, default=Path("analysis/output/spacecharge_diagnostics.png"), help="Output plot path.")
    p.add_argument("--out-csv", type=Path, default=Path("analysis/output/spacecharge_diagnostics.csv"), help="Output CSV path.")
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--max-ions", type=int, default=5000, help="Maximum ions to analyze.")
    p.add_argument("--time-stride", type=int, default=1, help="Time stride.")
    p.add_argument("--rng-seed", type=int, default=0, help="Subsampling seed.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    with open_trajectory(args.traj) as h5:
        traj = h5["trajectory"]
        species_ids, ion_indices = select_ions_from_trajectory(
            traj,
            species_filter=normalize_species_filter(args.species),
            max_ions=args.max_ions,
            max_per_species=None,
            rng_seed=args.rng_seed,
        )
        ensure_nonempty_selection(ion_indices)
        time, positions = read_positions_subset(traj, ion_indices, args.time_stride, max_frames=None)

    centroid = np.nanmean(positions, axis=1)
    centered = positions - centroid[:, None, :]
    rms_x = np.sqrt(np.nanmean(centered[:, :, 0] ** 2, axis=1))
    rms_y = np.sqrt(np.nanmean(centered[:, :, 1] ** 2, axis=1))
    rms_z = np.sqrt(np.nanmean(centered[:, :, 2] ** 2, axis=1))
    rms_r = np.sqrt(rms_x**2 + rms_y**2)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(time, rms_r, label="rms radial")
    ax.plot(time, rms_z, label="rms axial")
    ax.set_xlabel("time [s]")
    ax.set_ylabel("RMS width [m]")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(args.out, dpi=180)
    plt.close(fig)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["time_s", "centroid_x_m", "centroid_y_m", "centroid_z_m", "rms_x_m", "rms_y_m", "rms_z_m", "rms_radial_m"])
        for i, t in enumerate(time):
            writer.writerow([t, centroid[i, 0], centroid[i, 1], centroid[i, 2], rms_x[i], rms_y[i], rms_z[i], rms_r[i]])
    print(f"Wrote {args.out}")
    print(f"Wrote {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

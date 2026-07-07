#!/usr/bin/env python3
"""Summarize trap transmission, survival, and radial envelopes."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.stability import summarize_stability, survival_fraction_series
from analysis.trajectory import (
    ensure_nonempty_selection,
    normalize_species_filter,
    open_trajectory,
    read_domain_bounds,
    read_positions_subset,
    select_ions_from_trajectory,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Analyze transmission and radial stability of trap/filter trajectories.")
    p.add_argument("--traj", "--input", dest="traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", "--output", dest="out", type=Path, default=Path("analysis/output/trap_stability.png"), help="Output plot path.")
    p.add_argument("--out-csv", type=Path, default=Path("analysis/output/trap_stability.csv"), help="Output CSV path.")
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--domain-index", type=int, default=0, help="Domain index used for geometry defaults.")
    p.add_argument("--radius", type=float, default=None, help="Override radial aperture [m].")
    p.add_argument("--z-min", type=float, default=None, help="Override transmitted z lower bound [m].")
    p.add_argument("--z-max", type=float, default=None, help="Override transmitted z upper bound [m].")
    p.add_argument("--max-ions", type=int, default=2000, help="Maximum ions to analyze.")
    p.add_argument("--time-stride", type=int, default=1, help="Time stride for envelope analysis.")
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
        bounds = read_domain_bounds(h5, args.domain_index, np.asarray(traj["positions"][-1], dtype=float))

    radius = args.radius if args.radius is not None else bounds.r_out
    z_min = args.z_min if args.z_min is not None else bounds.z_min
    z_max = args.z_max if args.z_max is not None else bounds.z_max
    summary = summarize_stability(positions, radius_m=radius, z_min_m=z_min, z_max_m=z_max)
    survival = survival_fraction_series(positions, radius_m=radius)
    r = np.sqrt(positions[:, :, 0] ** 2 + positions[:, :, 1] ** 2)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax1 = plt.subplots(figsize=(8, 4.8))
    ax1.plot(time, np.nanmean(r, axis=1), label="mean radius")
    ax1.plot(time, np.nanmax(r, axis=1), label="max radius")
    if radius is not None:
        ax1.axhline(radius, color="black", linestyle="--", linewidth=1, label="aperture")
    ax1.set_xlabel("time [s]")
    ax1.set_ylabel("radius [m]")
    ax2 = ax1.twinx()
    ax2.plot(time, survival, color="tab:green", label="survival")
    ax2.set_ylabel("survival fraction")
    ax1.grid(alpha=0.25)
    lines, labels = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines + lines2, labels + labels2, loc="best")
    fig.tight_layout()
    fig.savefig(args.out, dpi=180)
    plt.close(fig)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["metric", "value"])
        for key, value in summary.__dict__.items():
            writer.writerow([key, value])
        writer.writerow(["radius_m", radius])
        writer.writerow(["z_min_m", z_min])
        writer.writerow(["z_max_m", z_max])
    print(f"transmission_fraction: {summary.transmission_fraction:.6g}")
    print(f"Wrote {args.out}")
    print(f"Wrote {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

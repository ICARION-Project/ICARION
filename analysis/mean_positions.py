#!/usr/bin/env python3
"""
Plottet mittlere x/y/z-Positionen über die Zeit aus ICARION HDF5-Trajektorien.
Optional pro Species gefiltert und mit Ion-Downsampling/Time-Striding für große Dateien.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Allow running as a standalone script without installing the package
if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.common import (
    active_ion_mask,
    normalize_species_filter,
    read_death_times_for_indices,
    read_trajectory_selection,
    species_color_map,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Mean x/y/z positions over time.")
    p.add_argument("--traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", type=Path, default=Path("analysis/output/mean_positions.png"), help="Output PNG path.")
    p.add_argument("--csv", type=Path, default=None, help="Optional CSV output path for mean coordinates.")
    p.add_argument(
        "--active-only",
        action="store_true",
        help="Average only over ions that are active at each timepoint (requires ions/death_time_s).",
    )
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--max-ions", type=int, default=400, help="Cap on ions after filtering.")
    p.add_argument("--max-per-species", type=int, default=200, help="Cap per species before global cap.")
    p.add_argument("--rng-seed", type=int, default=0, help="Seed for ion subsampling.")
    p.add_argument("--time-stride", type=int, default=5, help="Read every Nth timestep.")
    p.add_argument("--max-frames", type=int, default=None, help="Optional frame cap after striding.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    species_filter = normalize_species_filter(args.species)
    selection = read_trajectory_selection(
        traj_path=args.traj,
        species_filter=species_filter,
        max_ions=args.max_ions,
        max_per_species=args.max_per_species,
        rng_seed=args.rng_seed,
        time_stride=args.time_stride,
        max_frames=args.max_frames,
    )

    time = selection.time
    positions = selection.positions
    ion_indices = selection.ion_indices
    ion_species = selection.species_ids_selected

    death_times = None
    if args.active_only:
        death_times = read_death_times_for_indices(args.traj, ion_indices)
    unique_species = sorted(set(ion_species))
    if not unique_species:
        raise RuntimeError("No species available after filtering.")

    active_mask = None
    if args.active_only:
        if death_times is None:
            raise RuntimeError("--active-only requested but ions/death_time_s not found in HDF5")
        active_mask = active_ion_mask(time, death_times)

    means_x, means_y, means_z = compute_stats(positions, ion_species, unique_species, active_mask)
    plot_means(time, means_x, means_y, means_z, args.out)
    if args.csv is not None:
        write_csv(time, means_x, means_y, means_z, args.csv)
        print(f"Wrote {args.csv}")
    print(f"Wrote {args.out}")
    return 0


def compute_stats(
    positions: np.ndarray,
    ion_species: np.ndarray,
    unique_species: list[str],
    active_mask: np.ndarray | None = None,
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
    x = positions[:, :, 0]
    y = positions[:, :, 1]
    z = positions[:, :, 2]

    means_x = {}
    means_y = {}
    means_z = {}
    for sp in unique_species:
        mask = ion_species == sp
        if not np.any(mask):
            continue
        subset_x = x[:, mask]
        subset_y = y[:, mask]
        subset_z = z[:, mask]

        if active_mask is None:
            means_x[sp] = subset_x.mean(axis=1)
            means_y[sp] = subset_y.mean(axis=1)
            means_z[sp] = subset_z.mean(axis=1)
        else:
            subset_active = active_mask[:, mask]
            count = subset_active.sum(axis=1).astype(float)

            sum_x = (subset_x * subset_active).sum(axis=1)
            sum_y = (subset_y * subset_active).sum(axis=1)
            sum_z = (subset_z * subset_active).sum(axis=1)

            mean_x = np.full_like(sum_x, np.nan, dtype=float)
            mean_y = np.full_like(sum_y, np.nan, dtype=float)
            mean_z = np.full_like(sum_z, np.nan, dtype=float)

            valid = count > 0.0
            mean_x[valid] = sum_x[valid] / count[valid]
            mean_y[valid] = sum_y[valid] / count[valid]
            mean_z[valid] = sum_z[valid] / count[valid]

            means_x[sp] = mean_x
            means_y[sp] = mean_y
            means_z[sp] = mean_z
    return means_x, means_y, means_z


def plot_means(
    time: np.ndarray,
    means_x: dict[str, np.ndarray],
    means_y: dict[str, np.ndarray],
    means_z: dict[str, np.ndarray],
    out_path: Path,
):
    colors = species_color_map(means_x.keys())
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.5))

    ax_x, ax_y, ax_z = axes

    for sp, vals in means_x.items():
        ax_x.plot(time, vals, label=sp, color=colors[sp])
    ax_x.set_xlabel("Time [s]")
    ax_x.set_ylabel("Mean x [m]")
    ax_x.set_title("Mean x vs time")
    ax_x.grid(True, alpha=0.3)

    for sp, vals in means_y.items():
        ax_y.plot(time, vals, label=sp, color=colors[sp])
    ax_y.set_xlabel("Time [s]")
    ax_y.set_ylabel("Mean y [m]")
    ax_y.set_title("Mean y vs time")
    ax_y.grid(True, alpha=0.3)

    for sp, vals in means_z.items():
        ax_z.plot(time, vals, label=sp, color=colors[sp])
    ax_z.set_xlabel("Time [s]")
    ax_z.set_ylabel("Mean z [m]")
    ax_z.set_title("Mean z vs time")
    ax_z.grid(True, alpha=0.3)

    handles, labels = ax_x.get_legend_handles_labels()
    if labels:
        fig.legend(handles, labels, loc="upper center", ncol=min(4, len(labels)), bbox_to_anchor=(0.5, 1.05))

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)


def write_csv(
    time: np.ndarray,
    means_x: dict[str, np.ndarray],
    means_y: dict[str, np.ndarray],
    means_z: dict[str, np.ndarray],
    out_path: Path,
) -> None:
    species = sorted(means_x.keys())
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        header = ["time_s"]
        for sp in species:
            header += [f"{sp}_mean_x_m", f"{sp}_mean_y_m", f"{sp}_mean_z_m"]
        writer.writerow(header)

        for i, ti in enumerate(time):
            row = [f"{ti:.12e}"]
            for sp in species:
                row += [
                    f"{means_x[sp][i]:.12e}",
                    f"{means_y[sp][i]:.12e}",
                    f"{means_z[sp][i]:.12e}",
                ]
            writer.writerow(row)


if __name__ == "__main__":
    raise SystemExit(main())

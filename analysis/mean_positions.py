#!/usr/bin/env python3
"""
Plottet mittlere radiale und axiale Positionen über die Zeit aus ICARION HDF5-Trajektorien.
Optional pro Species gefiltert und mit Ion-Downsampling/Time-Striding für große Dateien.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

# Allow running as a standalone script without installing the package
if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.common import (
    load_positions_subset,
    load_species_ids,
    open_trajectory,
    select_ion_indices,
    species_color_map,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Mean radial and axial positions over time.")
    p.add_argument("--traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", type=Path, default=Path("analysis/output/mean_positions.png"), help="Output PNG path.")
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--max-ions", type=int, default=400, help="Cap on ions after filtering.")
    p.add_argument("--max-per-species", type=int, default=200, help="Cap per species before global cap.")
    p.add_argument("--rng-seed", type=int, default=0, help="Seed for ion subsampling.")
    p.add_argument("--time-stride", type=int, default=5, help="Read every Nth timestep.")
    p.add_argument("--max-frames", type=int, default=None, help="Optional frame cap after striding.")
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.traj.exists():
        raise FileNotFoundError(f"Trajectory file not found: {args.traj}")

    with open_trajectory(args.traj) as h5:
        traj = h5["trajectory"]
        species_ids = load_species_ids(traj)
        ion_indices = select_ion_indices(
            species_ids,
            species_filter=set(args.species) if args.species else None,
            max_ions=args.max_ions,
            max_per_species=args.max_per_species,
            rng_seed=args.rng_seed,
        )
        if len(ion_indices) == 0:
            raise RuntimeError("No ions selected; adjust species filter or caps.")

        time, positions = load_positions_subset(
            traj, ion_indices=ion_indices, time_stride=args.time_stride, max_frames=args.max_frames
        )

    ion_species = np.array(species_ids)[ion_indices]
    unique_species = sorted(set(ion_species))
    if not unique_species:
        raise RuntimeError("No species available after filtering.")

    means_r, stds_r, means_z, stds_z = compute_stats(positions, ion_species, unique_species)
    plot_means(time, means_r, stds_r, means_z, stds_z, args.out)
    print(f"Wrote {args.out}")
    return 0


def compute_stats(
    positions: np.ndarray, ion_species: np.ndarray, unique_species: list[str]
) -> tuple[dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray], dict[str, np.ndarray]]:
    r = np.sqrt(positions[:, :, 0] ** 2 + positions[:, :, 1] ** 2)
    z = positions[:, :, 2]

    means_r = {}
    stds_r = {}
    means_z = {}
    stds_z = {}
    for sp in unique_species:
        mask = ion_species == sp
        if not np.any(mask):
            continue
        subset_r = r[:, mask]
        subset_z = z[:, mask]
        means_r[sp] = subset_r.mean(axis=1)
        stds_r[sp] = subset_r.std(axis=1)
        means_z[sp] = subset_z.mean(axis=1)
        stds_z[sp] = subset_z.std(axis=1)
    return means_r, stds_r, means_z, stds_z


def plot_means(
    time: np.ndarray,
    means_r: dict[str, np.ndarray],
    stds_r: dict[str, np.ndarray],
    means_z: dict[str, np.ndarray],
    stds_z: dict[str, np.ndarray],
    out_path: Path,
):
    colors = species_color_map(means_r.keys())
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    ax_r, ax_z = axes
    for sp, vals in means_r.items():
        ax_r.plot(time, vals, label=sp, color=colors[sp])
        ax_r.fill_between(time, vals - stds_r[sp], vals + stds_r[sp], color=colors[sp], alpha=0.15)
    ax_r.set_xlabel("Time [s]")
    ax_r.set_ylabel("Mean radial position [m]")
    ax_r.set_title("Mean radial vs time")
    ax_r.grid(True, alpha=0.3)

    for sp, vals in means_z.items():
        ax_z.plot(time, vals, label=sp, color=colors[sp])
        ax_z.fill_between(time, vals - stds_z[sp], vals + stds_z[sp], color=colors[sp], alpha=0.15)
    ax_z.set_xlabel("Time [s]")
    ax_z.set_ylabel("Mean axial position [m]")
    ax_z.set_title("Mean axial vs time")
    ax_z.grid(True, alpha=0.3)

    handles, labels = ax_r.get_legend_handles_labels()
    if labels:
        fig.legend(handles, labels, loc="upper center", ncol=min(4, len(labels)), bbox_to_anchor=(0.5, 1.05))

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Plot a handful of trajectories from an ICARION HDF5 file.
Default output shows XY and XZ projections plus a simple 3D view.
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
    p = argparse.ArgumentParser(description="Plot a subset of trajectories (XY/XZ/3D).")
    p.add_argument(
        "--traj",
        type=Path,
        default=Path("results/lqit/lqit_trajectories.h5"),
        help="Path to trajectory HDF5 file.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("analysis/output/trajectory_plot.png"),
        help="Output image path.",
    )
    p.add_argument(
        "--species",
        nargs="+",
        default=None,
        help="Optional list of species IDs to include.",
    )
    p.add_argument(
        "--max-ions",
        type=int,
        default=120,
        help="Hard cap on ions to plot after filtering.",
    )
    p.add_argument(
        "--max-per-species",
        type=int,
        default=40,
        help="Cap ions per species before the global cap is applied.",
    )
    p.add_argument(
        "--time-stride",
        type=int,
        default=5,
        help="Use every Nth timestep to reduce memory load.",
    )
    p.add_argument(
        "--max-frames",
        type=int,
        default=600,
        help="Optional limit on number of frames after striding.",
    )
    p.add_argument(
        "--rng-seed",
        type=int,
        default=0,
        help="Seed for ion subsampling.",
    )
    p.add_argument(
        "--no-3d",
        action="store_true",
        help="Skip the 3D panel (useful when running headless/slow).",
    )
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
            raise RuntimeError("No ions selected; broaden species filter or caps.")

        time, positions = load_positions_subset(
            traj, ion_indices=ion_indices, time_stride=args.time_stride, max_frames=args.max_frames
        )

    ion_species = np.array(species_ids)[ion_indices]
    color_map = species_color_map(ion_species)
    ion_colors = np.array([color_map[sp] for sp in ion_species])

    fig = plt.figure(figsize=(12, 4))
    ax_xy = fig.add_subplot(1, 3 if not args.no_3d else 2, 1)
    ax_xz = fig.add_subplot(1, 3 if not args.no_3d else 2, 2)
    ax_3d = None
    if not args.no_3d:
        ax_3d = fig.add_subplot(1, 3, 3, projection="3d")

    _plot_projection(ax_xy, positions, ion_colors, ion_species, proj=(0, 1), title="XY projection")
    _plot_projection(ax_xz, positions, ion_colors, ion_species, proj=(0, 2), title="XZ projection")
    if ax_3d is not None:
        _plot_3d(ax_3d, positions, ion_colors, ion_species)

    handles, labels = _legend_handles(color_map)
    fig.legend(handles, labels, loc="upper center", ncol=min(4, len(labels)), bbox_to_anchor=(0.5, 1.05))
    fig.suptitle(f"Trajectories from {args.traj.name} (T={positions.shape[0]}, ions={positions.shape[1]})")
    fig.tight_layout()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(args.out, dpi=200)
    print(f"Wrote {args.out}")
    return 0


def _plot_projection(ax, positions, colors, species, proj, title):
    x_idx, y_idx = proj
    for i in range(positions.shape[1]):
        xy = positions[:, i, [x_idx, y_idx]]
        ax.plot(xy[:, 0], xy[:, 1], color=colors[i], alpha=0.8, linewidth=0.9)
    ax.set_xlabel(_axis_label(x_idx))
    ax.set_ylabel(_axis_label(y_idx))
    ax.set_title(title)
    ax.grid(True, alpha=0.2)

    # Center the limits around observed min/max with a small margin.
    xy_all = positions[:, :, [x_idx, y_idx]].reshape(-1, 2)
    mins = xy_all.min(axis=0)
    maxs = xy_all.max(axis=0)
    span = (maxs - mins) * 0.05
    ax.set_xlim(mins[0] - span[0], maxs[0] + span[0])
    ax.set_ylim(mins[1] - span[1], maxs[1] + span[1])


def _plot_3d(ax, positions, colors, species):
    for i in range(positions.shape[1]):
        traj = positions[:, i, :]
        ax.plot(traj[:, 0], traj[:, 1], traj[:, 2], color=colors[i], alpha=0.7, linewidth=0.9)
    ax.set_xlabel("X [m]")
    ax.set_ylabel("Y [m]")
    ax.set_zlabel("Z [m]")
    ax.set_title("3D view")

    mins = positions.reshape(-1, 3).min(axis=0)
    maxs = positions.reshape(-1, 3).max(axis=0)
    span = (maxs - mins) * 0.05
    ax.set_xlim(mins[0] - span[0], maxs[0] + span[0])
    ax.set_ylim(mins[1] - span[1], maxs[1] + span[1])
    ax.set_zlim(mins[2] - span[2], maxs[2] + span[2])


def _legend_handles(color_map):
    import matplotlib.patches as mpatches

    handles = []
    labels = []
    for sp, color in color_map.items():
        handles.append(mpatches.Patch(color=color, label=sp))
        labels.append(sp)
    return handles, labels


def _axis_label(idx: int) -> str:
    return ["X [m]", "Y [m]", "Z [m]"][idx]


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Create a quick animation (GIF/MP4) of trajectory projections from an ICARION HDF5 file.
Uses ion/time subsampling to stay light-weight on large runs.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np

from analysis.common import (
    load_positions_subset,
    load_species_ids,
    open_trajectory,
    select_ion_indices,
    species_color_map,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Animate a projection of ICARION trajectories.")
    p.add_argument(
        "--traj",
        type=Path,
        default=Path("results/lqit/lqit_trajectories.h5"),
        help="Path to trajectory HDF5 file.",
    )
    p.add_argument(
        "--out",
        type=Path,
        default=Path("analysis/output/trajectory.gif"),
        help="Output file (.gif recommended unless ffmpeg is available).",
    )
    p.add_argument(
        "--projection",
        choices=["xy", "xz", "yz"],
        default="xy",
        help="Plane to animate.",
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
        default=80,
        help="Hard cap on ions to animate after filtering.",
    )
    p.add_argument(
        "--max-per-species",
        type=int,
        default=25,
        help="Cap ions per species before the global cap is applied.",
    )
    p.add_argument(
        "--time-stride",
        type=int,
        default=4,
        help="Use every Nth timestep when reading from disk.",
    )
    p.add_argument(
        "--frame-step",
        type=int,
        default=2,
        help="Use every Nth frame after loading (reduces total frames).",
    )
    p.add_argument(
        "--max-frames",
        type=int,
        default=400,
        help="Optional max frames after striding; shorten long simulations.",
    )
    p.add_argument(
        "--rng-seed",
        type=int,
        default=0,
        help="Seed for ion subsampling.",
    )
    p.add_argument(
        "--trail",
        type=int,
        default=4,
        help="Number of past frames to show as a faint trail.",
    )
    p.add_argument(
        "--fps",
        type=int,
        default=15,
        help="Frames per second for the animation.",
    )
    p.add_argument(
        "--writer",
        choices=["pillow", "ffmpeg"],
        default="pillow",
        help="Animation writer; pillow -> GIF, ffmpeg -> MP4 (requires ffmpeg).",
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

    # Apply in-memory frame skipping after the initial stride.
    frame_step = max(1, args.frame_step)
    time = time[::frame_step]
    positions = positions[::frame_step, :, :]

    if positions.shape[0] == 0 or positions.shape[1] == 0:
        raise RuntimeError("No frames or ions available after filtering/striding.")

    proj_axes = _projection_axes(args.projection)
    axis_labels = ["X", "Y", "Z"]
    proj_label = f"{axis_labels[proj_axes[0]]}{axis_labels[proj_axes[1]]}"
    xy = positions[:, :, proj_axes]
    mins = xy.reshape(-1, 2).min(axis=0)
    maxs = xy.reshape(-1, 2).max(axis=0)
    pad = (maxs - mins) * 0.05

    ion_species = np.array(species_ids)[ion_indices]
    color_map = species_color_map(ion_species)
    ion_colors = np.array([color_map[sp] for sp in ion_species])

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_xlabel(f"{axis_labels[proj_axes[0]]} [m]")
    ax.set_ylabel(f"{axis_labels[proj_axes[1]]} [m]")
    ax.set_xlim(mins[0] - pad[0], maxs[0] + pad[0])
    ax.set_ylim(mins[1] - pad[1], maxs[1] + pad[1])
    ax.set_title(f"{proj_label} projection ({positions.shape[1]} ions)")
    ax.grid(True, alpha=0.2)
    ax.set_aspect("equal", adjustable="box")

    trail_scat = ax.scatter([], [], s=8, alpha=0.35, linewidths=0)
    scat = ax.scatter([], [], s=12, linewidths=0, alpha=0.9)

    def init():
        trail_scat.set_offsets(np.empty((0, 2)))
        scat.set_offsets(np.empty((0, 2)))
        return trail_scat, scat

    def update(frame_idx):
        coords = xy[frame_idx]
        scat.set_offsets(coords)
        scat.set_color(ion_colors)
        scat.set_alpha(0.9)

        trail_window = max(1, args.trail)
        start = max(0, frame_idx - trail_window + 1)
        trail_coords = xy[start : frame_idx + 1].reshape(-1, 2)
        trail_colors = np.tile(ion_colors, (frame_idx - start + 1, 1))
        trail_scat.set_offsets(trail_coords)
        trail_scat.set_color(trail_colors)
        trail_scat.set_alpha(0.35)

        ax.set_title(
            f"{proj_label} projection t={time[frame_idx]*1e6:.2f} us "
            f"({frame_idx + 1}/{xy.shape[0]} frames, {positions.shape[1]} ions)"
        )
        return trail_scat, scat

    anim = animation.FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=xy.shape[0],
        interval=1000 / args.fps,
        blit=True,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    writer = _resolve_writer(args.writer, fps=args.fps)
    anim.save(args.out, writer=writer, dpi=150)
    print(f"Wrote animation to {args.out}")
    return 0


def _projection_axes(proj: str) -> tuple[int, int]:
    if proj == "xy":
        return (0, 1)
    if proj == "xz":
        return (0, 2)
    if proj == "yz":
        return (1, 2)
    raise ValueError(f"Unsupported projection: {proj}")


def _resolve_writer(name: str, fps: int):
    if name == "ffmpeg":
        return animation.FFMpegWriter(fps=fps, bitrate=1800)
    return animation.PillowWriter(fps=fps)


if __name__ == "__main__":
    raise SystemExit(main())

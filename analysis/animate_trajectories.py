#!/usr/bin/env python3
"""
Create a quick animation (GIF/MP4) of trajectory projections from an ICARION HDF5 file.
Uses ion/time subsampling to stay light-weight on large runs.
"""

from __future__ import annotations

import argparse
import os
from pathlib import Path

import matplotlib


def _configure_matplotlib_backend() -> None:
    if os.environ.get("MPLBACKEND"):
        return
    parser = argparse.ArgumentParser(add_help=False)
    parser.add_argument("--backend", default=None)
    args, _ = parser.parse_known_args()
    if args.backend:
        matplotlib.use(args.backend)
    else:
        matplotlib.use("Agg")


_configure_matplotlib_backend()
import matplotlib.pyplot as plt
from matplotlib import animation
from matplotlib import colors as mcolors
import numpy as np

# Allow running as a standalone script without installing the package
if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.common import (
    load_positions_subset,
    open_trajectory,
    select_ion_indices,
)


_DARK_THEME = {
    "bg": "#14181C",
    "grid": "white",
    "text": "white",
    "grid_alpha": 0.2,
    "legend_face": "#14181C",
}
_LIGHT_THEME = {
    "bg": "white",
    "grid": "#333333",
    "text": "#111111",
    "grid_alpha": 0.2,
    "legend_face": "white",
}
_DEFAULT_COLOR_SCHEME = {
    "default": "#00D9FF",
    "H3O+": "#FF006E",
    "NH4+": "#FFBE0B",
    "NO+": "#8338EC",
    "O2+": "#3A86FF",
}
_SPECIES_DISPLAY = {
    "H3O+": r"$\mathrm{H_3O^+}$",
    "PentanalH+": r"$\mathrm{PentanalH^+}$",
}


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
        choices=["xy", "xz", "yz", "multi"],
        default="xy",
        help="Plane to animate (or 'multi' for XY/XZ/YZ panels).",
    )
    p.add_argument(
        "--style",
        choices=["modern", "trail", "glow"],
        default="modern",
        help="Visual style (default: modern).",
    )
    p.add_argument(
        "--theme",
        choices=["dark", "light"],
        default="dark",
        help="Background/theme style.",
    )
    p.add_argument(
        "--z-right",
        dest="z_right",
        action="store_true",
        default=True,
        help="Plot XZ/YZ with Z axis on the horizontal axis (default).",
    )
    p.add_argument(
        "--z-up",
        dest="z_right",
        action="store_false",
        help="Plot XZ/YZ with Z axis on the vertical axis (legacy behavior).",
    )
    p.add_argument(
        "--species",
        nargs="+",
        default=None,
        help="Optional list of species IDs to include.",
    )
    p.add_argument(
        "--species-color",
        action="append",
        default=None,
        help="Override species color as SPECIES=COLOR (e.g., H3O+=#ff006e). Repeatable.",
    )
    p.add_argument(
        "--species-frame",
        choices=["last", "first"],
        default="last",
        help="Frame used for species-based selection (default: last; better for reactions).",
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
        "--time-range-us",
        default=None,
        help="Time range to animate in microseconds, e.g. '0:50' or '10:'.",
    )
    p.add_argument(
        "--time-range",
        dest="time_range_us",
        default=None,
        help="Alias for --time-range-us.",
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
        help="Number of past frames to show as a faint trail (trail style only).",
    )
    p.add_argument(
        "--fps",
        type=int,
        default=15,
        help="Frames per second for the animation.",
    )
    p.add_argument(
        "--writer",
        choices=["auto", "pillow", "ffmpeg"],
        default="auto",
        help="Animation writer; auto picks ffmpeg for .mp4 when available, else pillow.",
    )
    p.add_argument(
        "--backend",
        default=None,
        help="Matplotlib backend (default: Agg unless MPLBACKEND is set).",
    )
    p.add_argument(
        "--dpi",
        type=int,
        default=150,
        help="Output DPI (lower = faster).",
    )
    p.add_argument(
        "--figscale",
        type=float,
        default=1.0,
        help="Scale factor for figure size (e.g., 1.5 or 2.0 for sharper GIFs).",
    )
    p.add_argument("--label-fontsize", type=float, default=11.0, help="Axis label font size.")
    p.add_argument("--tick-fontsize", type=float, default=10.0, help="Tick label font size.")
    p.add_argument("--title-fontsize", type=float, default=12.0, help="Animation title font size.")
    p.add_argument("--legend-fontsize", type=float, default=10.0, help="Legend font size.")
    p.add_argument("--spine-width", type=float, default=1.1, help="Axis spine line width.")
    p.add_argument("--grid-width", type=float, default=0.7, help="Grid line width.")
    p.add_argument("--marker-size", type=float, default=12.0, help="Core marker size.")
    p.add_argument("--trail-size", type=float, default=7.0, help="Trail marker size.")
    p.add_argument("--glow-size", type=float, default=60.0, help="Glow marker size.")
    p.add_argument(
        "--display-scale-x",
        "--x-scale",
        dest="display_scale_x",
        type=float,
        default=1.0,
        help="Multiply X coordinates for display only (default: 1.0).",
    )
    p.add_argument(
        "--display-scale-y",
        dest="display_scale_y",
        type=float,
        default=1.0,
        help="Multiply Y coordinates for display only (default: 1.0).",
    )
    p.add_argument(
        "--display-scale-z",
        dest="display_scale_z",
        type=float,
        default=1.0,
        help="Multiply Z coordinates for display only (default: 1.0).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.traj.exists():
        raise FileNotFoundError(f"Trajectory file not found: {args.traj}")

    with open_trajectory(args.traj) as h5:
        traj = h5["trajectory"]
        species_for_selection = _load_species_for_selection(traj, args.species_frame)
        ion_indices = select_ion_indices(
            species_for_selection,
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
        species_series = _load_species_series(
            traj,
            ion_indices=ion_indices,
            time_stride=args.time_stride,
            max_frames=args.max_frames,
        )

    if args.time_range_us:
        time, positions, species_series = _apply_time_range(time, positions, args.time_range_us, species_series)

    # Apply in-memory frame skipping after the initial stride.
    frame_step = max(1, args.frame_step)
    time = time[::frame_step]
    positions = positions[::frame_step, :, :]
    species_series = species_series[::frame_step, :]

    if positions.shape[0] == 0 or positions.shape[1] == 0:
        raise RuntimeError("No frames or ions available after filtering/striding.")
    if species_series.shape[:2] != positions.shape[:2]:
        raise RuntimeError(
            "Species series shape does not match positions after slicing: "
            f"{species_series.shape[:2]} vs {positions.shape[:2]}"
        )

    color_overrides = _parse_species_color_overrides(args.species_color)
    color_map = _build_color_map(species_series, color_overrides)
    colors_over_time = _colors_over_time(species_series, color_map)
    initial_colors = colors_over_time[0]

    theme = _theme_config(args.theme)
    use_trail = args.style == "trail" and args.trail > 0
    use_glow = args.style == "glow"
    trail_window = max(0, args.trail) if use_trail else 0
    num_ions = positions.shape[1]

    projections = ["xy", "xz", "yz"] if args.projection == "multi" else [args.projection]
    proj_views = []
    axis_meta = []
    for proj in projections:
        x_idx, y_idx, x_label, y_label = _projection_axes(proj, args.z_right)
        scales = np.array(
            [
                _axis_display_scale(x_idx, args),
                _axis_display_scale(y_idx, args),
            ],
            dtype=float,
        )
        view = positions[:, :, [x_idx, y_idx]].astype(float)
        view *= scales[None, None, :]
        x_label = _scaled_axis_label(x_label, scales[0])
        y_label = _scaled_axis_label(y_label, scales[1])
        mins = view.reshape(-1, 2).min(axis=0)
        maxs = view.reshape(-1, 2).max(axis=0)
        pad = (maxs - mins) * 0.05
        proj_views.append(view)
        axis_meta.append((x_label, y_label, mins, maxs, pad))

    figscale = max(0.25, float(args.figscale))
    if args.projection == "multi":
        fig, axes = plt.subplots(1, 3, figsize=(14 * figscale, 5 * figscale))
        fig.subplots_adjust(top=0.85, wspace=0.25)
        axes = list(axes)
    else:
        fig, ax = plt.subplots(figsize=(6.5 * figscale, 6 * figscale))
        axes = [ax]

    _apply_theme(fig, axes, theme, args)

    scatters = []
    trail_scatters = []
    glow_scatters = []
    for ax, meta in zip(axes, axis_meta):
        x_label, y_label, mins, maxs, pad = meta
        proj_label = _projection_label(x_label, y_label)
        ax.set_xlabel(x_label, fontsize=args.label_fontsize)
        ax.set_ylabel(y_label, fontsize=args.label_fontsize)
        ax.set_xlim(mins[0] - pad[0], maxs[0] + pad[0])
        ax.set_ylim(mins[1] - pad[1], maxs[1] + pad[1])
        ax.set_aspect("equal", adjustable="box")

        if use_glow:
            glow = ax.scatter([], [], s=args.glow_size, alpha=0.18, linewidths=0)
            glow.set_color(initial_colors)
            glow_scatters.append(glow)
        else:
            glow_scatters.append(None)

        if use_trail:
            trail = ax.scatter([], [], s=args.trail_size, alpha=0.35, linewidths=0)
            trail_scatters.append(trail)
        else:
            trail_scatters.append(None)

        core = ax.scatter([], [], s=args.marker_size, linewidths=0, alpha=0.9)
        core.set_color(initial_colors)
        scatters.append(core)

    if args.projection == "multi":
        title_text = fig.suptitle("", color=theme["text"], fontsize=args.title_fontsize, weight="bold")
    else:
        title_text = axes[0].text(
            0.5,
            1.04,
            "",
            transform=axes[0].transAxes,
            ha="center",
            color=theme["text"],
            fontsize=args.title_fontsize,
            weight="bold",
        )

    _maybe_add_legend(axes[0], color_map, theme, args)

    def init():
        artists = []
        for idx in range(len(axes)):
            if use_trail:
                trail_scatters[idx].set_offsets(np.empty((0, 2)))
                artists.append(trail_scatters[idx])
            if use_glow:
                glow_scatters[idx].set_offsets(np.empty((0, 2)))
                artists.append(glow_scatters[idx])
            scatters[idx].set_offsets(np.empty((0, 2)))
            artists.append(scatters[idx])
        artists.append(title_text)
        return artists

    def update(frame_idx):
        artists = []
        colors_frame = colors_over_time[frame_idx]
        for idx, view in enumerate(proj_views):
            coords = view[frame_idx]
            scatters[idx].set_offsets(coords)
            scatters[idx].set_color(colors_frame)
            artists.append(scatters[idx])

            if use_glow:
                glow_scatters[idx].set_offsets(coords)
                glow_scatters[idx].set_color(colors_frame)
                artists.append(glow_scatters[idx])

            if use_trail:
                start = max(0, frame_idx - trail_window + 1)
                trail_len = frame_idx - start + 1
                trail_coords = view[start : frame_idx + 1].reshape(-1, 2)
                trail_scatters[idx].set_offsets(trail_coords)
                trail_colors = colors_over_time[start : frame_idx + 1].reshape(-1, 4)
                trail_scatters[idx].set_color(trail_colors)
                artists.append(trail_scatters[idx])

        title_text.set_text(
            f"{_projection_title(args.projection, args.z_right)} "
            f"t={time[frame_idx]*1e6:.2f} us "
            f"({frame_idx + 1}/{proj_views[0].shape[0]} frames, {positions.shape[1]} ions)"
        )
        artists.append(title_text)
        return artists

    anim = animation.FuncAnimation(
        fig,
        update,
        init_func=init,
        frames=proj_views[0].shape[0],
        interval=1000 / args.fps,
        blit=True,
    )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    writer = _resolve_writer(args.writer, fps=args.fps, out_path=args.out)
    anim.save(args.out, writer=writer, dpi=args.dpi)
    print(f"Wrote animation to {args.out}")
    return 0


def _projection_label(x_label: str, y_label: str) -> str:
    return f"{x_label.split()[0]}{y_label.split()[0]}"


def _apply_time_range(
    time: np.ndarray,
    positions: np.ndarray,
    time_range_us: str,
    species_series: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray | None]:
    parts = time_range_us.split(":")
    if len(parts) != 2:
        raise ValueError("time range must look like 'start:end' in microseconds")
    start_str, end_str = (part.strip() for part in parts)
    t_us = time * 1e6
    start_idx = 0
    end_idx = len(t_us)
    if start_str:
        start_us = float(start_str)
        start_idx = int(np.searchsorted(t_us, start_us, side="left"))
    if end_str:
        end_us = float(end_str)
        end_idx = int(np.searchsorted(t_us, end_us, side="right"))
    time_sel = time[start_idx:end_idx]
    pos_sel = positions[start_idx:end_idx, :, :]
    if species_series is None:
        return time_sel, pos_sel, None
    sp_sel = species_series[start_idx:end_idx, :]
    return time_sel, pos_sel, sp_sel


def _projection_title(proj: str, z_right: bool) -> str:
    if proj == "multi":
        return "Multi-view"
    _, _, x_label, y_label = _projection_axes(proj, z_right)
    return f"{_projection_label(x_label, y_label)} projection"


def _projection_axes(proj: str, z_right: bool) -> tuple[int, int, str, str]:
    if proj == "xy":
        return (0, 1, "X [m]", "Y [m]")
    if proj == "xz":
        if z_right:
            return (2, 0, "Z [m]", "X [m]")
        return (0, 2, "X [m]", "Z [m]")
    if proj == "yz":
        if z_right:
            return (2, 1, "Z [m]", "Y [m]")
        return (1, 2, "Y [m]", "Z [m]")
    raise ValueError(f"Unsupported projection: {proj}")


def _theme_config(name: str) -> dict[str, str | float]:
    if name == "light":
        return _LIGHT_THEME
    return _DARK_THEME


def _apply_theme(fig, axes, theme, args: argparse.Namespace):
    fig.patch.set_facecolor(theme["bg"])
    for ax in axes:
        ax.set_facecolor(theme["bg"])
        ax.tick_params(colors=theme["text"], labelsize=args.tick_fontsize)
        for spine in ax.spines.values():
            spine.set_color(theme["text"])
            spine.set_linewidth(args.spine_width)
        ax.grid(True, alpha=theme["grid_alpha"], color=theme["grid"], linewidth=args.grid_width)


def _build_color_map(
    species: np.ndarray,
    overrides: dict[str, tuple[float, float, float, float]] | None = None,
) -> dict[str, tuple[float, float, float, float]]:
    flat = species.reshape(-1)
    unique: list[str] = []
    for sp in flat:
        sp_str = _decode_species(sp)
        if sp_str not in unique:
            unique.append(sp_str)

    color_map: dict[str, tuple[float, float, float, float]] = {}
    overrides = overrides or {}
    unknown = []
    for sp in unique:
        if sp in overrides:
            color_map[sp] = overrides[sp]
            continue
        if sp in _DEFAULT_COLOR_SCHEME:
            color_map[sp] = mcolors.to_rgba(_DEFAULT_COLOR_SCHEME[sp])
        else:
            unknown.append(sp)

    if unknown:
        cmap = plt.cm.get_cmap("tab20", len(unknown))
        for idx, sp in enumerate(unknown):
            color_map[sp] = cmap(idx)

    return color_map


def _parse_species_color_overrides(entries: list[str] | None) -> dict[str, tuple[float, float, float, float]]:
    if not entries:
        return {}
    overrides: dict[str, tuple[float, float, float, float]] = {}
    for entry in entries:
        if "=" not in entry:
            raise ValueError(f"--species-color must look like SPECIES=COLOR (got {entry!r})")
        species, color = entry.split("=", 1)
        species = species.strip()
        color = color.strip()
        if not species or not color:
            raise ValueError(f"--species-color must look like SPECIES=COLOR (got {entry!r})")
        overrides[species] = mcolors.to_rgba(color)
    return overrides


def _decode_species(val) -> str:
    if isinstance(val, (bytes, bytearray)):
        return val.decode("utf-8")
    return str(val)


def _decode_species_array(values: np.ndarray) -> np.ndarray:
    return np.asarray([_decode_species(v) for v in values], dtype=object)


def _species_row_for_selection(species_ds, which: str) -> np.ndarray:
    if species_ds.ndim == 1:
        return _decode_species_array(np.asarray(species_ds[()]))
    if which == "first":
        return _decode_species_array(np.asarray(species_ds[0]))
    return _decode_species_array(np.asarray(species_ds[-1]))


def _load_species_for_selection(traj, which: str) -> np.ndarray:
    if "species_ids" not in traj:
        raise KeyError("trajectory/species_ids is required for species-aware selection.")
    species_ds = traj["species_ids"]
    return _species_row_for_selection(species_ds, which)


def _time_indices(n_time: int, stride: int, max_frames: int | None) -> np.ndarray:
    stride = max(1, int(stride))
    indices = np.arange(n_time, dtype=int)[::stride]
    if max_frames is not None and indices.size > int(max_frames):
        indices = indices[: int(max_frames)]
    return indices


def _load_species_series(
    traj,
    ion_indices: np.ndarray,
    time_stride: int,
    max_frames: int | None,
) -> np.ndarray:
    if "species_ids" not in traj:
        # Fallback: unknown species -> empty strings.
        n_time = int(traj["time"].shape[0])
        time_idx = _time_indices(n_time, stride=time_stride, max_frames=max_frames)
        return np.full((time_idx.size, int(len(ion_indices))), "", dtype=object)

    species_ds = traj["species_ids"]
    n_time = int(traj["time"].shape[0])
    time_idx = _time_indices(n_time, stride=time_stride, max_frames=max_frames)

    if species_ds.ndim == 1:
        row = _decode_species_array(np.asarray(species_ds[()]))[ion_indices]
        return np.tile(row[None, :], (time_idx.size, 1))

    species_subset = np.asarray(species_ds[time_idx][:, ion_indices])
    return _decode_species_array(species_subset.reshape(-1)).reshape(species_subset.shape)


def _colors_over_time(species_series: np.ndarray, color_map: dict[str, tuple[float, float, float, float]]) -> np.ndarray:
    colors = np.empty(species_series.shape + (4,), dtype=float)
    for sp, rgba in color_map.items():
        mask = species_series == sp
        colors[mask] = rgba
    return colors


def _axis_display_scale(axis_idx: int, args: argparse.Namespace) -> float:
    if axis_idx == 0:
        return float(args.display_scale_x)
    if axis_idx == 1:
        return float(args.display_scale_y)
    if axis_idx == 2:
        return float(args.display_scale_z)
    return 1.0


def _scaled_axis_label(label: str, scale: float) -> str:
    if np.isclose(scale, 1.0):
        return label
    return f"{label} (x{scale:g} display)"


def _maybe_add_legend(ax, color_map, theme, args: argparse.Namespace):
    import matplotlib.patches as mpatches

    if len(color_map) > 8:
        return
    handles = []
    labels = []
    for sp, color in color_map.items():
        disp = _display_species(sp)
        handles.append(mpatches.Patch(color=color, label=disp))
        labels.append(disp)
    legend = ax.legend(
        handles=handles,
        labels=labels,
        loc="upper right",
        framealpha=0.85,
        fontsize=args.legend_fontsize,
    )
    legend.get_frame().set_facecolor(theme["legend_face"])
    legend.get_frame().set_edgecolor(theme["text"])
    for text in legend.get_texts():
        text.set_color(theme["text"])


def _display_species(species: str) -> str:
    return _SPECIES_DISPLAY.get(species, species)


def _resolve_writer(name: str, fps: int, out_path: Path):
    suffix = out_path.suffix.lower()
    ffmpeg_available = animation.FFMpegWriter.isAvailable()

    if name == "auto":
        if suffix in {".mp4", ".m4v", ".mov"} and ffmpeg_available:
            return animation.FFMpegWriter(fps=fps, bitrate=1800)
        return animation.PillowWriter(fps=fps)

    if name == "ffmpeg":
        if not ffmpeg_available:
            raise RuntimeError("ffmpeg writer requested but ffmpeg is not available.")
        return animation.FFMpegWriter(fps=fps, bitrate=1800)

    return animation.PillowWriter(fps=fps)


if __name__ == "__main__":
    raise SystemExit(main())

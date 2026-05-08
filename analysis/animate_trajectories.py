#!/usr/bin/env python3
"""
Create a quick animation (GIF/MP4) of trajectory projections from an ICARION HDF5 file.
Uses ion/time subsampling to stay light-weight on large runs.
"""

from __future__ import annotations

import argparse
import json
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
    "H5O2+": r"$\mathrm{H_3O^+(H_2O)}$",
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
        choices=["dark", "light", "style"],
        default="dark",
        help="Background/theme style. Use 'style' to keep mplstyle colors.",
    )
    p.add_argument(
        "--mplstyle",
        type=Path,
        default=None,
        help="Optional matplotlib style file (.mplstyle) to load before plotting.",
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
        "--legend",
        choices=["auto", "on", "off"],
        default="auto",
        help="Legend mode (auto=show only for <=8 species).",
    )
    p.add_argument(
        "--legend-loc",
        default="upper left",
        help="Legend location (matplotlib loc string).",
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
        "--axis-fit",
        choices=["data", "domain"],
        default="data",
        help="Axis extents from observed data or full domain geometry.",
    )
    p.add_argument(
        "--axis-padding-frac",
        type=float,
        default=0.05,
        help="Relative axis padding for data fit (ignored for domain fit).",
    )
    p.add_argument(
        "--domain-index",
        type=int,
        default=0,
        help="Domain index used for geometry extents in HDF5.",
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
    p.add_argument(
        "--title",
        choices=["on", "off"],
        default="on",
        help="Show or hide dynamic title/time overlay.",
    )
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
    p.add_argument(
        "--twims-config",
        dest="twims_config",
        type=Path,
        default=None,
        help="Path to runtime IMS JSON with TWIMS field_array_terms; enables electrode zone shading.",
    )
    p.add_argument(
        "--tims-config",
        dest="tims_config",
        type=Path,
        default=None,
        help=(
            "Path to runtime IMS JSON with TIMS axial field program; enables axial gradient background shading. "
            "If omitted, a neighboring '*.config.json' next to --traj is used when available."
        ),
    )
    p.add_argument(
        "--twims-ring-z0",
        dest="twims_ring_z0",
        type=float,
        default=0.010,
        help="Z position [m] of the first ring electrode centre (default: 0.010).",
    )
    p.add_argument(
        "--twims-ring-pitch",
        dest="twims_ring_pitch",
        type=float,
        default=0.0013,
        help="Electrode pitch [m] centre-to-centre (default: 0.0013).",
    )
    p.add_argument(
        "--twims-ring-count",
        dest="twims_ring_count",
        type=int,
        default=40,
        help="Total number of ring electrodes (default: 40).",
    )
    p.add_argument(
        "--twims-shading-alpha",
        dest="twims_shading_alpha",
        type=float,
        default=0.18,
        help="Alpha for active TWIMS zone background shading (default: 0.18).",
    )
    p.add_argument(
        "--twims-shading-color",
        dest="twims_shading_color",
        default="#ffcc00",
        help="Color for active TWIMS zone shading (default: #ffcc00).",
    )
    p.add_argument(
        "--tims-shading-alpha-max",
        dest="tims_shading_alpha_max",
        type=float,
        default=0.4,
        help="Maximum alpha for TIMS axial field-gradient shading (default: 0.4).",
    )
    p.add_argument(
        "--tims-shading-color",
        dest="tims_shading_color",
        default="#2a9d8f",
        help="Color for TIMS axial field-gradient shading (default: #2a9d8f).",
    )
    p.add_argument(
        "--tims-shading-slices",
        dest="tims_shading_slices",
        type=int,
        default=96,
        help="Number of background slices used for TIMS alpha gradient shading (default: 96).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.traj.exists():
        raise FileNotFoundError(f"Trajectory file not found: {args.traj}")
    if args.mplstyle is not None:
        if not args.mplstyle.exists():
            raise FileNotFoundError(f"Style file not found: {args.mplstyle}")
        plt.style.use(str(args.mplstyle))

    with open_trajectory(args.traj) as h5:
        traj = h5["trajectory"]
        domain_geometry = _read_domain_geometry(h5, domain_index=args.domain_index)
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

    color_overrides, positional_colors = _parse_species_color_overrides(args.species_color)
    color_map = _build_color_map(species_series, color_overrides, positional_colors)
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

        if domain_geometry is not None:
            # Preserve legacy behavior for "data": always include full Z domain.
            if args.axis_fit == "data":
                mins, maxs = _merge_axis_limits_with_domain(
                    mins=mins,
                    maxs=maxs,
                    x_idx=x_idx,
                    y_idx=y_idx,
                    scales=scales,
                    geometry=domain_geometry,
                    include_radial=False,
                    force_domain_limits=False,
                )
            else:
                mins, maxs = _merge_axis_limits_with_domain(
                    mins=mins,
                    maxs=maxs,
                    x_idx=x_idx,
                    y_idx=y_idx,
                    scales=scales,
                    geometry=domain_geometry,
                    include_radial=True,
                    force_domain_limits=True,
                )

        pad = _axis_padding(mins, maxs, 0.0 if args.axis_fit == "domain" else args.axis_padding_frac)
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
    title_color = _theme_text_color(theme)

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

    title_text = None
    if args.title == "on":
        if args.projection == "multi":
            title_text = fig.suptitle("", color=title_color, fontsize=args.title_fontsize, weight="bold")
        else:
            title_text = axes[0].text(
                0.5,
                1.04,
                "",
                transform=axes[0].transAxes,
                ha="center",
                color=title_color,
                fontsize=args.title_fontsize,
                weight="bold",
            )

    _maybe_add_legend(axes[0], color_map, theme, args)

    adjacent_config = _find_adjacent_runtime_config(args.traj)

    # --- TWIMS high-potential zone background shading ---
    twims_params: list[dict] = []
    twims_patches_per_axis: list[list] = [[] for _ in range(len(axes))]
    ring_to_group: list[int] = []
    twims_config_path = args.twims_config if args.twims_config is not None else adjacent_config
    if twims_config_path is not None:
        twims_params = _load_twims_shading_params(twims_config_path)
        if twims_params:
            n_pg = len(twims_params)
            ring_z = np.array(
                [args.twims_ring_z0 + i * args.twims_ring_pitch for i in range(args.twims_ring_count)]
            )
            ring_to_group = [i % n_pg for i in range(args.twims_ring_count)]
            z_scale = _axis_display_scale(2, args)
            for ax_panel_idx, (ax_p, proj_p) in enumerate(zip(axes, projections)):
                x_idx_proj, y_idx_proj, _, _ = _projection_axes(proj_p, args.z_right)
                if x_idx_proj == 2:
                    z_is_horiz = True
                elif y_idx_proj == 2:
                    z_is_horiz = False
                else:
                    continue  # xy projection has no Z axis
                twims_patches_per_axis[ax_panel_idx] = _make_twims_ring_patches(
                    ax_p, ring_z, args.twims_ring_pitch, z_scale, z_is_horiz,
                    args.twims_shading_color, args.twims_shading_alpha,
                )

    # --- TIMS axial field background shading ---
    tims_shading: dict | None = None
    tims_patches_per_axis: list[list] = [[] for _ in range(len(axes))]
    tims_config_path = args.tims_config if args.tims_config is not None else adjacent_config
    if tims_config_path is not None:
        tims_shading = _load_tims_shading_params(tims_config_path)
        if tims_shading is not None:
            z_scale = _axis_display_scale(2, args)
            slice_edges_global = np.linspace(
                float(tims_shading["z_min_global_m"]),
                float(tims_shading["z_max_global_m"]),
                max(2, int(args.tims_shading_slices)) + 1,
            )
            tims_shading["slice_edges_global_m"] = slice_edges_global
            tims_shading["slice_centers_local_m"] = (
                0.5 * (slice_edges_global[:-1] + slice_edges_global[1:]) - float(tims_shading["origin_z_m"])
            )
            for ax_panel_idx, (ax_p, proj_p) in enumerate(zip(axes, projections)):
                x_idx_proj, y_idx_proj, _, _ = _projection_axes(proj_p, args.z_right)
                if x_idx_proj == 2:
                    z_is_horiz = True
                elif y_idx_proj == 2:
                    z_is_horiz = False
                else:
                    continue
                tims_patches_per_axis[ax_panel_idx] = _make_tims_gradient_patches(
                    ax=ax_p,
                    z_edges_global_m=slice_edges_global,
                    z_scale=z_scale,
                    z_is_horiz=z_is_horiz,
                    color=args.tims_shading_color,
                )

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
        if title_text is not None:
            artists.append(title_text)
        for patches in tims_patches_per_axis:
            for patch in patches:
                artists.append(patch)
        for patches in twims_patches_per_axis:
            for patch in patches:
                artists.append(patch)
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

        if tims_shading is not None:
            alpha_profile = _tims_shading_alpha_profile(
                t_s=float(time[frame_idx]),
                tims=tims_shading,
                alpha_max=float(args.tims_shading_alpha_max),
            )
            for patches in tims_patches_per_axis:
                for patch, alpha in zip(patches, alpha_profile):
                    patch.set_alpha(float(alpha))
                    artists.append(patch)

        if twims_params:
            t = float(time[frame_idx])
            high_state = [
                _twims_is_high(t, tp["phase_rad"], tp["freq_hz"], tp["duty_cycle"])
                for tp in twims_params
            ]
            for patches in twims_patches_per_axis:
                for ring_idx, patch in enumerate(patches):
                    patch.set_visible(high_state[ring_to_group[ring_idx]])
                    artists.append(patch)
        if title_text is not None:
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


def _load_twims_shading_params(config_path: Path) -> list[dict]:
    """Parse TWIMS field_array_terms from a runtime IMS JSON and return a list of
    dicts with keys: phase_rad, freq_hz, duty_cycle."""
    with open(config_path) as f:
        cfg = json.load(f)
    terms = []
    for domain in cfg.get("domains", []):
        for term in domain.get("fields", {}).get("field_array_terms", []):
            if term.get("scale_type") == "TWIMS":
                terms.append({
                    "phase_rad": float(term["twims_phase_rad"]),
                    "freq_hz": float(term["twims_frequency_Hz"]),
                    "duty_cycle": float(term.get("twims_duty_cycle", 0.25)),
                })
    return terms


def _find_adjacent_runtime_config(traj_path: Path) -> Path | None:
    candidates = [
        traj_path.with_suffix(".config.json"),
        traj_path.parent / "trajectories.config.json",
        traj_path.parent / f"{traj_path.stem}.config.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


def _load_tims_shading_params(config_path: Path) -> dict | None:
    with open(config_path) as f:
        cfg = json.load(f)

    for domain in cfg.get("domains", []):
        fields = domain.get("fields", {})
        tims = fields.get("TIMS") or fields.get("tims")
        if not isinstance(tims, dict) or not bool(tims.get("enabled", True)):
            continue

        geometry = domain.get("geometry", {})
        origin = geometry.get("origin_m", [0.0, 0.0, 0.0])
        if len(origin) < 3:
            origin = [0.0, 0.0, 0.0]
        origin_z_m = float(origin[2])
        length_m = float(geometry.get("length_m", 0.0))

        z_positions_local = np.asarray(tims.get("z_positions_m", []), dtype=float)
        if z_positions_local.size < 2:
            return None

        ez_accum_profile = _profile_from_config(tims.get("Ez_accum_profile_V_m"), z_positions_local.size)
        ez_delta_profile = _profile_from_config(tims.get("Ez_delta_profile_V_m"), z_positions_local.size)

        return {
            "origin_z_m": origin_z_m,
            "z_min_global_m": origin_z_m,
            "z_max_global_m": origin_z_m + length_m if length_m > 0.0 else origin_z_m + float(z_positions_local[-1]),
            "z_positions_local_m": z_positions_local,
            "ez_accum_profile_V_m": ez_accum_profile,
            "ez_delta_profile_V_m": ez_delta_profile,
            "ez_accum_uniform_V_m": float(tims.get("Ez_accum_uniform_V_m", 0.0)),
            "ez_delta_uniform_V_m": float(tims.get("Ez_delta_uniform_V_m", 0.0)),
            "ramp_start_s": float(tims.get("ramp_start_s", 0.0)),
            "ramp_end_s": float(tims.get("ramp_end_s", 0.0)),
            "ramp_mode": str(tims.get("ramp_mode", "linear")),
            "ramp_tau_s": float(tims.get("ramp_tau_s", 0.0)),
            "ramp_fraction": tims.get("ramp_fraction"),
        }

    return None


def _profile_from_config(values, expected_len: int) -> np.ndarray | None:
    if values is None:
        return None
    arr = np.asarray(values, dtype=float)
    if arr.size != expected_len:
        return None
    return arr


def _tims_ramp_fraction(t_s: float, tims: dict) -> float:
    raw_fraction = tims.get("ramp_fraction")
    if isinstance(raw_fraction, (int, float)):
        return float(np.clip(raw_fraction, 0.0, 1.0))

    t0 = float(tims["ramp_start_s"])
    t1 = float(tims["ramp_end_s"])
    if t_s <= t0:
        return 0.0
    if t_s >= t1:
        return 1.0

    duration = t1 - t0
    if duration <= 0.0:
        return 1.0

    rel = t_s - t0
    if str(tims.get("ramp_mode", "linear")) == "exponential":
        tau = max(float(tims.get("ramp_tau_s", 0.0)), 1e-30)
        denom = 1.0 - np.exp(-duration / tau)
        if abs(denom) < 1e-15:
            return float(np.clip(rel / duration, 0.0, 1.0))
        frac = (1.0 - np.exp(-rel / tau)) / denom
        return float(np.clip(frac, 0.0, 1.0))

    return float(np.clip(rel / duration, 0.0, 1.0))


def _sample_tims_profile(
    z_local_m: np.ndarray,
    z_nodes_m: np.ndarray,
    values_V_m: np.ndarray | None,
    uniform_fallback_V_m: float,
) -> np.ndarray:
    if values_V_m is None or z_nodes_m.size == 0 or values_V_m.size != z_nodes_m.size:
        return np.full_like(z_local_m, float(uniform_fallback_V_m), dtype=float)
    return np.interp(z_local_m, z_nodes_m, values_V_m, left=values_V_m[0], right=values_V_m[-1])


def _tims_shading_alpha_profile(t_s: float, tims: dict, alpha_max: float) -> np.ndarray:
    z_local = np.asarray(tims["slice_centers_local_m"], dtype=float)
    z_nodes = np.asarray(tims["z_positions_local_m"], dtype=float)
    ramp_fraction = _tims_ramp_fraction(t_s, tims)
    e_accum = _sample_tims_profile(
        z_local_m=z_local,
        z_nodes_m=z_nodes,
        values_V_m=tims.get("ez_accum_profile_V_m"),
        uniform_fallback_V_m=float(tims.get("ez_accum_uniform_V_m", 0.0)),
    )
    e_delta = _sample_tims_profile(
        z_local_m=z_local,
        z_nodes_m=z_nodes,
        values_V_m=tims.get("ez_delta_profile_V_m"),
        uniform_fallback_V_m=float(tims.get("ez_delta_uniform_V_m", 0.0)),
    )

    # Keep a fixed normalization so ramp-driven amplitude changes remain visible.
    if "_alpha_norm_max" not in tims:
        e_t0 = e_accum
        e_t1 = e_accum - e_delta
        tims["_alpha_norm_max"] = float(max(np.max(np.abs(e_t0)), np.max(np.abs(e_t1)), 1e-30))

    magnitude = np.abs(e_accum - ramp_fraction * e_delta)
    norm = float(tims.get("_alpha_norm_max", 0.0))
    if norm <= 0.0:
        return np.zeros_like(magnitude)
    return np.clip(magnitude / norm, 0.0, 1.0) * max(0.0, float(alpha_max))


def _twims_is_high(t_s: float, phase_rad: float, freq_hz: float, duty_cycle: float) -> bool:
    """Return True if the asymmetric-square TWIMS waveform is in the high state at time t_s."""
    wrapped = (2.0 * np.pi * freq_hz * t_s + phase_rad) % (2.0 * np.pi)
    return bool(wrapped < duty_cycle * 2.0 * np.pi)


def _make_twims_ring_patches(
    ax,
    ring_z_positions: np.ndarray,
    pitch: float,
    z_scale: float,
    z_is_horiz: bool,
    color: str,
    alpha: float,
) -> list:
    """Create one invisible Rectangle patch per ring electrode on *ax*.
    Patches span the full perpendicular axis via a blended transform so they
    stay correct after axis-limit changes."""
    from matplotlib.patches import Rectangle
    from matplotlib.transforms import blended_transform_factory

    patches = []
    for z in ring_z_positions:
        z_lo = (float(z) - pitch * 0.5) * z_scale
        z_hi = (float(z) + pitch * 0.5) * z_scale
        if z_is_horiz:
            transform = blended_transform_factory(ax.transData, ax.transAxes)
            rect = Rectangle(
                (z_lo, 0.0), width=z_hi - z_lo, height=1.0,
                transform=transform,
                facecolor=color, alpha=alpha, linewidth=0,
                visible=False, zorder=0.5,
            )
        else:
            transform = blended_transform_factory(ax.transAxes, ax.transData)
            rect = Rectangle(
                (0.0, z_lo), width=1.0, height=z_hi - z_lo,
                transform=transform,
                facecolor=color, alpha=alpha, linewidth=0,
                visible=False, zorder=0.5,
            )
        ax.add_patch(rect)
        patches.append(rect)
    return patches


def _make_tims_gradient_patches(
    ax,
    z_edges_global_m: np.ndarray,
    z_scale: float,
    z_is_horiz: bool,
    color: str,
) -> list:
    from matplotlib.patches import Rectangle
    from matplotlib.transforms import blended_transform_factory

    patches = []
    for z_lo_global, z_hi_global in zip(z_edges_global_m[:-1], z_edges_global_m[1:]):
        z_lo = float(z_lo_global) * z_scale
        z_hi = float(z_hi_global) * z_scale
        if z_is_horiz:
            transform = blended_transform_factory(ax.transData, ax.transAxes)
            rect = Rectangle(
                (z_lo, 0.0), width=z_hi - z_lo, height=1.0,
                transform=transform,
                facecolor=color, alpha=0.0, linewidth=0,
                visible=True, zorder=0.1,
            )
        else:
            transform = blended_transform_factory(ax.transAxes, ax.transData)
            rect = Rectangle(
                (0.0, z_lo), width=1.0, height=z_hi - z_lo,
                transform=transform,
                facecolor=color, alpha=0.0, linewidth=0,
                visible=True, zorder=0.1,
            )
        ax.add_patch(rect)
        patches.append(rect)
    return patches


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
    if name == "style":
        return {}
    if name == "light":
        return _LIGHT_THEME
    return _DARK_THEME


def _apply_theme(fig, axes, theme, args: argparse.Namespace):
    if not theme:
        for ax in axes:
            ax.tick_params(labelsize=args.tick_fontsize)
            for spine in ax.spines.values():
                spine.set_linewidth(args.spine_width)
        return
    fig.patch.set_facecolor(theme["bg"])
    for ax in axes:
        ax.set_facecolor(theme["bg"])
        ax.tick_params(colors=theme["text"], labelsize=args.tick_fontsize)
        for spine in ax.spines.values():
            spine.set_color(theme["text"])
            spine.set_linewidth(args.spine_width)
        ax.grid(True, alpha=theme["grid_alpha"], color=theme["grid"], linewidth=args.grid_width)


def _theme_text_color(theme: dict[str, str | float]) -> str:
    if theme and "text" in theme:
        return str(theme["text"])
    return str(plt.rcParams.get("text.color", "#111111"))


def _build_color_map(
    species: np.ndarray,
    overrides: dict[str, tuple[float, float, float, float]] | None = None,
    positional_colors: list[tuple[float, float, float, float]] | None = None,
) -> dict[str, tuple[float, float, float, float]]:
    flat = species.reshape(-1)
    unique: list[str] = []
    for sp in flat:
        sp_str = _decode_species(sp)
        if sp_str not in unique:
            unique.append(sp_str)

    color_map: dict[str, tuple[float, float, float, float]] = {}
    overrides = overrides or {}
    positional_colors = positional_colors or []
    unknown = []
    positional_idx = 0
    for sp in unique:
        if sp in overrides:
            color_map[sp] = overrides[sp]
            continue
        if positional_idx < len(positional_colors):
            color_map[sp] = positional_colors[positional_idx]
            positional_idx += 1
            continue
        if sp in _DEFAULT_COLOR_SCHEME:
            color_map[sp] = mcolors.to_rgba(_DEFAULT_COLOR_SCHEME[sp])
        else:
            unknown.append(sp)

    if unknown:
        cmap = plt.get_cmap("tab20", len(unknown))
        for idx, sp in enumerate(unknown):
            color_map[sp] = cmap(idx)

    return color_map


def _parse_species_color_overrides(
    entries: list[str] | None,
) -> tuple[dict[str, tuple[float, float, float, float]], list[tuple[float, float, float, float]]]:
    """Parse --species-color entries.

    Supports two formats (may be mixed):
      - ``SPECIES=COLOR``  — named override for a specific species
      - ``COLOR[,COLOR…]`` — positional colors assigned to species in encounter order
    """
    if not entries:
        return {}, []
    overrides: dict[str, tuple[float, float, float, float]] = {}
    positional: list[tuple[float, float, float, float]] = []
    for entry in entries:
        if "=" not in entry:
            # Treat as comma-separated positional color(s)
            for part in entry.split(","):
                part = part.strip()
                if not part:
                    continue
                hex_part = part if part.startswith("#") else f"#{part}"
                positional.append(mcolors.to_rgba(hex_part))
        else:
            species, color = entry.split("=", 1)
            species = species.strip()
            color = color.strip()
            if not species or not color:
                raise ValueError(f"--species-color must look like SPECIES=COLOR (got {entry!r})")
            overrides[species] = mcolors.to_rgba(color)
    return overrides, positional


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


def _read_domain_z_range(h5, domain_index: int = 0) -> tuple[float, float] | None:
    """Read domain z-range from geometry as (z_min, z_max)."""
    base = f"/domains/domain_{domain_index}/geometry"
    origin_key = f"{base}/origin_m"
    length_key = f"{base}/length_m"
    if origin_key not in h5 or length_key not in h5:
        return None
    origin = np.asarray(h5[origin_key][()], dtype=float)
    if origin.shape[0] < 3:
        return None
    z0 = float(origin[2])
    length = float(h5[length_key][()])
    return (z0, z0 + length)


def _read_domain_geometry(h5, domain_index: int = 0) -> dict[str, float] | None:
    base = f"/domains/domain_{domain_index}/geometry"
    origin_key = f"{base}/origin_m"
    if origin_key not in h5:
        return None
    origin = np.asarray(h5[origin_key][()], dtype=float)
    if origin.shape[0] < 3:
        return None
    length = None
    length_key = f"{base}/length_m"
    if length_key in h5:
        length = float(h5[length_key][()])
    radius = None
    for key in (f"{base}/radius_m", f"{base}/radius_out_m"):
        if key in h5:
            radius = float(h5[key][()])
            break
    return {
        "x0": float(origin[0]),
        "y0": float(origin[1]),
        "z0": float(origin[2]),
        "length_m": float(length) if length is not None else np.nan,
        "radius_m": float(radius) if radius is not None else np.nan,
    }


def _axis_padding(mins: np.ndarray, maxs: np.ndarray, padding_frac: float) -> np.ndarray:
    frac = max(0.0, float(padding_frac))
    span = np.maximum(maxs - mins, 0.0)
    pad = span * frac
    # Avoid singular axes when all points have same coordinate.
    for idx in range(2):
        if span[idx] <= 0.0:
            pad[idx] = max(1e-12, abs(mins[idx]) * 1e-6)
    return pad


def _axis_domain_extent(
    axis_idx: int,
    geometry: dict[str, float],
    include_radial: bool,
) -> tuple[float, float] | None:
    if axis_idx == 2:
        length = float(geometry["length_m"])
        if np.isfinite(length):
            z0 = float(geometry["z0"])
            return (z0, z0 + length)
        return None
    if include_radial and axis_idx in (0, 1):
        radius = float(geometry["radius_m"])
        if np.isfinite(radius):
            center = float(geometry["x0"] if axis_idx == 0 else geometry["y0"])
            return (center - radius, center + radius)
    return None


def _merge_axis_limits_with_domain(
    mins: np.ndarray,
    maxs: np.ndarray,
    x_idx: int,
    y_idx: int,
    scales: np.ndarray,
    geometry: dict[str, float],
    include_radial: bool,
    force_domain_limits: bool,
) -> tuple[np.ndarray, np.ndarray]:
    mins_out = mins.copy()
    maxs_out = maxs.copy()

    for arr_idx, axis_idx in enumerate((x_idx, y_idx)):
        extent = _axis_domain_extent(axis_idx, geometry, include_radial=include_radial)
        if extent is None:
            continue
        lo, hi = sorted(extent)
        scaled = np.array([lo, hi], dtype=float) * float(scales[arr_idx])
        ext_min = float(np.min(scaled))
        ext_max = float(np.max(scaled))
        if force_domain_limits:
            mins_out[arr_idx] = ext_min
            maxs_out[arr_idx] = ext_max
        else:
            mins_out[arr_idx] = min(mins_out[arr_idx], ext_min)
            maxs_out[arr_idx] = max(maxs_out[arr_idx], ext_max)
    return mins_out, maxs_out


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

    if args.legend == "off":
        return
    if args.legend == "auto" and len(color_map) > 8:
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
        loc=args.legend_loc,
        framealpha=0.85,
        fontsize=args.legend_fontsize,
    )
    if theme:
        legend.get_frame().set_facecolor(theme["legend_face"])
        legend.get_frame().set_edgecolor(theme["text"])
    for text in legend.get_texts():
        if theme:
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

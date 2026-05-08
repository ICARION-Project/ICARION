#!/usr/bin/env python3
"""
Plot arrival-time distributions per species with optional Gaussian peak fits.

This script focuses on arrival-time distributions only (no mobility calculation).
In TIMS mode, the CSV export is extended with ejection voltage values read from H5.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

# Allow running as a standalone script without installing the package
if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.common import (
    ensure_nonempty_selection,
    find_arrival_times,
    fit_gaussian_histogram,
    gaussian,
    histogram_bin_edges,
    normalize_species_filter,
    open_trajectory,
    read_positions_subset,
    read_domain_geometry,
    select_ions_from_trajectory,
    species_color_map,
    species_for_indices,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot arrival-time distributions per species with Gaussian fits.")
    p.add_argument("--traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", type=Path, default=Path("analysis/output/arrival_time_dist.png"), help="Output plot path.")
    p.add_argument(
        "--out-csv",
        "--csv-out",
        dest="out_csv",
        type=Path,
        default=Path("analysis/output/arrival_time_dist.csv"),
        help="CSV output path.",
    )
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--max-ions", type=int, default=800, help="Cap on ions after filtering.")
    p.add_argument("--max-per-species", type=int, default=800, help="Cap per species before global cap.")
    p.add_argument("--rng-seed", type=int, default=0, help="Seed for ion subsampling.")
    p.add_argument("--bins", type=int, default=50, help="Bins for arrival-time histogram.")
    p.add_argument("--domain-index", type=int, default=0, help="Domain index for geometry.")
    p.add_argument("--tol-frac", type=float, default=0.02, help="Tolerance as fraction of drift length/radius.")
    p.add_argument("--arrival-z-m", type=float, default=None, help="Optional absolute z-position [m] of a custom arrival plane.")
    p.add_argument("--arrival-z-mm", type=float, default=None, help="Optional absolute z-position [mm] of a custom arrival plane.")
    p.add_argument("--time-stride", type=int, default=1, help="Stride when scanning positions for arrival (1 = full).")
    p.add_argument("--per-species", dest="per_species", action="store_true", default=True, help="Plot per-species breakdown (enabled by default).")
    p.add_argument("--no-per-species", dest="per_species", action="store_false", help="Disable per-species breakdown plots.")
    p.add_argument("--out-per-species", type=Path, default=None, help="Custom output path for per-species plot.")
    p.add_argument("--fit-gaussian", action="store_true", default=True, help="Fit Gaussian to each species (enabled by default).")
    p.add_argument("--tims-mode", action="store_true", help="Write ejection voltage in CSV (TIMS mode).")
    p.add_argument("--tims-voltage-path", type=str, default=None, help="Optional explicit H5 dataset path for TIMS ejection voltage.")
    return p.parse_args()


def _normalize_plot_path(path: Path) -> Path:
    # If user passes path without extension, default to PNG.
    if path.suffix:
        return path
    return path.with_suffix(".png")


def main() -> int:
    args = parse_args()
    if not args.traj.exists():
        raise FileNotFoundError(f"Trajectory file not found: {args.traj}")
    species_filter = normalize_species_filter(args.species)

    out_plot = _normalize_plot_path(args.out)

    with open_trajectory(args.traj) as h5:
        traj = h5["trajectory"]
        total_time = float(traj["time"][-1]) if "time" in traj else float(np.nan)
        species_ids, ion_indices = select_ions_from_trajectory(
            traj,
            species_filter=species_filter,
            max_ions=args.max_ions,
            max_per_species=args.max_per_species,
            rng_seed=args.rng_seed,
        )
        ensure_nonempty_selection(ion_indices)

        length_m, radius_m, z_origin_m = read_domain_geometry(h5, args.domain_index, traj, ion_indices)

        times, positions = read_positions_subset(
            traj,
            ion_indices=ion_indices,
            time_stride=max(1, args.time_stride),
            max_frames=None,
        )
        arrival_plane_override_m = args.arrival_z_m
        if args.arrival_z_mm is not None:
            arrival_plane_override_m = args.arrival_z_mm * 1e-3

        arrivals, mask_valid, _arrival_plane_z = find_arrival_times(
            times=times,
            positions=positions,
            z_origin_m=z_origin_m,
            length_m=length_m,
            radius_m=radius_m,
            tol_frac=args.tol_frac,
            arrival_plane_z_m=arrival_plane_override_m,
        )
        if not np.any(mask_valid):
            raise RuntimeError("No ions reached the end of the drift region.")

        species_selected = species_for_indices(species_ids, ion_indices)
        species_valid = species_selected[mask_valid]

        ejection_voltage_valid: Optional[np.ndarray] = None
        ejection_voltage_source: Optional[str] = None
        if args.tims_mode:
            ejection_voltage_valid, ejection_voltage_source = _extract_tims_ejection_voltage(
                h5,
                domain_idx=args.domain_index,
                arrival_times=arrivals,
                trajectory_times=np.asarray(traj["time"], dtype=float),
                explicit_path=args.tims_voltage_path,
            )

    _plot(
        arrivals,
        species_valid,
        out_plot,
        bins=args.bins,
        fit_gaussian=args.fit_gaussian,
        total_time=total_time,
        num_ions=len(arrivals),
    )

    if args.per_species:
        per_species_path = (
            _normalize_plot_path(args.out_per_species)
            if args.out_per_species is not None
            else out_plot.with_name(f"{out_plot.stem}_per_species{out_plot.suffix}")
        )
        _plot_per_species(
            arrivals,
            species_valid,
            per_species_path,
            bins=args.bins,
            fit_gaussian=args.fit_gaussian,
        )
        print(f"Wrote {per_species_path}")

    _write_csv(
        arrivals,
        species_valid,
        args.out_csv,
        ejection_voltage=ejection_voltage_valid,
        ejection_voltage_source=ejection_voltage_source,
    )

    if args.tims_mode:
        if ejection_voltage_source is None:
            print("TIMS mode: no ejection-voltage dataset found; CSV values are NaN")
        else:
            print(f"TIMS mode: ejection voltage source = {ejection_voltage_source}")

    print(f"Wrote {out_plot}")
    print(f"Wrote {args.out_csv}")
    return 0


def _extract_tims_ejection_voltage(
    h5,
    domain_idx: int,
    arrival_times: np.ndarray,
    trajectory_times: np.ndarray,
    explicit_path: Optional[str],
) -> Tuple[np.ndarray, Optional[str]]:
    """
    Extract ejection voltage for TIMS mode.

    Supported dataset formats:
      - scalar: same voltage for all ions
      - vector shape (T,): interpolated at per-ion arrival time
      - vector shape (N_arrivals,): used directly
      - matrix shape (T,1): interpolated at per-ion arrival time
    """

    def _to_arrival_values(arr: np.ndarray) -> Optional[np.ndarray]:
        arr = np.asarray(arr, dtype=float)
        if arr.ndim == 0:
            return np.full(arrival_times.shape, float(arr), dtype=float)
        if arr.ndim == 1:
            if arr.size == trajectory_times.size:
                return np.interp(arrival_times, trajectory_times, arr)
            if arr.size == arrival_times.size:
                return arr.astype(float)
        if arr.ndim == 2 and arr.shape[1] == 1:
            vec = arr[:, 0]
            if vec.size == trajectory_times.size:
                return np.interp(arrival_times, trajectory_times, vec)
        return None

    candidate_paths = []
    if explicit_path:
        candidate_paths.append(explicit_path)
    candidate_paths += [
        f"/domains/domain_{domain_idx}/fields/tims/ejection_voltage_V",
        f"/domains/domain_{domain_idx}/fields/tims/elution_voltage_V",
        f"/domains/domain_{domain_idx}/fields/dc/ejection_voltage_V",
        f"/domains/domain_{domain_idx}/fields/DC/ejection_voltage_V",
        f"/domains/domain_{domain_idx}/fields/dc/elution_voltage_V",
        f"/domains/domain_{domain_idx}/fields/DC/elution_voltage_V",
        "/tims/ejection_voltage_V",
        "/tims/elution_voltage_V",
    ]

    seen = set()
    for path in candidate_paths:
        if path in seen:
            continue
        seen.add(path)
        if path not in h5:
            continue
        vals = _to_arrival_values(np.asarray(h5[path]))
        if vals is not None:
            return vals, path

    cfg_vals = _extract_tims_voltage_from_config_json(
        h5=h5,
        domain_idx=domain_idx,
        arrival_times=arrival_times,
    )
    if cfg_vals is not None:
        return cfg_vals

    fallback_paths = [
        f"/domains/domain_{domain_idx}/fields/dc/axial_V",
        f"/domains/domain_{domain_idx}/fields/DC/axial_V",
        f"/domains/domain_{domain_idx}/fields/dc/voltage_V",
        f"/domains/domain_{domain_idx}/fields/DC/voltage_V",
    ]
    for path in fallback_paths:
        if path not in h5:
            continue
        vals = _to_arrival_values(np.asarray(h5[path]))
        if vals is not None:
            return vals, path

    return np.full(arrival_times.shape, np.nan, dtype=float), None


def _extract_tims_voltage_from_config_json(
    h5,
    domain_idx: int,
    arrival_times: np.ndarray,
) -> Optional[Tuple[np.ndarray, str]]:
    cfg_path = "/metadata/config/config_json"
    if cfg_path not in h5:
        return None

    cfg_raw = h5[cfg_path][()]
    if isinstance(cfg_raw, bytes):
        cfg_text = cfg_raw.decode("utf-8")
    else:
        cfg_text = str(cfg_raw)
    try:
        cfg = json.loads(cfg_text)
    except Exception:
        return None

    domains = cfg.get("domains")
    if not isinstance(domains, list) or not (0 <= domain_idx < len(domains)):
        return None
    domain = domains[domain_idx]
    if not isinstance(domain, dict):
        return None

    fields = domain.get("fields", {})
    if not isinstance(fields, dict):
        return None
    tims = fields.get("TIMS") or fields.get("tims")
    if not isinstance(tims, dict) or not bool(tims.get("enabled", True)):
        return None

    z_nodes = np.asarray(tims.get("z_positions_m", []), dtype=float)
    ez_accum_profile = np.asarray(tims.get("Ez_accum_profile_V_m", []), dtype=float)
    ez_delta_profile = np.asarray(tims.get("Ez_delta_profile_V_m", []), dtype=float)

    length_m = None
    length_path = f"/domains/domain_{domain_idx}/geometry/length_m"
    if length_path in h5:
        length_m = float(h5[length_path][()])
    if length_m is None:
        geom = domain.get("geometry", {})
        if isinstance(geom, dict) and "length_m" in geom:
            length_m = float(geom["length_m"])
    if length_m is None:
        if z_nodes.size >= 2:
            length_m = float(z_nodes[-1] - z_nodes[0])
        else:
            length_m = 0.0

    accum_uniform = float(tims.get("Ez_accum_uniform_V_m", 0.0))
    delta_uniform = float(tims.get("Ez_delta_uniform_V_m", 0.0))

    def _integrated_voltage(values_profile: np.ndarray, uniform_fallback: float) -> float:
        # TIMSAxialFieldModel defines axial electric field Ez(z,t). Convert field to
        # potential difference along +z via ΔV = -∫ Ez dz.
        if z_nodes.size >= 2 and values_profile.size == z_nodes.size:
            trapz_fn = getattr(np, "trapezoid", None)
            if trapz_fn is None:
                trapz_fn = np.trapz
            return -float(trapz_fn(values_profile, z_nodes))
        return -float(uniform_fallback) * float(length_m)

    v_accum = _integrated_voltage(ez_accum_profile, accum_uniform)
    v_delta = _integrated_voltage(ez_delta_profile, delta_uniform)
    frac = _tims_ramp_fraction_for_times(arrival_times, tims)
    values = v_accum - frac * v_delta
    source = f"{cfg_path}:domains[{domain_idx}].fields.TIMS"
    return values.astype(float), source


def _tims_ramp_fraction_for_times(times_s: np.ndarray, tims: dict) -> np.ndarray:
    times = np.asarray(times_s, dtype=float)
    raw_fraction = tims.get("ramp_fraction")
    if isinstance(raw_fraction, (int, float)):
        return np.full(times.shape, float(np.clip(raw_fraction, 0.0, 1.0)), dtype=float)

    t0 = float(tims.get("ramp_start_s", 0.0))
    t1 = float(tims.get("ramp_end_s", 0.0))
    duration = t1 - t0
    if duration <= 0.0:
        return np.ones(times.shape, dtype=float)

    frac = np.empty(times.shape, dtype=float)
    before = times <= t0
    after = times >= t1
    mid = ~(before | after)
    frac[before] = 0.0
    frac[after] = 1.0

    if np.any(mid):
        rel = times[mid] - t0
        if str(tims.get("ramp_mode", "linear")).lower() == "exponential":
            tau = max(float(tims.get("ramp_tau_s", 0.0)), 1e-30)
            denom = 1.0 - np.exp(-duration / tau)
            if abs(denom) < 1e-15:
                frac[mid] = rel / duration
            else:
                frac[mid] = (1.0 - np.exp(-rel / tau)) / denom
        else:
            frac[mid] = rel / duration

    return np.clip(frac, 0.0, 1.0)
def _plot(
    arrivals: np.ndarray,
    species: np.ndarray,
    out_path: Path,
    bins: int,
    fit_gaussian: bool = True,
    total_time: float = np.nan,
    num_ions: int = 0,
):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(10, 6))

    edges = histogram_bin_edges(arrivals, bins)

    ax.hist(arrivals, bins=edges, color="#4c78a8", alpha=0.7, edgecolor="black", linewidth=1.0, label="Observed")

    if fit_gaussian:
        fit_params, _, _hist = fit_gaussian_histogram(arrivals, bins)
        if fit_params is not None:
            mu = fit_params["mu"]
            sigma = fit_params["sigma"]
            x_fit = np.linspace(edges[0], edges[-1], 200)
            y_fit = gaussian(x_fit, fit_params["a"], mu, sigma)
            ax.plot(x_fit, y_fit, "r-", linewidth=2.0, label=f"Gaussian fit (mu={mu:.3e}s, sigma={sigma:.3e}s)")

    mean_t = arrivals.mean()
    median_t = np.median(arrivals)
    ax.axvline(mean_t, color="#e45756", linestyle="--", linewidth=2.0, label=f"Mean: {mean_t:.3e}s")
    ax.axvline(median_t, color="#54a24b", linestyle=":", linewidth=2.0, label=f"Median: {median_t:.3e}s")

    ax.set_xlabel("Arrival time [s]", fontsize=11, fontweight="bold")
    ax.set_ylabel("Ion count", fontsize=11, fontweight="bold")
    ax.set_title("Arrival-time distribution", fontsize=13, fontweight="bold")
    ax.grid(True, alpha=0.2, linestyle="-", linewidth=0.5)
    ax.legend(loc="upper right", fontsize=10, framealpha=0.95)

    unique_species = len(set(species))
    title_extra = f"{unique_species} species, {num_ions} ions"
    if not np.isnan(total_time):
        title_extra += f", sim T={total_time:.3e}s"
    fig.suptitle(title_extra, fontsize=10, style="italic", y=0.98)

    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _plot_per_species(
    arrivals: np.ndarray,
    species: np.ndarray,
    out_path: Path,
    bins: int,
    fit_gaussian: bool = True,
):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    unique_species = sorted(set(species))
    n = len(unique_species)

    if n == 0:
        return

    ncols = min(2, n)
    nrows = (n + ncols - 1) // ncols
    fig, axes = plt.subplots(nrows, ncols, figsize=(12, 4 * nrows), squeeze=False)
    axes = axes.flatten()

    color_map = species_color_map(unique_species)

    for idx, sp in enumerate(unique_species):
        mask = species == sp
        t_sp = arrivals[mask]
        ax = axes[idx]

        edges = histogram_bin_edges(t_sp, bins)

        ax.hist(t_sp, bins=edges, color=color_map[sp], alpha=0.7, edgecolor="black", linewidth=0.8)

        if fit_gaussian:
            fit_params, _, _hist = fit_gaussian_histogram(t_sp, bins)
            if fit_params is not None:
                mu = fit_params["mu"]
                sigma = fit_params["sigma"]
                x_fit = np.linspace(edges[0], edges[-1], 200)
                y_fit = gaussian(x_fit, fit_params["a"], mu, sigma)
                ax.plot(x_fit, y_fit, "r-", linewidth=2.0, alpha=0.8)

        mean_t = t_sp.mean()
        median_t = np.median(t_sp)
        ax.axvline(mean_t, color="#e45756", linestyle="--", linewidth=1.5, alpha=0.9)
        ax.axvline(median_t, color="#54a24b", linestyle=":", linewidth=1.5, alpha=0.9)

        ax.set_xlabel("Arrival time [s]", fontsize=10)
        ax.set_ylabel("Count", fontsize=10)
        ax.set_title(f"{sp} (n={len(t_sp)}, mu={mean_t:.3e}s)", fontsize=11, fontweight="bold")
        ax.grid(True, alpha=0.2, linestyle="-", linewidth=0.5)

    for idx in range(n, len(axes)):
        axes[idx].axis("off")

    fig.suptitle("Arrival-time distributions per species", fontsize=13, fontweight="bold")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def _write_csv(
    arrivals: np.ndarray,
    species: np.ndarray,
    out_path: Path,
    ejection_voltage: Optional[np.ndarray] = None,
    ejection_voltage_source: Optional[str] = None,
):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    use_tims_cols = ejection_voltage is not None

    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        if use_tims_cols:
            writer.writerow(["ion_index", "species", "arrival_time_s", "ejection_voltage_V", "ejection_voltage_source"])
            for i, (t, sp, v) in enumerate(zip(arrivals, species, ejection_voltage)):
                src = ejection_voltage_source if ejection_voltage_source is not None else ""
                writer.writerow([i, sp, f"{float(t):.6e}", f"{float(v):.6e}", src])
        else:
            writer.writerow(["ion_index", "species", "arrival_time_s"])
            for i, (t, sp) in enumerate(zip(arrivals, species)):
                writer.writerow([i, sp, f"{float(t):.6e}"])


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Histogramme für Eliminationszeiten (radial/axial) aus ICARION HDF5-Trajektorien.
Unterstützt Orbitrap (inner/outer) und generische Zylinder.
"""

from __future__ import annotations

import argparse
from pathlib import Path
from typing import Dict, Iterable, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np

# Allow running as a standalone script without installing the package
if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.common import load_species_ids, open_trajectory, select_ion_indices


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Plot elimination-time histograms (radial/axial).")
    p.add_argument("--traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", type=Path, default=Path("analysis/output/elimination_hist.png"), help="Output image path.")
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--max-ions", type=int, default=500, help="Max ions to include after filtering.")
    p.add_argument("--max-per-species", type=int, default=200, help="Cap per species before global cap.")
    p.add_argument("--rng-seed", type=int, default=0, help="Seed for ion subsampling.")
    p.add_argument("--bins", type=int, default=60, help="Number of histogram bins.")
    p.add_argument("--log-bins", action="store_true", help="Use logarithmic time bins (requires t>0).")
    p.add_argument("--tol-frac", type=float, default=0.05, help="Boundary tolerance as fraction of radius/length.")
    p.add_argument(
        "--orbitrap",
        choices=["auto", "yes", "no"],
        default="auto",
        help="Orbitrap-specific radial classification (inner/outer). Default: auto-detect via geometry.",
    )
    p.add_argument("--domain-index", type=int, default=0, help="Domain index (default: 0).")
    p.add_argument(
        "--per-species",
        dest="per_species",
        action="store_true",
        default=True,
        help="Plot per-species breakdown (counts + time hist). Enabled by default.",
    )
    p.add_argument(
        "--no-per-species",
        dest="per_species",
        action="store_false",
        help="Disable per-species breakdown plots.",
    )
    p.add_argument(
        "--out-per-species",
        type=Path,
        default=None,
        help="Custom output path for per-species plot (default: <out> with _per_species suffix).",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    traj_path = args.traj
    if not traj_path.exists():
        raise FileNotFoundError(f"Trajectory file not found: {traj_path}")

    with open_trajectory(traj_path) as h5:
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
            raise RuntimeError("No ions selected; broaden filters or caps.")

        death_times = _read_death_times(h5, ion_indices)
        last_positions = traj["positions"][-1, ion_indices, :]  # Shape: (n, 3)
        total_time = float(traj["time"][-1]) if "time" in traj else float(np.nan)

        geom = _read_geometry(h5, args.domain_index, last_positions)
        orbitrap_mode = _orbitrap_mode(args.orbitrap, geom)

        reasons = classify_ions(last_positions, death_times, geom, args.tol_frac, orbitrap_mode)
        selected_species = np.array(species_ids)[ion_indices]
        plot_histograms(
            reasons,
            death_times,
            species=selected_species,
            bins=args.bins,
            log_bins=args.log_bins,
            out_path=args.out,
            total_time=total_time,
            orbitrap_mode=orbitrap_mode,
        )

        if args.per_species:
            per_species_path = (
                args.out_per_species
                if args.out_per_species is not None
                else args.out.with_name(f"{args.out.stem}_per_species{args.out.suffix}")
            )
            plot_per_species(
                reasons,
                death_times,
                species=selected_species,
                bins=args.bins,
                log_bins=args.log_bins,
                out_path=per_species_path,
                orbitrap_mode=orbitrap_mode,
            )
            print(f"Wrote {per_species_path}")

    print(f"Wrote {args.out}")
    return 0


def _read_death_times(h5, ion_indices) -> np.ndarray:
    if "ions" not in h5 or "death_time_s" not in h5["ions"]:
        raise KeyError("HDF5 file missing '/ions/death_time_s'")
    dt = h5["ions/death_time_s"][ion_indices]
    return np.asarray(dt, dtype=float)


class Geometry:
    def __init__(
        self,
        r_in: Optional[float],
        r_out: Optional[float],
        z_min: Optional[float],
        z_max: Optional[float],
    ):
        self.r_in = r_in
        self.r_out = r_out
        self.z_min = z_min
        self.z_max = z_max


def _read_geometry(h5, domain_idx: int, positions_last: np.ndarray) -> Geometry:
    dom_name = f"domain_{domain_idx}"
    geom_group = h5.get(f"domains/{dom_name}/geometry")

    def _scalar(name: str) -> Optional[float]:
        if geom_group and name in geom_group:
            return float(np.array(geom_group[name]).astype(float))
        return None

    origin = _scalar("origin_m")
    if origin is not None:
        try:
            origin_z = float(np.array(geom_group["origin_m"])[2])
        except Exception:
            origin_z = 0.0
    else:
        origin_z = 0.0

    length = _scalar("length_m")
    r_out = _scalar("radius_out_m") or _scalar("radius_m")
    r_in = _scalar("radius_in_m")

    if length is not None:
        z_min = origin_z
        z_max = origin_z + length
    else:
        # Fallback: infer rough bounds from last positions.
        z_min = float(positions_last[:, 2].min())
        z_max = float(positions_last[:, 2].max())

    return Geometry(r_in=r_in, r_out=r_out, z_min=z_min, z_max=z_max)


def _orbitrap_mode(flag: str, geom: Geometry) -> bool:
    if flag == "yes":
        return True
    if flag == "no":
        return False
    return geom.r_in is not None and geom.r_out is not None


def classify_ions(
    positions_last: np.ndarray,
    death_times: np.ndarray,
    geom: Geometry,
    tol_frac: float,
    orbitrap_mode: bool,
) -> np.ndarray:
    r = np.sqrt(positions_last[:, 0] ** 2 + positions_last[:, 1] ** 2)
    z = positions_last[:, 2]

    tol_r = (geom.r_out or np.nan) * tol_frac if geom.r_out else np.nan
    length = (geom.z_max - geom.z_min) if geom.z_min is not None and geom.z_max is not None else np.nan
    tol_z = length * tol_frac if not np.isnan(length) else tol_r

    reasons = np.full(r.shape, fill_value="alive", dtype=object)

    for i in range(len(r)):
        if death_times[i] < 0:
            reasons[i] = "alive"
            continue

        if orbitrap_mode:
            if geom.r_in is not None and r[i] <= geom.r_in + tol_r:
                reasons[i] = "inner"
                continue
            if geom.r_out is not None and r[i] >= geom.r_out - tol_r:
                reasons[i] = "outer"
                continue
        else:
            if geom.r_out is not None and r[i] >= geom.r_out - tol_r:
                reasons[i] = "radial"
                continue

        if geom.z_min is not None and z[i] <= geom.z_min - tol_z:
            reasons[i] = "axial_low"
            continue
        if geom.z_max is not None and z[i] >= geom.z_max + tol_z:
            reasons[i] = "axial_high"
            continue

        reasons[i] = "other"
    return reasons


def plot_histograms(
    reasons: np.ndarray,
    death_times: np.ndarray,
    species: np.ndarray,
    bins: int,
    log_bins: bool,
    out_path: Path,
    total_time: float,
    orbitrap_mode: bool,
):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Count plot (including alive).
    ax_counts = axes[0]
    order = ["alive"]
    order += ["inner", "outer"] if orbitrap_mode else ["radial"]
    order += ["axial_low", "axial_high", "other"]
    counts = {key: int(np.sum(reasons == key)) for key in order}
    ax_counts.bar(range(len(order)), [counts[k] for k in order], color="steelblue", alpha=0.8)
    ax_counts.set_xticks(range(len(order)), order, rotation=30, ha="right")
    ax_counts.set_ylabel("Ion count")
    ax_counts.set_title("Elimination channels")

    # Time histogram (only death_time >= 0).
    ax_hist = axes[1]
    mask_alive = death_times < 0
    t_elim = death_times[~mask_alive]
    reasons_elim = reasons[~mask_alive]

    if len(t_elim) > 0:
        t_min = t_elim[t_elim > 0].min() if np.any(t_elim > 0) else 1e-12
        t_max = t_elim.max()
        if log_bins and t_max > 0:
            edges = np.geomspace(t_min, t_max, num=bins + 1)
        else:
            edges = np.linspace(0, t_max, num=bins + 1)

        categories = ["inner", "outer"] if orbitrap_mode else ["radial"]
        categories += ["axial_low", "axial_high", "other"]
        colors: Dict[str, str] = {
            "inner": "#4c78a8",
            "outer": "#f58518",
            "radial": "#4c78a8",
            "axial_low": "#54a24b",
            "axial_high": "#e45756",
            "other": "#bab0ac",
        }
        bottoms = np.zeros_like(edges[:-1])
        for cat in categories:
            vals = t_elim[reasons_elim == cat]
            if len(vals) == 0:
                continue
            hist, _ = np.histogram(vals, bins=edges)
            ax_hist.bar(edges[:-1], hist, width=np.diff(edges), bottom=bottoms, align="edge", color=colors.get(cat, None), label=cat, alpha=0.8)
            bottoms += hist
        ax_hist.set_xlabel("Death time [s]")
        ax_hist.set_ylabel("Ion count")
        ax_hist.set_title("Elimination time histogram")
        if log_bins:
            ax_hist.set_xscale("log")
        ax_hist.legend()
    else:
        ax_hist.text(0.5, 0.5, "No eliminated ions", ha="center", va="center", transform=ax_hist.transAxes)
        ax_hist.set_axis_off()

    if not np.isnan(total_time):
        fig.suptitle(f"Elimination summary – {out_path.name} (sim T={total_time:.3e}s, ions={len(reasons)})")
    fig.tight_layout()
    fig.savefig(out_path, dpi=200)


def plot_per_species(
    reasons: np.ndarray,
    death_times: np.ndarray,
    species: np.ndarray,
    bins: int,
    log_bins: bool,
    out_path: Path,
    orbitrap_mode: bool,
):
    unique_species = sorted(set(species))
    categories = ["inner", "outer"] if orbitrap_mode else ["radial"]
    categories += ["axial_low", "axial_high", "other"]
    colors: Dict[str, str] = {
        "inner": "#4c78a8",
        "outer": "#f58518",
        "radial": "#4c78a8",
        "axial_low": "#54a24b",
        "axial_high": "#e45756",
        "other": "#bab0ac",
    }

    n = len(unique_species)
    if n == 0:
        return
    fig, axes = plt.subplots(n, 2, figsize=(12, 3 * n), squeeze=False)

    for row, sp in enumerate(unique_species):
        idx = species == sp
        sp_reasons = reasons[idx]
        sp_dt = death_times[idx]
        mask_alive = sp_dt < 0
        sp_elim = sp_dt[~mask_alive]
        sp_reasons_elim = sp_reasons[~mask_alive]

        # Counts per category
        ax_counts = axes[row, 0]
        counts = [int(np.sum(sp_reasons == cat)) for cat in (["alive"] + categories)]
        labels = ["alive"] + categories
        ax_counts.bar(range(len(labels)), counts, color=[colors.get(cat, "gray") for cat in labels], alpha=0.85)
        ax_counts.set_xticks(range(len(labels)), labels, rotation=30, ha="right")
        ax_counts.set_ylabel("Count")
        ax_counts.set_title(f"{sp} – elimination channels")

        # Time histogram per category
        ax_hist = axes[row, 1]
        if len(sp_elim) > 0:
            t_min = sp_elim[sp_elim > 0].min() if np.any(sp_elim > 0) else 1e-12
            t_max = sp_elim.max()
            if log_bins and t_max > 0:
                edges = np.geomspace(t_min, t_max, num=bins + 1)
            else:
                edges = np.linspace(0, t_max, num=bins + 1)

            bottoms = np.zeros_like(edges[:-1])
            for cat in categories:
                vals = sp_elim[sp_reasons_elim == cat]
                if len(vals) == 0:
                    continue
                hist, _ = np.histogram(vals, bins=edges)
                ax_hist.bar(
                    edges[:-1],
                    hist,
                    width=np.diff(edges),
                    bottom=bottoms,
                    align="edge",
                    color=colors.get(cat, None),
                    label=cat,
                    alpha=0.85,
                )
                bottoms += hist
            if log_bins:
                ax_hist.set_xscale("log")
            ax_hist.set_ylabel("Count")
            ax_hist.set_xlabel("Death time [s]")
            ax_hist.set_title(f"{sp} – elimination time")
        else:
            ax_hist.text(0.5, 0.5, "No eliminated ions", ha="center", va="center", transform=ax_hist.transAxes)
            ax_hist.set_axis_off()

    # Add legend to the last histogram axis
    axes[0, 1].legend(loc="upper right")
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)


if __name__ == "__main__":
    raise SystemExit(main())

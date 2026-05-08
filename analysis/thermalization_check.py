#!/usr/bin/env python3
"""
Quick thermalization check:
  a) Temperature vs. time from velocities compared to domain temperature
  b) Speed distribution averaged over the last N timesteps vs. Maxwell-Boltzmann at that temperature

Assumes one dominant species; if multiple species remain after filtering, it will use the first species' mass for MB overlay and emit a warning.
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
    ensure_nonempty_selection,
    maxwell_speed_pdf,
    normalize_species_filter,
    open_trajectory,
    read_domain_environment,
    read_species_mass_map,
    select_ions_from_trajectory,
    species_masses_for_selection,
    species_for_indices,
    temperature_from_velocities,
    temperature_series_from_velocities,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Thermalization temperature and velocity distribution check.")
    p.add_argument("--traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", type=Path, default=Path("analysis/output/thermalization.png"), help="Output plot path.")
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--max-ions", type=int, default=800, help="Cap on ions after filtering.")
    p.add_argument("--max-per-species", type=int, default=800, help="Cap per species before global cap.")
    p.add_argument("--rng-seed", type=int, default=0, help="Seed for ion subsampling.")
    p.add_argument("--time-stride", type=int, default=10, help="Read every Nth timestep for temperature curve.")
    p.add_argument("--window", type=int, default=10, help="Number of last timesteps to average for the speed histogram.")
    p.add_argument("--bins", type=int, default=70, help="Histogram bins for speed distribution.")
    p.add_argument("--domain-index", type=int, default=0, help="Domain index for target temperature.")
    p.add_argument(
        "--tol-frac",
        type=float,
        default=0.05,
        help="Relative tolerance for thermalization time estimate (|T-T_target|/T_target <= tol).",
    )
    p.add_argument(
        "--steady-window",
        type=int,
        default=5,
        help="Consecutive points within tolerance required to declare thermalized.",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    if not args.traj.exists():
        raise FileNotFoundError(f"Trajectory file not found: {args.traj}")
    species_filter = normalize_species_filter(args.species)

    with open_trajectory(args.traj) as h5:
        traj = h5["trajectory"]
        species_ids, ion_indices = select_ions_from_trajectory(
            traj,
            species_filter=species_filter,
            max_ions=args.max_ions,
            max_per_species=args.max_per_species,
            rng_seed=args.rng_seed,
        )
        ensure_nonempty_selection(ion_indices)

        species_selected = species_for_indices(species_ids, ion_indices)
        mass_map = read_species_mass_map(h5)
        masses = species_masses_for_selection(species_selected, mass_map)

        target_temp, _pressure = read_domain_environment(
            h5,
            args.domain_index,
            require_pressure=False,
            allow_missing=True,
        )

        time_ds = traj["time"]
        vel_ds = traj["velocities"]  # shape: (T, N, 3)

        # Temperature vs time (strided)
        stride = max(1, args.time_stride)
        time = time_ds[::stride]
        velocities = vel_ds[::stride, ion_indices, :]
        temp_series = temperature_series_from_velocities(velocities, masses)
        t_therm = estimate_thermalization_time(
            time=time,
            temp_series=temp_series,
            target_temp=target_temp,
            tol_frac=args.tol_frac,
            steady_window=args.steady_window,
        )

        # Last window for MB comparison (unstrided to get true last frames)
        window = max(1, min(args.window, vel_ds.shape[0]))
        vel_tail = vel_ds[-window:, ion_indices, :]
        speeds_tail = np.linalg.norm(vel_tail, axis=-1).reshape(-1)
        T_tail = temperature_from_velocities(vel_tail, masses)

    # MB distribution uses a single mass (best effort)
    unique_masses = np.unique(masses)
    mb_mass = unique_masses[0]
    if len(unique_masses) > 1:
        print(f"⚠️  Multiple species detected; using mass of first species ({mb_mass:.3e} kg) for MB overlay.")

    bins = args.bins
    counts, edges = np.histogram(speeds_tail, bins=bins, density=True)
    centers = 0.5 * (edges[:-1] + edges[1:])
    mb_tail = maxwell_speed_pdf(centers, T_tail, mb_mass)
    mb_target = maxwell_speed_pdf(centers, target_temp, mb_mass) if target_temp is not None else None

    _plot(
        time=time,
        temp_series=temp_series,
        target_temp=target_temp,
        t_therm=t_therm,
        centers=centers,
        counts=counts,
        mb_tail=mb_tail,
        mb_target=mb_target,
        T_tail=T_tail,
        out_path=args.out,
        n_ions=len(ion_indices),
        window=window,
    )
    print(f"Wrote {args.out}")
    return 0


def estimate_thermalization_time(
    time: np.ndarray,
    temp_series: np.ndarray,
    target_temp: float | None,
    tol_frac: float,
    steady_window: int,
) -> float | None:
    if target_temp is None or np.isnan(target_temp):
        return None
    tol = abs(target_temp) * tol_frac
    within = np.abs(temp_series - target_temp) <= tol
    w = max(1, steady_window)
    for i in range(len(time) - w + 1):
        if np.all(within[i : i + w]):
            return float(time[i])
    return None


def _plot(
    time: np.ndarray,
    temp_series: np.ndarray,
    target_temp: float | None,
    t_therm: float | None,
    centers: np.ndarray,
    counts: np.ndarray,
    mb_tail: np.ndarray,
    mb_target: np.ndarray | None,
    T_tail: float,
    out_path: Path,
    n_ions: int,
    window: int,
):
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    ax_t = axes[0]
    ax_t.plot(time, temp_series, label="Simulated T", color="#1f77b4")
    if target_temp is not None:
        ax_t.axhline(target_temp, color="#e45756", linestyle="--", label=f"Domain T={target_temp:.1f} K")
    if t_therm is not None:
        ax_t.axvline(t_therm, color="#54a24b", linestyle=":", label=f"T therm ~ {t_therm:.3e}s")
    ax_t.set_xlabel("Time [s]")
    ax_t.set_ylabel("Temperature [K]")
    ax_t.set_title("Temperature vs time")
    ax_t.grid(True, alpha=0.3)
    ax_t.legend()

    ax_mb = axes[1]
    ax_mb.bar(centers, counts, width=np.diff(np.hstack([centers, centers[-1]]))[0], alpha=0.6, label="Sim speeds (last window)")
    ax_mb.plot(centers, mb_tail, color="#1f77b4", label=f"MB @ T_sim={T_tail:.1f} K")
    if mb_target is not None:
        ax_mb.plot(centers, mb_target, color="#e45756", linestyle="--", label="MB @ domain T")
    ax_mb.set_xlabel("Speed [m/s]")
    ax_mb.set_ylabel("Density")
    ax_mb.set_title(f"Speed distribution (last {window} steps, {n_ions} ions)")
    ax_mb.grid(True, alpha=0.3)
    ax_mb.legend()

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=200)


if __name__ == "__main__":
    raise SystemExit(main())

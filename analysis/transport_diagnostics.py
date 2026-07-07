#!/usr/bin/env python3
"""Transport diagnostics for IMS/TIMS/TWIMS-style trajectory outputs."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.trajectory import (
    ensure_nonempty_selection,
    find_arrival_times,
    normalize_species_filter,
    open_trajectory,
    read_domain_environment,
    read_domain_geometry,
    read_positions_subset,
    select_ions_from_trajectory,
    species_for_indices,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute arrival, width, resolving-power, and drift-speed diagnostics.")
    p.add_argument("--traj", "--input", dest="traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", "--output", dest="out", type=Path, default=Path("analysis/output/transport_diagnostics.png"), help="Output plot path.")
    p.add_argument("--out-csv", type=Path, default=Path("analysis/output/transport_diagnostics.csv"), help="Output CSV path.")
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--domain-index", type=int, default=0, help="Domain index for geometry/environment.")
    p.add_argument("--arrival-z-m", type=float, default=None, help="Optional arrival plane [m].")
    p.add_argument("--tol-frac", type=float, default=0.02, help="Arrival tolerance fraction.")
    p.add_argument("--max-ions", type=int, default=2000, help="Maximum ions to analyze.")
    p.add_argument("--rng-seed", type=int, default=0, help="Subsampling seed.")
    p.add_argument("--bins", type=int, default=60, help="Arrival histogram bins.")
    return p.parse_args()


def _width_stats(values: np.ndarray) -> tuple[float, float, float]:
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        return float("nan"), float("nan"), float("nan")
    mean = float(np.mean(vals))
    sigma = float(np.std(vals))
    resolving_power = float(mean / (2.355 * sigma)) if sigma > 0.0 else float("inf")
    return mean, sigma, resolving_power


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
        length_m, radius_m, z_origin_m = read_domain_geometry(h5, args.domain_index, traj, ion_indices)
        temperature_K, pressure_Pa = read_domain_environment(h5, args.domain_index, require_pressure=False, allow_missing=True)
        time, positions = read_positions_subset(traj, ion_indices, time_stride=1, max_frames=None)

    arrivals, valid_mask, _ = find_arrival_times(
        time,
        positions,
        z_origin_m=z_origin_m,
        length_m=length_m,
        radius_m=radius_m,
        tol_frac=args.tol_frac,
        arrival_plane_z_m=args.arrival_z_m,
    )
    if arrivals.size == 0:
        raise RuntimeError("No arrivals detected.")
    species_valid = species_for_indices(species_ids, ion_indices)[valid_mask]

    rows = []
    for species in sorted(set(species_valid)):
        vals = arrivals[species_valid == species]
        mean, sigma, resolving_power = _width_stats(vals)
        drift_speed = length_m / mean if mean > 0.0 else float("nan")
        rows.append(
            {
                "species": species,
                "n_arrivals": int(vals.size),
                "arrival_mean_s": mean,
                "arrival_sigma_s": sigma,
                "resolving_power_t_over_fwhm": resolving_power,
                "drift_speed_m_s": drift_speed,
                "temperature_K": temperature_K,
                "pressure_Pa": pressure_Pa,
            }
        )

    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8, 4.8))
    for species in sorted(set(species_valid)):
        ax.hist(arrivals[species_valid == species], bins=args.bins, alpha=0.55, label=species)
    ax.set_xlabel("arrival time [s]")
    ax.set_ylabel("count")
    ax.grid(alpha=0.25)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(args.out, dpi=180)
    plt.close(fig)

    args.out_csv.parent.mkdir(parents=True, exist_ok=True)
    with args.out_csv.open("w", newline="", encoding="utf-8") as f:
        fieldnames = list(rows[0].keys())
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {args.out}")
    print(f"Wrote {args.out_csv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

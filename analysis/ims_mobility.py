#!/usr/bin/env python3
"""
Estimate ion mobility (K) and reduced mobility (K0) from IMS trajectories.

Workflow:
  - Select ions (species filter, subsampling)
  - Find arrival times at drift tube end (first timestep with z >= z_end - tol and within radial bounds)
  - Compute drift velocity v = L_eff / t_arrival, mobility K = v / E
    - Convert to reduced mobility K0 using the shared reference-state conversion helper
  - Plot arrival-time histogram with mean/median markers and save CSV with per-ion values

Inputs pulled from HDF5 when available:
  - Drift length: /domains/domain_<i>/geometry/length_m
  - Drift origin: /domains/domain_<i>/geometry/origin_m (for physically correct end-plane)
  - Electric field: derived from /domains/domain_<i>/fields/DC/* (voltage) divided by length, if present.
    Fallback: requires --field-V or --field-VM.
  - Temperature, pressure: /domains/domain_<i>/environment/temperature_K, pressure_Pa
"""

from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Optional

import matplotlib.pyplot as plt
import numpy as np

# Allow running as a standalone script without installing the package
if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.common import (
    ensure_nonempty_selection,
    field_strength_from_en_td,
    find_arrival_times,
    histogram_bin_edges,
    normalize_species_filter,
    open_trajectory,
    read_domain_environment,
    read_domain_geometry,
    read_positions_subset,
    reduced_mobility_from_mobility,
    select_ions_from_trajectory,
    species_for_indices,
)


@dataclass(frozen=True)
class FieldStrength:
    value_V_m: float
    source: str
    en_td: float | None = None
    voltage_V: float | None = None


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Compute mobility and reduced mobility from IMS trajectories.")
    p.add_argument("--traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--out", type=Path, default=Path("analysis/output/ims_mobility.png"), help="Output plot path.")
    p.add_argument("--out-csv", type=Path, default=Path("analysis/output/ims_mobility.csv"), help="CSV output path.")
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--max-ions", type=int, default=800, help="Cap on ions after filtering.")
    p.add_argument("--max-per-species", type=int, default=800, help="Cap per species before global cap.")
    p.add_argument("--rng-seed", type=int, default=0, help="Seed for ion subsampling.")
    p.add_argument("--bins", type=int, default=60, help="Bins for arrival-time histogram.")
    p.add_argument("--domain-index", type=int, default=0, help="Domain index for geometry/fields/env.")
    p.add_argument("--tol-frac", type=float, default=0.02, help="Tolerance as fraction of drift length/radius.")
    p.add_argument("--arrival-z-m", type=float, default=None, help="Optional absolute z-position [m] of a custom arrival plane.")
    p.add_argument("--arrival-z-mm", type=float, default=None, help="Optional absolute z-position [mm] of a custom arrival plane.")
    p.add_argument("--field-V", dest="field_V", type=float, default=None, help="Axial voltage (V) across drift tube.")
    p.add_argument("--field-VM", dest="field_Vm", type=float, default=None, help="Axial field strength (V/m). Overrides voltage/length if set.")
    p.add_argument("--time-stride", type=int, default=1, help="Stride when scanning positions for arrival (1 = full).")
    p.add_argument(
        "--length-mode",
        choices=("effective", "full", "path", "window20_80"),
        default="effective",
        help=(
            "Drift-length definition for mobility: "
            "'effective' uses per-ion distance from initial z to arrival plane (recommended), "
            "'full' uses domain length_m, "
            "'path' uses the actual trajectory arc length up to the arrival timestep, "
            "'window20_80' uses only the transit from 20%% to 80%% of the drift length."
        ),
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

        # Env parameters
        temp_K, pressure_Pa = read_domain_environment(h5, args.domain_index)
        length_m, radius_m, z_origin_m = read_domain_geometry(h5, args.domain_index, traj, ion_indices)
        field = _read_field_strength(h5, args.domain_index, length_m, args.field_V, args.field_Vm)
        field_Vm = field.value_V_m

        # Arrival detection
        times, positions = read_positions_subset(
            traj,
            ion_indices=ion_indices,
            time_stride=max(1, args.time_stride),
            max_frames=None,
        )
        arrival_plane_override_m = args.arrival_z_m
        if args.arrival_z_mm is not None:
            arrival_plane_override_m = args.arrival_z_mm * 1e-3

        arrivals, mask_valid, arrival_plane_z = find_arrival_times(
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

        z_start = positions[0, :, 2]
        effective_length_all = np.maximum(1e-30, arrival_plane_z - z_start)
        effective_length = effective_length_all[mask_valid]

        if args.length_mode == "window20_80":
            window_times, window_mask, window_length = compute_fractional_window_mobility_inputs(
                times=np.asarray(times, dtype=float),
                positions=np.asarray(positions, dtype=float),
                z_origin_m=z_origin_m,
                length_m=length_m,
                radius_m=radius_m,
                tol_frac=args.tol_frac,
                z_fraction_start=0.2,
                z_fraction_end=0.8,
            )
            if not np.any(window_mask):
                raise RuntimeError("No ions crossed both the 20% and 80% drift-length planes.")
            arrivals = window_times
            mask_valid = window_mask
            length_for_mobility = window_length
        elif args.length_mode == "effective":
            length_for_mobility = effective_length
        elif args.length_mode == "path":
            length_for_mobility = compute_path_lengths_to_arrival(
                positions=positions,
                mask_valid=mask_valid,
                arrival_times=arrivals,
                sampled_times=np.asarray(times, dtype=float),
            )
        else:
            length_for_mobility = np.full_like(arrivals, length_m, dtype=float)
        drift_vel = length_for_mobility / arrivals  # m/s
        mobility = drift_vel / field_Vm  # m2/Vs
        K0 = reduced_mobility_from_mobility(mobility, pressure_Pa, temp_K)

        species_selected = species_for_indices(species_ids, ion_indices)
        species_valid = species_selected[mask_valid]

    # Plot and CSV
    _plot(
        arrivals,
        mobility,
        K0,
        species_valid,
        args.out,
        bins=args.bins,
        field_Vm=field_Vm,
        temp=temp_K,
        pressure=pressure_Pa,
        length=length_m,
        length_mode=args.length_mode,
        effective_length=effective_length,
        arrival_plane_z=arrival_plane_z,
        field_source=field.source,
    )
    _write_csv(arrivals, mobility, K0, species_valid, args.out_csv, length_for_mobility, field)
    print(f"Wrote {args.out}")
    print(f"Wrote {args.out_csv}")
    return 0
def _read_field_strength(
    h5,
    domain_idx: int,
    length_m: float,
    field_V: Optional[float],
    field_Vm_override: Optional[float],
) -> FieldStrength:
    field = _field_strength_details(h5, domain_idx, length_m, field_V, field_Vm_override)
    if field.source.startswith("manual"):
        print(f"Warning: using manual field override ({field.source}) = {field.value_V_m:.6g} V/m")
    return field


def _field_strength_details(
    h5,
    domain_idx: int,
    length_m: float,
    field_V: Optional[float],
    field_Vm_override: Optional[float],
) -> FieldStrength:
    if field_Vm_override is not None:
        return FieldStrength(float(field_Vm_override), "manual_field_Vm")

    en_paths = [
        f"/domains/domain_{domain_idx}/fields/dc/EN_Td",
        f"/domains/domain_{domain_idx}/fields/DC/EN_Td",
    ]
    for path in en_paths:
        if path in h5:
            en_td = float(h5[path][()])
            temp, pressure = read_domain_environment(h5, domain_idx)
            return FieldStrength(
                field_strength_from_en_td(en_td, pressure, temp),
                f"hdf5:{path}",
                en_td=en_td,
            )

    dc_paths = [
        f"/domains/domain_{domain_idx}/fields/dc/axial_V",
        f"/domains/domain_{domain_idx}/fields/DC/axial_V",
        f"/domains/domain_{domain_idx}/fields/dc/voltage_V",
        f"/domains/domain_{domain_idx}/fields/DC/voltage_V",
        f"/domains/domain_{domain_idx}/fields/dc/drift_V",
        f"/domains/domain_{domain_idx}/fields/DC/drift_V",
    ]
    voltage = None
    for path in dc_paths:
        if path in h5:
            voltage = float(h5[path][()])
            return FieldStrength(voltage / length_m, f"hdf5:{path}/length_m", voltage_V=voltage)

    if field_V is None:
        raise KeyError("No EN_Td or DC voltage found in HDF5; provide --field-V or --field-VM.")
    voltage = float(field_V)
    return FieldStrength(voltage / length_m, "manual_field_V_over_length", voltage_V=voltage)


def compute_fractional_window_mobility_inputs(
    times: np.ndarray,
    positions: np.ndarray,
    z_origin_m: float,
    length_m: float,
    radius_m: float,
    tol_frac: float,
    z_fraction_start: float,
    z_fraction_end: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    z_plane_start = z_origin_m + z_fraction_start * length_m
    z_plane_end = z_origin_m + z_fraction_end * length_m
    t_start, mask_start = _first_plane_crossing_times(
        times=times,
        positions=positions,
        z_plane_m=z_plane_start,
        radius_m=radius_m,
        tol_frac=tol_frac,
    )
    t_end, mask_end = _first_plane_crossing_times(
        times=times,
        positions=positions,
        z_plane_m=z_plane_end,
        radius_m=radius_m,
        tol_frac=tol_frac,
    )
    mask_valid = mask_start & mask_end & (t_end > t_start)
    window_time = t_end[mask_valid] - t_start[mask_valid]
    window_length = np.full(window_time.shape, max((z_fraction_end - z_fraction_start) * length_m, 1e-30), dtype=float)
    return window_time, mask_valid, window_length


def _first_plane_crossing_times(
    times: np.ndarray,
    positions: np.ndarray,
    z_plane_m: float,
    radius_m: float,
    tol_frac: float,
) -> tuple[np.ndarray, np.ndarray]:
    z = positions[:, :, 2]
    r = np.sqrt(positions[:, :, 0] ** 2 + positions[:, :, 1] ** 2)
    tol_r = radius_m * tol_frac
    reached = (z >= z_plane_m) & (r <= radius_m + tol_r)

    crossing_times = np.full(z.shape[1], np.nan, dtype=float)
    for ion in range(z.shape[1]):
        if not np.any(reached[:, ion]):
            continue
        idx = int(np.argmax(reached[:, ion]))
        if idx == 0:
            crossing_times[ion] = float(times[idx])
            continue
        z0 = float(z[idx - 1, ion])
        z1 = float(z[idx, ion])
        t0 = float(times[idx - 1])
        t1 = float(times[idx])
        if z1 <= z0:
            crossing_times[ion] = t1
            continue
        frac = np.clip((z_plane_m - z0) / max(z1 - z0, 1e-30), 0.0, 1.0)
        crossing_times[ion] = t0 + frac * (t1 - t0)
    mask_valid = ~np.isnan(crossing_times)
    return crossing_times, mask_valid


def compute_path_lengths_to_arrival(
    positions: np.ndarray,
    mask_valid: np.ndarray,
    arrival_times: np.ndarray,
    sampled_times: np.ndarray,
) -> np.ndarray:
    """Compute per-ion trajectory arc length up to the sampled arrival timestep."""
    valid_indices = np.flatnonzero(mask_valid)
    lengths = np.full(arrival_times.shape, np.nan, dtype=float)
    for out_idx, ion_idx in enumerate(valid_indices):
        arrival_t = arrival_times[out_idx]
        time_idx = int(np.searchsorted(sampled_times, arrival_t, side="left"))
        time_idx = min(max(time_idx, 0), positions.shape[0] - 1)
        coords = positions[: time_idx + 1, ion_idx, :]
        if coords.shape[0] < 2:
            lengths[out_idx] = 0.0
            continue
        step_lengths = np.linalg.norm(np.diff(coords, axis=0), axis=1)
        lengths[out_idx] = float(np.sum(step_lengths))
    return lengths


def _plot(
    arrivals,
    mobility,
    K0,
    species,
    out_path: Path,
    bins: int,
    field_Vm: float,
    temp: float,
    pressure: float,
    length: float,
    length_mode: str,
    effective_length: np.ndarray,
    arrival_plane_z: float,
    field_source: str,
):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Arrival histogram
    ax_arr = axes[0]
    edges = histogram_bin_edges(arrivals, bins)
    ax_arr.hist(arrivals, bins=edges, color="#4c78a8", alpha=0.8)
    ax_arr.axvline(np.median(arrivals), color="#e45756", linestyle="--", label=f"Median t={np.median(arrivals):.3e}s")
    ax_arr.set_xlabel("Arrival time [s]")
    ax_arr.set_ylabel("Count")
    ax_arr.set_title("Arrival-time distribution")
    ax_arr.legend()
    ax_arr.grid(True, alpha=0.2)

    # Mobility scatter
    ax_mob = axes[1]
    ax_mob.scatter(arrivals, mobility, alpha=0.7, s=12, label="K per ion")
    ax_mob.axhline(
        np.median(mobility),
        color="#e45756",
        linestyle="--",
        label=f"Median K={np.median(mobility):.3e} m²/Vs",
    )
    ax_mob.set_xlabel("Arrival time [s]")
    ax_mob.set_ylabel("Mobility K [m²/(V·s)]")
    ax_mob.set_title(f"K, K0 stats (N={len(arrivals)})")
    ax_mob.grid(True, alpha=0.2)
    ax_mob.legend()

    K0_cm2 = K0 * 1e4
    if effective_length.size:
        leff_med = np.median(effective_length)
        leff_mean = np.mean(effective_length)
    else:
        leff_med = np.nan
        leff_mean = np.nan
    fig.suptitle(
        f"IMS mobility: E={field_Vm:.3e} V/m ({field_source}), L={length:.3e} m, T={temp:.1f} K, p={pressure:.1f} Pa\n"
        f"K median={np.median(mobility):.3e} m²/Vs, K0 median={np.median(K0_cm2):.3e} cm²/Vs | "
        f"length_mode={length_mode}, L_eff_med={leff_med:.3e} m, L_eff_mean={leff_mean:.3e} m, z_arr={arrival_plane_z:.3e} m"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=200)


def _write_csv(arrivals, mobility, K0, species, out_path: Path, drift_lengths_m: np.ndarray, field: FieldStrength):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "species",
            "arrival_time_s",
            "drift_length_m",
            "field_source",
            "E_V_m",
            "EN_Td",
            "voltage_V",
            "K_m2_Vs",
            "K0_m2_Vs",
            "K0_cm2_Vs",
        ])
        for sp, t, Lm, k, k0 in zip(species, arrivals, drift_lengths_m, mobility, K0):
            writer.writerow([
                sp,
                f"{t:.6e}",
                f"{Lm:.6e}",
                field.source,
                f"{field.value_V_m:.6e}",
                "" if field.en_td is None else f"{field.en_td:.6e}",
                "" if field.voltage_V is None else f"{field.voltage_V:.6e}",
                f"{k:.6e}",
                f"{k0:.6e}",
                f"{k0 * 1e4:.6e}",
            ])


if __name__ == "__main__":
    raise SystemExit(main())

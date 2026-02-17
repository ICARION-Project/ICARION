#!/usr/bin/env python3
"""
Estimate ion mobility (K) and reduced mobility (K0) from IMS trajectories.

Workflow:
  - Select ions (species filter, subsampling)
  - Find arrival times at drift tube end (first timestep with z >= z_max - tol and within radial bounds)
  - Compute drift velocity v = L / t_arrival, mobility K = v / E
  - Convert to reduced mobility K0 = K * (P/P0) * (T0/T)
  - Plot arrival-time histogram with mean/median markers and save CSV with per-ion values

Inputs pulled from HDF5 when available:
  - Drift length: /domains/domain_<i>/geometry/length_m
  - Electric field: derived from /domains/domain_<i>/fields/DC/* (voltage) divided by length, if present.
    Fallback: requires --field-V or --field-VM.
  - Temperature, pressure: /domains/domain_<i>/environment/temperature_K, pressure_Pa
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
from typing import Dict, Optional

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
)

P0 = 101325.0  # Pa
T0 = 273.15  # K


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
    p.add_argument("--field-V", dest="field_V", type=float, default=None, help="Axial voltage (V) across drift tube.")
    p.add_argument("--field-VM", dest="field_Vm", type=float, default=None, help="Axial field strength (V/m). Overrides voltage/length if set.")
    p.add_argument("--time-stride", type=int, default=1, help="Stride when scanning positions for arrival (1 = full).")
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

        # Env parameters
        temp_K, pressure_Pa = _read_env(h5, args.domain_index)
        length_m, radius_m = _read_geometry(h5, args.domain_index, traj, ion_indices)
        field_Vm = _read_field_strength(h5, args.domain_index, length_m, args.field_V, args.field_Vm)

        # Arrival detection
        times = traj["time"][:: max(1, args.time_stride)]
        positions = traj["positions"][:: max(1, args.time_stride), ion_indices, :]
        arrivals, mask_valid = find_arrival_times(
            times=times,
            positions=positions,
            length_m=length_m,
            radius_m=radius_m,
            tol_frac=args.tol_frac,
        )
        if not np.any(mask_valid):
            raise RuntimeError("No ions reached the end of the drift region.")

        # Compute mobility
        drift_vel = length_m / arrivals  # m/s
        mobility = drift_vel / field_Vm  # m2/Vs
        K0 = mobility * (pressure_Pa / P0) * (T0 / temp_K)

        species_selected = np.array(species_ids)[ion_indices]
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
    )
    _write_csv(arrivals, mobility, K0, species_valid, args.out_csv)
    print(f"Wrote {args.out}")
    print(f"Wrote {args.out_csv}")
    return 0


def _read_env(h5, domain_idx: int) -> tuple[float, float]:
    base = f"/domains/domain_{domain_idx}/environment"
    temp = float(h5[f"{base}/temperature_K"][()]) if f"{base}/temperature_K" in h5 else np.nan
    pressure = float(h5[f"{base}/pressure_Pa"][()]) if f"{base}/pressure_Pa" in h5 else np.nan
    if np.isnan(temp) or np.isnan(pressure):
        raise KeyError("Missing temperature_K or pressure_Pa in environment.")
    return temp, pressure


def _read_geometry(h5, domain_idx: int, traj_group, ion_indices) -> tuple[float, float]:
    base = f"/domains/domain_{domain_idx}/geometry"
    length = float(h5[base + "/length_m"][()]) if base + "/length_m" in h5 else None
    radius = float(h5[base + "/radius_m"][()]) if base + "/radius_m" in h5 else None
    if length is None:
        # fallback: estimate from trajectory z-span
        z = traj_group["positions"][:, ion_indices, 2]
        length = float(z.max() - z.min())
        print(f"⚠️  length_m missing; using span from trajectories: {length:.3e} m")
    if radius is None and "positions" in traj_group:
        r = np.sqrt(traj_group["positions"][:, ion_indices, 0] ** 2 + traj_group["positions"][:, ion_indices, 1] ** 2)
        radius = float(np.nanmax(r))
        print(f"⚠️  radius_m missing; using max observed radius: {radius:.3e} m")
    return length, radius


def _read_field_strength(
    h5,
    domain_idx: int,
    length_m: float,
    field_V: Optional[float],
    field_Vm_override: Optional[float],
) -> float:
    if field_Vm_override is not None:
        return float(field_Vm_override)
    # try to find DC voltage in HDF5
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
            break
    if voltage is None:
        # Fallback: if EN_Td is stored, compute E from gas density.
        en_paths = [
            f"/domains/domain_{domain_idx}/fields/dc/EN_Td",
            f"/domains/domain_{domain_idx}/fields/DC/EN_Td",
        ]
        en_td = None
        for path in en_paths:
            if path in h5:
                en_td = float(h5[path][()])
                break
        if en_td is not None:
            env_path = f"/domains/domain_{domain_idx}/environment"
            temp = float(h5[f"{env_path}/temperature_K"][()]) if f"{env_path}/temperature_K" in h5 else None
            pressure = float(h5[f"{env_path}/pressure_Pa"][()]) if f"{env_path}/pressure_Pa" in h5 else None
            if temp is None or pressure is None:
                raise KeyError("Missing temperature_K or pressure_Pa for EN_Td conversion.")
            kB = 1.380649e-23  # J/K
            n = pressure / (kB * temp)
            return en_td * 1e-21 * n
        if field_V is None:
            raise KeyError("No DC voltage or EN_Td found in HDF5; provide --field-V or --field-VM.")
        voltage = float(field_V)
    return voltage / length_m


def find_arrival_times(
    times: np.ndarray,
    positions: np.ndarray,
    length_m: float,
    radius_m: float,
    tol_frac: float,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Determine arrival times per ion:
      - condition: z >= z_max - tol and radial <= radius + tol
      - returns (arrival_times, mask_valid)
    """
    tol_z = length_m * tol_frac
    tol_r = radius_m * tol_frac
    z = positions[:, :, 2]  # (T, N)
    r = np.sqrt(positions[:, :, 0] ** 2 + positions[:, :, 1] ** 2)
    z_max = z.max()
    mask_reached = (z >= z_max - tol_z) & (r <= radius_m + tol_r)

    # earliest timestep per ion where condition true
    arrivals = np.full(z.shape[1], np.nan)
    for ion in range(z.shape[1]):
        idx = np.argmax(mask_reached[:, ion]) if np.any(mask_reached[:, ion]) else None
        if idx is not None and mask_reached[idx, ion]:
            arrivals[ion] = times[idx]
    mask_valid = ~np.isnan(arrivals)
    return arrivals[mask_valid], mask_valid


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
):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5))

    # Arrival histogram
    ax_arr = axes[0]
    if len(arrivals) > 1:
        edges = np.linspace(arrivals.min(), arrivals.max(), bins + 1)
    else:
        edges = np.linspace(arrivals.min() * 0.99, arrivals.max() * 1.01 + 1e-12, bins + 1)
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
    fig.suptitle(
        f"IMS mobility: E={field_Vm:.3e} V/m, L={length:.3e} m, T={temp:.1f} K, p={pressure:.1f} Pa\n"
        f"K median={np.median(mobility):.3e} m²/Vs, K0 median={np.median(K0_cm2):.3e} cm²/Vs"
    )
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    fig.savefig(out_path, dpi=200)


def _write_csv(arrivals, mobility, K0, species, out_path: Path):
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["species", "arrival_time_s", "K_m2_Vs", "K0_m2_Vs", "K0_cm2_Vs"])
        for sp, t, k, k0 in zip(species, arrivals, mobility, K0):
            writer.writerow([sp, f"{t:.6e}", f"{k:.6e}", f"{k0:.6e}", f"{k0 * 1e4:.6e}"])


if __name__ == "__main__":
    raise SystemExit(main())

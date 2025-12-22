#!/usr/bin/env python3
"""
Compute mean arrival times and mobilities for an IMS E/N sweep.

- Expects files named like ims_en_<EN>Td.h5 (EN parsed from filename).
- Computes arrival per ion (first time z >= z_max - tol and within radial tol).
- Outputs CSV with per-species stats: mean/median arrival, fraction arrived,
  mobility K (using actual E = (E/N)*N), and reduced mobility K0 using Loschmidt constant.
"""

from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import numpy as np

KB = 1.380649e-23  # J/K
LOSCHMIDT_N0 = 2.686780111e25  # m^-3 at 273.15 K, 1 atm
TD_TO_V_M2 = 1e-21  # 1 Td = 1e-21 V*m^2

FILENAME_RE = re.compile(r"ims_en_(?P<en>[-\deE\.]+)Td\.h5$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Analyze IMS E/N sweep for arrival times and mobility.")
    p.add_argument("--input-dir", type=Path, default=Path("results/ims_en_sweep"),
                   help="Directory with trajectory files.")
    p.add_argument("--pattern", type=str, default="ims_en_*Td.h5",
                   help="Glob pattern for input files.")
    p.add_argument("--out-csv", type=Path, default=Path("analysis/output/ims_en_sweep.csv"),
                   help="CSV output path.")
    p.add_argument("--species", nargs="+", default=None,
                   help="Optional species filter.")
    p.add_argument("--domain-index", type=int, default=0,
                   help="Domain index for geometry/env.")
    p.add_argument("--tol-frac", type=float, default=0.02,
                   help="Tolerance fraction for arrival (of length/radius).")
    p.add_argument("--time-stride", type=int, default=1,
                   help="Stride when scanning times/positions (1 = full).")
    return p.parse_args()


def find_files(input_dir: Path, pattern: str) -> List[Path]:
    return sorted(input_dir.glob(pattern))


def parse_en_td(path: Path) -> float:
    m = FILENAME_RE.search(path.name)
    if not m:
        raise ValueError(f"Cannot parse E/N from filename: {path.name}")
    return float(m.group("en"))


def load_species_ids(traj_group) -> np.ndarray:
    species_ds = traj_group["species_ids"]
    first_row = species_ds[0]
    return np.array([s.decode("utf-8") if isinstance(s, (bytes, bytearray)) else str(s) for s in first_row])


def select_indices(species_ids: np.ndarray, species_filter: set | None) -> np.ndarray:
    if not species_filter:
        return np.arange(len(species_ids))
    mask = np.array([sp in species_filter for sp in species_ids])
    return np.where(mask)[0]


def read_env(h5, domain_idx: int) -> Tuple[float, float]:
    base = f"/domains/domain_{domain_idx}/environment"
    t = float(h5[f"{base}/temperature_K"][()]) if f"{base}/temperature_K" in h5 else np.nan
    p = float(h5[f"{base}/pressure_Pa"][()]) if f"{base}/pressure_Pa" in h5 else np.nan
    if np.isnan(t) or np.isnan(p):
        raise KeyError("Missing temperature_K or pressure_Pa in environment.")
    return t, p


def read_length_radius(h5, domain_idx: int) -> Tuple[float, float]:
    base = f"/domains/domain_{domain_idx}/geometry"
    length = float(h5[f"{base}/length_m"][()]) if f"{base}/length_m" in h5 else None
    radius = float(h5[f"{base}/radius_m"][()]) if f"{base}/radius_m" in h5 else None
    if length is None or radius is None:
        raise KeyError("Missing geometry length_m or radius_m.")
    return length, radius


def find_arrivals(times: np.ndarray, positions: np.ndarray,
                  length_m: float, radius_m: float, tol_frac: float) -> np.ndarray:
    tol_z = length_m * tol_frac
    tol_r = radius_m * tol_frac
    z = positions[:, :, 2]
    r = np.sqrt(positions[:, :, 0] ** 2 + positions[:, :, 1] ** 2)
    z_max = z.max()
    mask = (z >= z_max - tol_z) & (r <= radius_m + tol_r)
    arrivals = np.full(z.shape[1], np.nan)
    for i in range(z.shape[1]):
        idx = np.argmax(mask[:, i]) if np.any(mask[:, i]) else None
        if idx is not None and mask[idx, i]:
            arrivals[i] = times[idx]
    return arrivals


def per_species_stats(arrivals: np.ndarray, species: np.ndarray,
                      en_td: float, length_m: float,
                      temp_K: float, pressure_Pa: float) -> List[Dict]:
    stats = []
    # Number density
    N = pressure_Pa / (KB * temp_K)
    EN_SI = en_td * TD_TO_V_M2
    for sp in sorted(set(species)):
        mask = species == sp
        arr_sp = arrivals[mask]
        valid = arr_sp[~np.isnan(arr_sp)]
        if valid.size == 0:
            continue
        mean_t = float(np.mean(valid))
        median_t = float(np.median(valid))
        frac = float(valid.size) / float(arr_sp.size)
        drift_v = length_m / mean_t
        K = drift_v / (EN_SI * N)
        K0 = drift_v / (EN_SI * LOSCHMIDT_N0)
        stats.append({
            "species": sp,
            "en_td": en_td,
            "mean_arrival_s": mean_t,
            "median_arrival_s": median_t,
            "fraction_arrived": frac,
            "mobility_m2_Vs": K,
            "reduced_mobility_m2_Vs": K0,
        })
    return stats


def process_file(path: Path, species_filter: set | None,
                 domain_idx: int, tol_frac: float, stride: int) -> List[Dict]:
    en_td = parse_en_td(path)
    with h5py.File(path, "r") as h5:
        if "trajectory" not in h5:
            raise KeyError("Missing /trajectory group")
        traj = h5["trajectory"]
        species_ids = load_species_ids(traj)
        idx = select_indices(species_ids, species_filter)
        if idx.size == 0:
            return []
        times = traj["time"][:: max(1, stride)]
        positions = traj["positions"][:: max(1, stride), idx, :]
        length_m, radius_m = read_length_radius(h5, domain_idx)
        temp_K, pressure_Pa = read_env(h5, domain_idx)
        arrivals = find_arrivals(times, positions, length_m, radius_m, tol_frac)
        species_sel = species_ids[idx]
        return per_species_stats(arrivals, species_sel, en_td, length_m, temp_K, pressure_Pa)


def write_csv(rows: List[Dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        print("No data to write.")
        return
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {path}")


def main() -> int:
    args = parse_args()
    files = find_files(args.input_dir, args.pattern)
    if not files:
        raise FileNotFoundError(f"No files matching {args.pattern} in {args.input_dir}")
    species_filter = set(args.species) if args.species else None
    all_rows: List[Dict] = []
    for f in files:
        try:
            rows = process_file(
                f, species_filter,
                domain_idx=args.domain_index,
                tol_frac=args.tol_frac,
                stride=max(1, args.time_stride),
            )
            all_rows.extend(rows)
        except Exception as exc:
            print(f"⚠️  Skip {f}: {exc}")
    # sort by species then EN
    all_rows.sort(key=lambda r: (r["species"], r["en_td"]))
    write_csv(all_rows, args.out_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

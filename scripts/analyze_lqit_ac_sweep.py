#!/usr/bin/env python3
"""
Aggregate LQIT AC sweep results into per-species heatmaps:
- For each (pressure, frequency) point: compute max over time of mean radial position per species.
- Save one heatmap per species (pressure vs frequency).

Example:
  python tmp/analyze_lqit_ac_sweep.py \
      --input-dir results/lqit_ac_sweep \
      --out-dir tmp/lqit_ac_sweep/analysis
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List, Tuple

import h5py
import matplotlib.pyplot as plt
import numpy as np


FILENAME_RE = re.compile(r"p(?P<p>[-\deE\.]+)Pa_f(?P<f>\d+)Hz\.h5$")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="LQIT AC sweep heatmaps (max mean radial per species).")
    p.add_argument("--input-dir", type=Path, default=Path("results/lqit_ac_sweep"),
                   help="Directory with trajectory HDF5 files.")
    p.add_argument("--out-dir", type=Path, default=Path("tmp/lqit_ac_sweep/analysis"),
                   help="Directory to write heatmaps.")
    p.add_argument("--pattern", type=str, default="lqit_ac_p*Pa_f*Hz.h5",
                   help="Glob pattern for files.")
    p.add_argument("--species", nargs="+", default=None,
                   help="Optional list of species to include (default: all present).")
    p.add_argument("--chunk", type=int, default=128,
                   help="Timesteps per chunk when reading positions.")
    return p.parse_args()


def find_files(input_dir: Path, pattern: str) -> List[Path]:
    return sorted(input_dir.glob(pattern))


def parse_pf_from_name(path: Path) -> Tuple[float, float]:
    m = FILENAME_RE.search(path.name)
    if not m:
        raise ValueError(f"Cannot parse pressure/frequency from filename: {path.name}")
    return float(m.group("p")), float(m.group("f"))


def load_species_ids(traj_group) -> np.ndarray:
    species_ds = traj_group["species_ids"]
    first_row = species_ds[0]
    return np.array([s.decode("utf-8") if isinstance(s, (bytes, bytearray)) else str(s) for s in first_row])


def max_mean_radial_for_file(path: Path, species_filter: set | None, chunk_size: int) -> Tuple[float, float, Dict[str, float]]:
    p, f = parse_pf_from_name(path)
    with h5py.File(path, "r") as h5:
        if "trajectory" not in h5:
            raise KeyError(f"{path} missing /trajectory group")
        traj = h5["trajectory"]
        time = traj["time"]
        positions = traj["positions"]  # (T, N, 3)
        species_ids = load_species_ids(traj)

        if species_filter:
            selected = [s for s in species_ids if s in species_filter]
            if not selected:
                return p, f, {}
            unique_species = sorted(set(selected))
        else:
            unique_species = sorted(set(species_ids))

        masks = {sp: species_ids == sp for sp in unique_species}
        max_means = {sp: -np.inf for sp in unique_species}

        n_steps = time.shape[0]
        for start in range(0, n_steps, chunk_size):
            end = min(n_steps, start + chunk_size)
            pos_chunk = positions[start:end, :, :]  # (chunk, N, 3)
            r_chunk = np.sqrt(pos_chunk[:, :, 0] ** 2 + pos_chunk[:, :, 1] ** 2)
            for sp, mask in masks.items():
                subset = r_chunk[:, mask]
                mean_r = subset.mean(axis=1)
                max_means[sp] = max(max_means[sp], float(mean_r.max()))

        return p, f, max_means


def aggregate(files: List[Path], species_filter: set | None, chunk_size: int):
    pressures = set()
    freqs = set()
    data: Dict[str, Dict[Tuple[float, float], float]] = {}

    for path in files:
        try:
            p, f, max_means = max_mean_radial_for_file(path, species_filter, chunk_size)
        except Exception as exc:
            print(f"⚠️  Skip {path}: {exc}")
            continue
        if not max_means:
            continue
        pressures.add(p)
        freqs.add(f)
        for sp, val in max_means.items():
            data.setdefault(sp, {})[(p, f)] = val

    pressures = sorted(pressures)
    freqs = sorted(freqs)
    return pressures, freqs, data


def plot_heatmaps(pressures: List[float], freqs: List[float], data, out_dir: Path):
    out_dir.mkdir(parents=True, exist_ok=True)
    P, F = np.meshgrid(pressures, freqs, indexing="ij")  # shape (P,F)

    for sp, vals in data.items():
        grid = np.full_like(P, fill_value=np.nan, dtype=float)
        for (p, f), v in vals.items():
            i = pressures.index(p)
            j = freqs.index(f)
            grid[i, j] = v

        plt.figure(figsize=(8, 5))
        im = plt.pcolormesh(F, P, grid, shading="auto", cmap="viridis")
        plt.colorbar(im, label="Max mean radial position [m]")
        plt.xlabel("AC frequency [Hz]")
        plt.ylabel("Pressure [Pa]")
        plt.title(f"LQIT AC sweep – {sp}")
        plt.tight_layout()

        out_path = out_dir / f"heatmap_{sp}.png"
        plt.savefig(out_path, dpi=200)
        plt.close()
        print(f"Wrote {out_path}")


def main():
    args = parse_args()
    files = find_files(args.input_dir, args.pattern)
    if not files:
        raise FileNotFoundError(f"No files matching {args.pattern} in {args.input_dir}")

    species_filter = set(args.species) if args.species else None
    pressures, freqs, data = aggregate(files, species_filter, args.chunk)
    if not data:
        raise RuntimeError("No data aggregated (species filter too strict or files invalid).")

    plot_heatmaps(pressures, freqs, data, args.out_dir)


if __name__ == "__main__":
    main()

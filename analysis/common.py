#!/usr/bin/env python3
"""
Shared helpers for analysis scripts.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable, Sequence

import h5py
import numpy as np


def open_trajectory(path: Path) -> h5py.File:
    """Open trajectory file and verify required groups."""
    handle = h5py.File(path, "r")
    if "trajectory" not in handle:
        handle.close()
        raise KeyError(f"{path} missing '/trajectory' group")
    return handle


def load_species_ids(traj_group) -> np.ndarray:
    """Decode species IDs from the first timestep."""
    species_ds = traj_group["species_ids"]
    first_row = species_ds[0]
    return np.array([_decode(v) for v in first_row])


def _decode(val) -> str:
    if isinstance(val, (bytes, bytearray)):
        return val.decode("utf-8")
    return str(val)


def select_ion_indices(
    species_ids: Sequence[str],
    species_filter: Iterable[str] | None,
    max_ions: int | None,
    max_per_species: int | None,
    rng_seed: int | None,
) -> np.ndarray:
    """
    Pick ion indices with optional filtering and per-species caps.
    Strategy:
      1) Filter to requested species (if any)
      2) Cap per-species counts (if set)
      3) Cap total ions (if set) via random sampling
    """
    species_ids = np.asarray(species_ids)
    if species_filter:
        filt = set(species_filter)
        mask = np.array([s in filt for s in species_ids])
    else:
        mask = np.ones(len(species_ids), dtype=bool)

    candidates = np.nonzero(mask)[0]
    rng = np.random.default_rng(rng_seed)

    if max_per_species is not None:
        selected = []
        for sp in sorted(set(species_ids[candidates])):
            sp_indices = np.nonzero((species_ids == sp) & mask)[0]
            rng.shuffle(sp_indices)
            selected.extend(sp_indices[:max_per_species])
        candidates = np.array(sorted(selected))

    if max_ions is not None and len(candidates) > max_ions:
        candidates = rng.choice(candidates, size=max_ions, replace=False)
        candidates.sort()

    return candidates


def load_positions_subset(
    traj_group,
    ion_indices: Sequence[int],
    time_stride: int,
    max_frames: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Load time and positions for a subset of ions with optional time stride/frame cap.
    Returns (time[T], positions[T, n_ions, 3]).
    """
    positions_ds = traj_group["positions"]
    time_ds = traj_group["time"]

    stride = max(1, int(time_stride))
    if max_frames is not None:
        time = time_ds[::stride][:max_frames]
        positions = positions_ds[0 : time.shape[0] * stride : stride, ion_indices, :]
    else:
        time = time_ds[::stride]
        positions = positions_ds[::stride, ion_indices, :]

    return time, positions


def species_color_map(species: Sequence[str]) -> dict[str, str]:
    """Assign distinct-ish matplotlib colors to species list."""
    import matplotlib.pyplot as plt  # Lazy import to keep dependencies local

    unique = sorted(set(species))
    cmap = plt.cm.get_cmap("tab20", len(unique))
    return {sp: cmap(i) for i, sp in enumerate(unique)}

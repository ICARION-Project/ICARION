#!/usr/bin/env python3
"""Trajectory HDF5 access and selection helpers for analysis scripts."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Sequence

import h5py
import numpy as np


@dataclass(frozen=True)
class TrajectorySelection:
    """In-memory subset loaded from a trajectory file for analysis scripts."""

    time: np.ndarray
    positions: np.ndarray
    species_ids_all: np.ndarray
    ion_indices: np.ndarray
    species_ids_selected: np.ndarray


@dataclass(frozen=True)
class DomainBounds:
    """Radial and axial bounds used by elimination-style analyses."""

    r_in: float | None
    r_out: float | None
    z_min: float | None
    z_max: float | None


def normalize_species_filter(species: Iterable[str] | None) -> set[str] | None:
    """Normalize optional species arguments to a trimmed set, or None."""
    if species is None:
        return None
    normalized = {str(sp).strip() for sp in species if str(sp).strip()}
    return normalized or None


def decode_hdf5_string(val) -> str:
    """Decode HDF5 bytes/scalars to Python strings."""
    if isinstance(val, (bytes, bytearray)):
        return val.decode("utf-8")
    return str(val)


def open_trajectory(path: Path) -> h5py.File:
    """Open trajectory file and verify required groups."""
    handle = h5py.File(path, "r")
    if "trajectory" not in handle:
        handle.close()
        raise KeyError(f"{path} missing '/trajectory' group")
    return handle


def read_species_ids(traj_group, frame_index: int = 0) -> np.ndarray:
    """Decode species IDs from a requested timestep or a 1D dataset."""
    species_ds = traj_group["species_ids"]
    row = species_ds[frame_index] if species_ds.ndim >= 2 else species_ds[()]
    return np.array([decode_hdf5_string(v) for v in row])


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
      1. Filter to requested species, if any.
      2. Cap per-species counts, if set.
      3. Cap total ions, if set, via random sampling.
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


def select_ions_from_trajectory(
    traj_group,
    species_filter: Iterable[str] | None,
    max_ions: int | None,
    max_per_species: int | None,
    rng_seed: int | None,
    species_frame_index: int = 0,
) -> tuple[np.ndarray, np.ndarray]:
    """Decode species IDs and select ion indices with standard filtering/caps."""
    species_ids = read_species_ids(traj_group, frame_index=species_frame_index)
    ion_indices = select_ion_indices(
        species_ids,
        species_filter=species_filter,
        max_ions=max_ions,
        max_per_species=max_per_species,
        rng_seed=rng_seed,
    )
    return species_ids, ion_indices


def ensure_nonempty_selection(ion_indices: Sequence[int], reason: str = "broaden species filter or caps") -> None:
    """Raise a consistent error when no ions are selected."""
    if len(ion_indices) == 0:
        raise RuntimeError(f"No ions selected; {reason}.")


def species_for_indices(species_ids: Sequence[str], ion_indices: Sequence[int]) -> np.ndarray:
    """Return species IDs corresponding to selected ion indices."""
    return np.asarray(species_ids)[np.asarray(ion_indices, dtype=int)]


def read_death_times_for_indices(traj_path: Path, ion_indices: Sequence[int]) -> np.ndarray | None:
    """Load ions/death_time_s for selected indices; returns None if dataset is missing."""
    with open_trajectory(traj_path) as h5:
        if "ions/death_time_s" not in h5:
            return None
        return np.asarray(h5["ions/death_time_s"][:], dtype=float)[np.asarray(ion_indices, dtype=int)]


def active_ion_mask(time: np.ndarray, death_times: np.ndarray) -> np.ndarray:
    """Build a [T, N] mask of ions active at each sampled timepoint."""
    time = np.asarray(time, dtype=float)
    death_times = np.asarray(death_times, dtype=float)
    mask = death_times[None, :] > time[:, None]
    mask[:, death_times <= 0.0] = True
    return mask


def read_positions_subset(
    traj_group,
    ion_indices: Sequence[int],
    time_stride: int,
    max_frames: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Load time and positions for a subset of ions with optional time stride/frame cap.

    Returns ``(time[T], positions[T, n_ions, 3])``.
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


def read_domain_geometry(
    h5,
    domain_idx: int,
    traj_group,
    ion_indices: Sequence[int],
    warning_prefix: str = "Warning:",
) -> tuple[float, float, float]:
    """Load drift geometry from HDF5 with trajectory-based fallbacks."""
    base = f"/domains/domain_{domain_idx}/geometry"
    length = float(h5[base + "/length_m"][()]) if base + "/length_m" in h5 else None
    radius = float(h5[base + "/radius_m"][()]) if base + "/radius_m" in h5 else None
    z_origin = None
    if base + "/origin_m" in h5:
        origin = np.asarray(h5[base + "/origin_m"][()], dtype=float)
        if origin.size >= 3:
            z_origin = float(origin[2])
    if length is None:
        z = traj_group["positions"][:, ion_indices, 2]
        length = float(z.max() - z.min())
        print(f"{warning_prefix} length_m missing; using trajectory span {length:.3e} m")
    if z_origin is None:
        z0 = np.asarray(traj_group["positions"][0, ion_indices, 2], dtype=float)
        z_origin = float(np.median(z0))
        print(f"{warning_prefix} origin_m missing; using median initial z {z_origin:.3e} m")
    if radius is None and "positions" in traj_group:
        r = np.sqrt(traj_group["positions"][:, ion_indices, 0] ** 2 + traj_group["positions"][:, ion_indices, 1] ** 2)
        radius = float(np.nanmax(r))
        print(f"{warning_prefix} radius_m missing; using max observed radius {radius:.3e} m")
    return length, radius, z_origin


def read_domain_environment(
    h5,
    domain_idx: int,
    require_pressure: bool = True,
    allow_missing: bool = False,
) -> tuple[float | None, float | None]:
    """Read domain temperature/pressure from HDF5 environment metadata."""
    base = f"/domains/domain_{domain_idx}/environment"
    temp = float(h5[f"{base}/temperature_K"][()]) if f"{base}/temperature_K" in h5 else None
    pressure = float(h5[f"{base}/pressure_Pa"][()]) if f"{base}/pressure_Pa" in h5 else None
    if allow_missing:
        return temp, pressure
    if temp is None or (require_pressure and pressure is None):
        if require_pressure:
            raise KeyError("Missing temperature_K or pressure_Pa in environment.")
        raise KeyError("Missing temperature_K in environment.")
    return temp, pressure


def read_species_mass_map(h5) -> dict[str, float]:
    """Read species-name to mass mapping from metadata."""
    names_path = "/metadata/species/names"
    masses_path = "/metadata/species/mass_kg"
    if names_path not in h5 or masses_path not in h5:
        raise KeyError("Missing /metadata/species datasets (names, mass_kg).")
    names_raw = h5[names_path][:]
    masses = np.asarray(h5[masses_path][:], dtype=float)
    names = [decode_hdf5_string(val) for val in names_raw]
    return dict(zip(names, masses))


def species_masses_for_selection(species_ids_selected: Sequence[str], mass_map: dict[str, float]) -> np.ndarray:
    """Map selected species IDs to per-ion masses and fail on missing entries."""
    masses = np.asarray([mass_map.get(sp, np.nan) for sp in species_ids_selected], dtype=float)
    if np.any(np.isnan(masses)):
        missing = [sp for sp, mass in zip(species_ids_selected, masses) if np.isnan(mass)]
        raise KeyError(f"Missing mass for species: {sorted(set(missing))}")
    return masses


def find_arrival_times(
    times: np.ndarray,
    positions: np.ndarray,
    z_origin_m: float,
    length_m: float,
    radius_m: float,
    tol_frac: float,
    arrival_plane_z_m: float | None = None,
) -> tuple[np.ndarray, np.ndarray, float]:
    """Return first arrival times for ions reaching the end plane within radial tolerance."""
    tol_z = length_m * tol_frac
    tol_r = radius_m * tol_frac
    z = positions[:, :, 2]
    r = np.sqrt(positions[:, :, 0] ** 2 + positions[:, :, 1] ** 2)
    z_end = z_origin_m + length_m if arrival_plane_z_m is None else arrival_plane_z_m
    z_target = z_end - tol_z if arrival_plane_z_m is None else z_end
    mask_reached = (z >= z_target) & (r <= radius_m + tol_r)

    arrivals = np.full(z.shape[1], np.nan, dtype=float)
    for ion in range(z.shape[1]):
        if np.any(mask_reached[:, ion]):
            arrivals[ion] = float(times[int(np.argmax(mask_reached[:, ion]))])
    mask_valid = ~np.isnan(arrivals)
    return arrivals[mask_valid], mask_valid, z_target


def read_domain_bounds(
    h5,
    domain_idx: int,
    positions_last: np.ndarray,
) -> DomainBounds:
    """Read radial and axial bounds for elimination-style analyses."""
    geom_group = h5.get(f"domains/domain_{domain_idx}/geometry")

    def _scalar(name: str) -> float | None:
        if geom_group is None or name not in geom_group:
            return None
        arr = np.asarray(geom_group[name], dtype=float)
        if arr.size == 0:
            return None
        return float(arr.flat[0])

    origin_z = 0.0
    if geom_group is not None and "origin_m" in geom_group:
        origin = np.asarray(geom_group["origin_m"], dtype=float)
        if origin.size >= 3:
            origin_z = float(origin[2])

    length = _scalar("length_m")
    r_out = _scalar("radius_out_m") or _scalar("radius_m")
    r_in = _scalar("radius_in_m")

    if length is not None:
        z_min = origin_z
        z_max = origin_z + length
    else:
        z_min = float(positions_last[:, 2].min())
        z_max = float(positions_last[:, 2].max())

    return DomainBounds(r_in=r_in, r_out=r_out, z_min=z_min, z_max=z_max)


def load_trajectory_selection(
    traj_path: Path,
    species_filter: Iterable[str] | None,
    max_ions: int | None,
    max_per_species: int | None,
    rng_seed: int | None,
    time_stride: int,
    max_frames: int | None,
) -> TrajectorySelection:
    """Load a filtered trajectory subset into memory using shared defaults/logic."""
    if not traj_path.exists():
        raise FileNotFoundError(f"Trajectory file not found: {traj_path}")

    with open_trajectory(traj_path) as h5:
        traj = h5["trajectory"]
        species_ids, ion_indices = select_ions_from_trajectory(
            traj,
            species_filter=species_filter,
            max_ions=max_ions,
            max_per_species=max_per_species,
            rng_seed=rng_seed,
        )
        ensure_nonempty_selection(ion_indices)

        time, positions = read_positions_subset(
            traj,
            ion_indices=ion_indices,
            time_stride=time_stride,
            max_frames=max_frames,
        )

    species_selected = species_for_indices(species_ids, ion_indices)
    return TrajectorySelection(
        time=np.asarray(time),
        positions=np.asarray(positions),
        species_ids_all=np.asarray(species_ids),
        ion_indices=np.asarray(ion_indices),
        species_ids_selected=np.asarray(species_selected),
    )


def read_trajectory_selection(
    traj_path: Path,
    species_filter: Iterable[str] | None,
    max_ions: int | None,
    max_per_species: int | None,
    rng_seed: int | None,
    time_stride: int,
    max_frames: int | None,
) -> TrajectorySelection:
    """Canonical name for loading filtered trajectory selections."""
    return load_trajectory_selection(
        traj_path=traj_path,
        species_filter=species_filter,
        max_ions=max_ions,
        max_per_species=max_per_species,
        rng_seed=rng_seed,
        time_stride=time_stride,
        max_frames=max_frames,
    )


def load_species_ids(traj_group, frame_index: int = 0) -> np.ndarray:
    """Backward-compatible alias for read_species_ids."""
    return read_species_ids(traj_group, frame_index=frame_index)


def load_death_times_for_indices(traj_path: Path, ion_indices: Sequence[int]) -> np.ndarray | None:
    """Backward-compatible alias for read_death_times_for_indices."""
    return read_death_times_for_indices(traj_path, ion_indices)


def load_positions_subset(
    traj_group,
    ion_indices: Sequence[int],
    time_stride: int,
    max_frames: int | None,
) -> tuple[np.ndarray, np.ndarray]:
    """Backward-compatible alias for read_positions_subset."""
    return read_positions_subset(
        traj_group,
        ion_indices=ion_indices,
        time_stride=time_stride,
        max_frames=max_frames,
    )

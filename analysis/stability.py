#!/usr/bin/env python3
"""Trap stability and trajectory envelope helpers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class StabilitySummary:
    n_ions: int
    n_transmitted: int
    transmission_fraction: float
    max_radial_m: float
    final_radial_mean_m: float
    final_axial_mean_m: float


def radial_distance(positions: np.ndarray) -> np.ndarray:
    """Return radial distance sqrt(x^2+y^2) for trajectory positions."""
    pos = np.asarray(positions, dtype=float)
    return np.sqrt(pos[..., 0] ** 2 + pos[..., 1] ** 2)


def transmission_mask(
    positions: np.ndarray,
    radius_m: float | None = None,
    z_min_m: float | None = None,
    z_max_m: float | None = None,
) -> np.ndarray:
    """Classify ions as transmitted from final position and optional radial/axial bounds."""
    final = np.asarray(positions, dtype=float)[-1]
    mask = np.ones(final.shape[0], dtype=bool)
    if radius_m is not None:
        mask &= np.sqrt(final[:, 0] ** 2 + final[:, 1] ** 2) <= radius_m
    if z_min_m is not None:
        mask &= final[:, 2] >= z_min_m
    if z_max_m is not None:
        mask &= final[:, 2] <= z_max_m
    return mask


def summarize_stability(
    positions: np.ndarray,
    radius_m: float | None = None,
    z_min_m: float | None = None,
    z_max_m: float | None = None,
) -> StabilitySummary:
    """Summarize survival/transmission and radial envelope for a trajectory block."""
    pos = np.asarray(positions, dtype=float)
    if pos.ndim != 3 or pos.shape[-1] != 3:
        raise ValueError("positions must have shape (time, ion, xyz)")
    r = radial_distance(pos)
    transmitted = transmission_mask(pos, radius_m=radius_m, z_min_m=z_min_m, z_max_m=z_max_m)
    n_ions = int(pos.shape[1])
    n_transmitted = int(np.count_nonzero(transmitted))
    return StabilitySummary(
        n_ions=n_ions,
        n_transmitted=n_transmitted,
        transmission_fraction=float(n_transmitted / n_ions) if n_ions else 0.0,
        max_radial_m=float(np.nanmax(r)) if r.size else float("nan"),
        final_radial_mean_m=float(np.nanmean(r[-1])) if r.size else float("nan"),
        final_axial_mean_m=float(np.nanmean(pos[-1, :, 2])) if n_ions else float("nan"),
    )


def survival_fraction_series(positions: np.ndarray, radius_m: float | None = None) -> np.ndarray:
    """Return fraction of ions still inside a radial aperture at each frame."""
    pos = np.asarray(positions, dtype=float)
    if radius_m is None:
        finite = np.all(np.isfinite(pos), axis=-1)
        return np.mean(finite, axis=1)
    return np.mean(radial_distance(pos) <= radius_m, axis=1)

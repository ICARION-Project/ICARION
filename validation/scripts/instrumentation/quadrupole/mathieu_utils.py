#!/usr/bin/env python3
"""Shared Mathieu stability helpers for quadrupole validation scripts."""

from __future__ import annotations

from functools import lru_cache
from typing import Optional, Tuple

from scipy.special import mathieu_a, mathieu_b

# The first stability zone of the Mathieu equation terminates at q_std ≈ 0.908.
_FIRST_ZONE_Q_MAX = 0.908
# Numerical tolerance for boundary checks (dimensionless parameters)
_EPS = 1e-8


def _normalized_params(a: float, q: float) -> Tuple[float, float]:
    """Return the Mathieu parameters in the form used by scipy.special."""

    return float(a), float(q)


@lru_cache(maxsize=1024)
def _first_zone_bounds(q_std: float) -> Optional[Tuple[float, float]]:
    """Return the (lower, upper) bounds of the first stability zone."""

    if q_std < 0.0 or q_std > _FIRST_ZONE_Q_MAX:
        return None

    a0 = float(mathieu_a(0, q_std))
    b1 = float(mathieu_b(1, q_std))

    # For a quadrupole both radial axes must be simultaneously stable.
    # The y-axis sees ``-a`` in the Mathieu equation, so the valid region is
    # the intersection of the stable zones for +a and -a.
    lower = max(a0, -b1)
    upper = min(b1, -a0)

    return lower, upper


def in_first_stability_region(a: float, q: float, tol: float = _EPS) -> bool:
    """Return True if the (a,q) point lies inside the first stability region."""

    a_std, q_std = _normalized_params(a, q)
    bounds = _first_zone_bounds(q_std)
    if bounds is None:
        return False

    lower, upper = bounds
    return (lower - tol) <= a_std <= (upper + tol)

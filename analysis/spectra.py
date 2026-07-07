#!/usr/bin/env python3
"""Spectrum and frequency helpers for ICARION analysis scripts."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class Peak:
    x: float
    intensity: float


@dataclass(frozen=True)
class TimeGridStats:
    n_samples: int
    dt_min_s: float
    dt_median_s: float
    dt_max_s: float
    dt_jitter_frac: float
    is_uniform: bool


@dataclass(frozen=True)
class ResampledSignal:
    times_s: np.ndarray
    signal: np.ndarray
    dt_s: float
    original_stats: TimeGridStats


def histogram_spectrum(values: np.ndarray, bins: int, density: bool = False) -> tuple[np.ndarray, np.ndarray]:
    """Return histogram bin centers and counts/intensity for scalar spectrum values."""
    vals = np.asarray(values, dtype=float)
    vals = vals[np.isfinite(vals)]
    if vals.size == 0:
        raise ValueError("Cannot build spectrum from empty values")
    counts, edges = np.histogram(vals, bins=bins, density=density)
    centers = 0.5 * (edges[:-1] + edges[1:])
    return centers, counts


def time_grid_stats(times_s: np.ndarray, max_jitter_frac: float = 1e-6) -> TimeGridStats:
    """Return basic sampling statistics and uniformity classification."""
    times = np.asarray(times_s, dtype=float)
    if times.ndim != 1:
        raise ValueError("times_s must be a 1D array")
    if times.size < 2:
        raise ValueError("Need at least two samples for time-grid statistics")
    dt = np.diff(times)
    if np.any(dt <= 0.0):
        raise ValueError("times_s must be strictly increasing")
    dt_min = float(np.min(dt))
    dt_med = float(np.median(dt))
    dt_max = float(np.max(dt))
    jitter = float((dt_max - dt_min) / dt_med) if dt_med > 0.0 else float("inf")
    return TimeGridStats(
        n_samples=int(times.size),
        dt_min_s=dt_min,
        dt_median_s=dt_med,
        dt_max_s=dt_max,
        dt_jitter_frac=jitter,
        is_uniform=bool(jitter <= max_jitter_frac),
    )


def resample_to_uniform_grid(times_s: np.ndarray, signal: np.ndarray, dt_s: float | None = None) -> ResampledSignal:
    """Interpolate a scalar signal onto a uniform time grid."""
    times = np.asarray(times_s, dtype=float)
    sig = np.asarray(signal, dtype=float)
    if times.ndim != 1 or sig.ndim != 1 or times.size != sig.size:
        raise ValueError("times_s and signal must be 1D arrays with the same length")
    if times.size < 3:
        raise ValueError("Need at least three samples for resampling")
    stats = time_grid_stats(times)
    if dt_s is None:
        dt_s = stats.dt_median_s
    dt = float(dt_s)
    if dt <= 0.0:
        raise ValueError("resample dt must be positive")
    n = int(np.floor((times[-1] - times[0]) / dt)) + 1
    if n < 3:
        raise ValueError("Resampled grid would have fewer than three samples")
    uniform_times = times[0] + np.arange(n, dtype=float) * dt
    uniform_times[-1] = min(uniform_times[-1], times[-1])
    finite = np.isfinite(sig)
    if np.count_nonzero(finite) < 2:
        raise ValueError("Need at least two finite signal samples for resampling")
    uniform_signal = np.interp(uniform_times, times[finite], sig[finite])
    return ResampledSignal(times_s=uniform_times, signal=uniform_signal, dt_s=dt, original_stats=stats)


def fft_frequency_spectrum(
    times_s: np.ndarray,
    signal: np.ndarray,
    min_frequency_hz: float = 0.0,
    max_jitter_frac: float | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute a one-sided FFT amplitude spectrum for a uniformly sampled signal."""
    times = np.asarray(times_s, dtype=float)
    sig = np.asarray(signal, dtype=float)
    if times.ndim != 1 or sig.ndim != 1 or times.size != sig.size:
        raise ValueError("times_s and signal must be 1D arrays with the same length")
    if times.size < 3:
        raise ValueError("Need at least three samples for FFT spectrum")
    if max_jitter_frac is not None:
        stats = time_grid_stats(times, max_jitter_frac=max_jitter_frac)
        if not stats.is_uniform:
            raise ValueError(
                "Nonuniform time grid for FFT: "
                f"dt_min={stats.dt_min_s:.6g}s, dt_median={stats.dt_median_s:.6g}s, "
                f"dt_max={stats.dt_max_s:.6g}s, jitter={stats.dt_jitter_frac:.6g}; "
                "use --resample uniform or increase --max-dt-jitter-frac"
            )
        dt = stats.dt_median_s
    else:
        dt = float(np.median(np.diff(times)))
        if dt <= 0.0:
            raise ValueError("times_s must be strictly increasing")
    windowed = (sig - np.nanmean(sig)) * np.hanning(sig.size)
    freq = np.fft.rfftfreq(sig.size, d=dt)
    amp = np.abs(np.fft.rfft(windowed))
    mask = freq >= min_frequency_hz
    return freq[mask], amp[mask]


def top_peaks(x: np.ndarray, y: np.ndarray, n: int = 5, min_separation_bins: int = 1) -> list[Peak]:
    """Return the strongest local maxima with a simple bin-separation guard."""
    xx = np.asarray(x, dtype=float)
    yy = np.asarray(y, dtype=float)
    if xx.size != yy.size:
        raise ValueError("x and y must have the same length")
    if xx.size == 0:
        return []

    if xx.size < 3:
        idx = np.argsort(yy)[::-1][:n]
        return [Peak(float(xx[i]), float(yy[i])) for i in idx]

    candidates = np.nonzero((yy[1:-1] >= yy[:-2]) & (yy[1:-1] >= yy[2:]))[0] + 1
    if candidates.size == 0:
        candidates = np.arange(yy.size)
    order = candidates[np.argsort(yy[candidates])[::-1]]

    selected: list[int] = []
    sep = max(1, int(min_separation_bins))
    for idx in order:
        if all(abs(int(idx) - int(existing)) >= sep for existing in selected):
            selected.append(int(idx))
        if len(selected) >= n:
            break
    return [Peak(float(xx[i]), float(yy[i])) for i in selected]

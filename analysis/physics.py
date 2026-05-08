#!/usr/bin/env python3
"""Physics and numerical helpers shared by analysis scripts."""

from __future__ import annotations

import numpy as np

BOLTZMANN_CONSTANT = 1.380649e-23  # J/K
STANDARD_PRESSURE_PA = 101325.0
STANDARD_TEMPERATURE_K = 273.15
TOWNSEND_TO_V_M2 = 1e-21
ELEMENTARY_CHARGE_C = 1.602176634e-19
ATOMIC_MASS_KG = 1.66053906660e-27


def temperature_series_from_velocities(velocities: np.ndarray, masses: np.ndarray) -> np.ndarray:
    """Compute instantaneous temperature series from velocity snapshots."""
    velocities = np.asarray(velocities, dtype=float)
    masses = np.asarray(masses, dtype=float)
    v2 = np.sum(velocities**2, axis=-1)
    energy = v2 * masses
    mean_energy = energy.mean(axis=1)
    return mean_energy / (3.0 * BOLTZMANN_CONSTANT)


def temperature_from_velocities(velocities: np.ndarray, masses: np.ndarray) -> float:
    """Compute aggregate temperature from velocity snapshots."""
    velocities = np.asarray(velocities, dtype=float)
    masses = np.asarray(masses, dtype=float)
    v2 = np.sum(velocities**2, axis=-1)
    energy = v2 * masses
    return float(np.mean(energy) / (3.0 * BOLTZMANN_CONSTANT))


def mass_to_charge_da_per_e(mass_kg: np.ndarray, charge_C: np.ndarray) -> np.ndarray:
    """Convert SI mass/charge values to nominal m/z in Da/e."""
    mass = np.asarray(mass_kg, dtype=float)
    charge = np.asarray(charge_C, dtype=float)
    z = np.abs(charge) / ELEMENTARY_CHARGE_C
    out = np.full(np.broadcast_shapes(mass.shape, z.shape), np.nan, dtype=float)
    return np.divide(mass / ATOMIC_MASS_KG, z, out=out, where=z != 0.0)


def kinetic_energy_ev(velocities: np.ndarray, masses_kg: np.ndarray) -> np.ndarray:
    """Return kinetic energy from velocity vectors in electron-volts."""
    velocities = np.asarray(velocities, dtype=float)
    masses = np.asarray(masses_kg, dtype=float)
    v2 = np.sum(velocities * velocities, axis=-1)
    return 0.5 * v2 * masses / ELEMENTARY_CHARGE_C


def tof_mass_to_charge_da_per_e(
    flight_time_s: np.ndarray,
    acceleration_voltage_V: float,
    flight_distance_m: float,
) -> np.ndarray:
    """Convert TOF flight times to m/z in Da/e using the constant-field TOF relation."""
    t = np.asarray(flight_time_s, dtype=float)
    if acceleration_voltage_V <= 0.0:
        raise ValueError("acceleration_voltage_V must be positive")
    if flight_distance_m <= 0.0:
        raise ValueError("flight_distance_m must be positive")
    return (2.0 * acceleration_voltage_V * ELEMENTARY_CHARGE_C * (t / flight_distance_m) ** 2) / ATOMIC_MASS_KG


def maxwell_speed_pdf(speed: np.ndarray, temperature_K: float, mass_kg: float) -> np.ndarray:
    """Maxwell-Boltzmann speed distribution."""
    speed = np.asarray(speed, dtype=float)
    coeff = np.sqrt(2.0 / np.pi) * (mass_kg / (BOLTZMANN_CONSTANT * temperature_K)) ** 1.5
    return coeff * speed**2 * np.exp(-mass_kg * speed**2 / (2.0 * BOLTZMANN_CONSTANT * temperature_K))


def field_strength_from_en_td(en_td: float, pressure_pa: float, temperature_K: float) -> float:
    """Convert reduced field E/N in Td to electric field strength E in V/m."""
    number_density = pressure_pa / (BOLTZMANN_CONSTANT * temperature_K)
    return float(en_td) * TOWNSEND_TO_V_M2 * number_density


def reduced_mobility_from_mobility(
    mobility_m2_Vs: np.ndarray,
    pressure_pa: float,
    temperature_K: float,
    p0_pa: float = STANDARD_PRESSURE_PA,
    t0_K: float = STANDARD_TEMPERATURE_K,
) -> np.ndarray:
    """Convert mobility K to reduced mobility K0."""
    mobility = np.asarray(mobility_m2_Vs, dtype=float)
    return mobility * (pressure_pa / p0_pa) * (t0_K / temperature_K)


def gaussian(x: np.ndarray, amplitude: float, mean: float, sigma: float) -> np.ndarray:
    """Simple Gaussian curve used for histogram overlays."""
    return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def histogram_bin_edges(values: np.ndarray, bins: int) -> np.ndarray:
    """Build stable histogram edges for one or more scalar observations."""
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        raise ValueError("Cannot build histogram edges for empty input.")
    if values.size > 1:
        return np.linspace(values.min(), values.max(), bins + 1)
    return np.linspace(values.min() * 0.99, values.max() * 1.01 + 1e-12, bins + 1)


def event_time_bin_edges(values: np.ndarray, bins: int, log_bins: bool) -> np.ndarray:
    """Build elimination/event-time histogram edges with optional log scaling."""
    values = np.asarray(values, dtype=float)
    if values.size == 0:
        raise ValueError("Cannot build time histogram edges for empty input.")
    t_max = float(values.max())
    if log_bins and t_max > 0.0:
        positive = values[values > 0.0]
        t_min = float(positive.min()) if positive.size else 1e-12
        return np.geomspace(t_min, t_max, num=bins + 1)
    return np.linspace(0.0, t_max, num=bins + 1)


def fit_gaussian_histogram(values: np.ndarray, bins: int) -> tuple[dict[str, float] | None, np.ndarray, np.ndarray]:
    """Fit a Gaussian to histogram counts derived from the provided values."""
    from scipy.optimize import curve_fit  # Lazy import; not all scripts need fitting

    values = np.asarray(values, dtype=float)
    if values.size < 3:
        return None, np.array([]), np.array([])

    edges = histogram_bin_edges(values, bins)
    hist, _ = np.histogram(values, bins=edges)
    bin_centers = 0.5 * (edges[:-1] + edges[1:])

    if hist.size < 3 or hist.max() == 0:
        return None, bin_centers, hist

    try:
        amplitude0 = float(hist.max())
        mean0 = float(values.mean())
        sigma0 = float(values.std())
        params, _ = curve_fit(gaussian, bin_centers, hist, p0=[amplitude0, mean0, sigma0], maxfev=10000)
        return {
            "a": float(params[0]),
            "mu": float(params[1]),
            "sigma": float(params[2]),
        }, bin_centers, hist
    except Exception:
        return None, bin_centers, hist

#!/usr/bin/env python3
"""Backward-compatible facade for shared analysis helpers.

New code should import from the focused modules:
- ``analysis.trajectory`` for trajectory HDF5 access and selection
- ``analysis.physics`` for constants, conversions, histograms, and fitting
"""

from __future__ import annotations

from typing import Sequence

from analysis.physics import (
    BOLTZMANN_CONSTANT,
    STANDARD_PRESSURE_PA,
    STANDARD_TEMPERATURE_K,
    TOWNSEND_TO_V_M2,
    event_time_bin_edges,
    field_strength_from_en_td,
    fit_gaussian_histogram,
    gaussian,
    histogram_bin_edges,
    maxwell_speed_pdf,
    reduced_mobility_from_mobility,
    temperature_from_velocities,
    temperature_series_from_velocities,
)
from analysis.trajectory import (
    DomainBounds,
    TrajectorySelection,
    active_ion_mask,
    decode_hdf5_string,
    ensure_nonempty_selection,
    find_arrival_times,
    load_death_times_for_indices,
    load_positions_subset,
    load_species_ids,
    load_trajectory_selection,
    normalize_species_filter,
    open_trajectory,
    read_death_times_for_indices,
    read_domain_bounds,
    read_domain_environment,
    read_domain_geometry,
    read_positions_subset,
    read_species_ids,
    read_species_mass_map,
    read_trajectory_selection,
    select_ion_indices,
    select_ions_from_trajectory,
    species_for_indices,
    species_masses_for_selection,
)


def _decode(val) -> str:
    """Backward-compatible private alias used by older scripts."""
    return decode_hdf5_string(val)


def species_color_map(species: Sequence[str]) -> dict[str, str]:
    """Assign distinct-ish matplotlib colors to a species list."""
    import matplotlib.pyplot as plt  # Lazy import to keep dependencies local

    unique = sorted(set(species))
    cmap = plt.cm.get_cmap("tab20", len(unique))
    return {sp: cmap(i) for i, sp in enumerate(unique)}


__all__ = [
    "BOLTZMANN_CONSTANT",
    "DomainBounds",
    "STANDARD_PRESSURE_PA",
    "STANDARD_TEMPERATURE_K",
    "TOWNSEND_TO_V_M2",
    "TrajectorySelection",
    "active_ion_mask",
    "decode_hdf5_string",
    "ensure_nonempty_selection",
    "event_time_bin_edges",
    "field_strength_from_en_td",
    "find_arrival_times",
    "fit_gaussian_histogram",
    "gaussian",
    "histogram_bin_edges",
    "load_death_times_for_indices",
    "load_positions_subset",
    "load_species_ids",
    "load_trajectory_selection",
    "maxwell_speed_pdf",
    "normalize_species_filter",
    "open_trajectory",
    "read_death_times_for_indices",
    "read_domain_bounds",
    "read_domain_environment",
    "read_domain_geometry",
    "read_positions_subset",
    "read_species_ids",
    "read_species_mass_map",
    "read_trajectory_selection",
    "reduced_mobility_from_mobility",
    "select_ion_indices",
    "select_ions_from_trajectory",
    "species_color_map",
    "species_for_indices",
    "species_masses_for_selection",
    "temperature_from_velocities",
    "temperature_series_from_velocities",
]

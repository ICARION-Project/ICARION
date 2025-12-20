#!/usr/bin/env python3
"""
Shared HDF5 helpers for validation/analysis scripts.

Handles schema changes (species_id_indices) and embedded input blobs
under /metadata/reproducibility/input_blobs.
"""

from __future__ import annotations

import json
import numpy as np


def load_species_ids(handle):
    """
    Return species labels as a numpy array of strings with the same shape as
    trajectory species datasets. Supports legacy string dataset and new
    species_id_indices + metadata/species/names.
    """
    if "trajectory/species_ids" in handle:
        raw = handle["trajectory/species_ids"][:]
        # bytes -> str
        if raw.dtype.kind in ("S", "O"):
            return np.vectorize(lambda x: x.decode("utf-8") if isinstance(x, (bytes, bytearray)) else str(x))(raw)
        # Already numeric: try to map through metadata
        if raw.dtype.kind in ("i", "u"):
            names = _read_species_names(handle)
            return np.vectorize(lambda idx: names[idx] if 0 <= idx < len(names) else str(idx))(raw)
        return raw.astype(str)

    if "trajectory/species_id_indices" in handle:
        idx = handle["trajectory/species_id_indices"][:]
        names = _read_species_names(handle)
        return np.vectorize(lambda i: names[i] if 0 <= i < len(names) else str(i))(idx)

    raise KeyError("No species IDs found (trajectory/species_ids or species_id_indices missing)")


def load_config_json(handle):
    """Return embedded config JSON dict if present, else None."""
    path = "/metadata/reproducibility/input_blobs/config_json"
    if path in handle:
        try:
            return json.loads(handle[path][()].decode("utf-8"))
        except Exception:
            return None
    return None


def load_species_db_json(handle):
    """Return embedded species DB JSON dict if present, else None."""
    path = "/metadata/reproducibility/input_blobs/species_db_json"
    if path in handle:
        try:
            return json.loads(handle[path][()].decode("utf-8"))
        except Exception:
            return None
    return None


def load_reaction_db_json(handle):
    """Return embedded reaction DB JSON dict if present, else None."""
    path = "/metadata/reproducibility/input_blobs/reaction_db_json"
    if path in handle:
        try:
            return json.loads(handle[path][()].decode("utf-8"))
        except Exception:
            return None
    return None


def _read_species_names(handle):
    if "metadata/species/names" in handle:
        return [s.decode("utf-8") if isinstance(s, (bytes, bytearray)) else str(s) for s in handle["metadata/species/names"][:]]
    return []

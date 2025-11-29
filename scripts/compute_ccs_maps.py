#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# SPDX-FileCopyrightText: 2025 ICARION Project Contributors

"""
Compute gas-specific CCS maps from a reference CCS for HSS/EHSS models.

Given a species database JSON and a reference gas + CCS, derive
gas-specific CCS values using kinetic diameters:

    CCS_gas = π * (r_ion + r_gas)^2

where r_ion is inferred from the reference CCS and gas radius.

Usage:
    python scripts/compute_ccs_maps.py --input species.json --output species_with_ccs.json \
        --ref-gas He --ref-ccs-A2 110.0

Flags:
    --model {HSS,EHSS}   Model label to store (default: HSS)
    --override            Overwrite existing CCS_* maps if present
"""

import argparse
import json
import math
import sys
from copy import deepcopy

# Gas radii (meters) matching utils/constants.h
GAS_RADII_M = {
    "He": 1.3e-10,
    "Ar": 1.70e-10,
    "CO2": 1.65e-10,
    "Ne": 1.38e-10,
    "N2": 1.82e-10,
    "O2": 1.73e-10,
    "H2O": 1.58e-10,
}

M2_TO_A2 = 1e20
A2_TO_M2 = 1e-20


def derive_sigma_for_gas(sigma_ref_m2: float, gas_ref: str, gas_target: str) -> float:
    if gas_ref not in GAS_RADII_M or gas_target not in GAS_RADII_M:
        return 0.0
    r_ref = GAS_RADII_M[gas_ref]
    r_ion = max(0.0, math.sqrt(max(sigma_ref_m2, 0.0) / math.pi) - r_ref)
    r_target = GAS_RADII_M[gas_target]
    return math.pi * (r_ion + r_target) ** 2


def process_species_db(data: dict, ref_gas: str, ref_ccs_a2: float, model: str, override: bool) -> dict:
    out = deepcopy(data)
    species_dict = out.get("species", {})

    sigma_ref_m2 = ref_ccs_a2 * A2_TO_M2
    supported_gases = list(GAS_RADII_M.keys())

    for sid, props in species_dict.items():
        # Skip neutrals without CCS
        if props.get("charge", 0) == 0:
            continue

        # Add reference info
        key_model = f"CCS_{model}"
        if not override and key_model in props:
            continue

        ccs_map = {}
        for g in supported_gases:
            sigma_g = derive_sigma_for_gas(sigma_ref_m2, ref_gas, g)
            if sigma_g > 0.0:
                ccs_map[g] = sigma_g * M2_TO_A2

        props["CCS_reference_gas"] = ref_gas
        props["CCS_model"] = model
        props[key_model] = ccs_map

    out["species"] = species_dict
    return out


def main():
    parser = argparse.ArgumentParser(description="Compute gas-specific CCS maps from reference CCS.")
    parser.add_argument("--input", required=True, help="Input species JSON (SpeciesLoader format)")
    parser.add_argument("--output", required=True, help="Output JSON path")
    parser.add_argument("--ref-gas", required=True, help="Reference gas ID (e.g., He, N2, O2)")
    parser.add_argument("--ref-ccs-A2", type=float, required=True, help="Reference CCS in Å^2")
    parser.add_argument("--model", choices=["HSS", "EHSS"], default="HSS", help="CCS model label")
    parser.add_argument("--override", action="store_true", help="Overwrite existing CCS maps if present")
    args = parser.parse_args()

    try:
        with open(args.input, "r") as f:
            data = json.load(f)
    except Exception as exc:
        print(f"Failed to read input JSON: {exc}", file=sys.stderr)
        sys.exit(1)

    result = process_species_db(data, args.ref_gas, args.ref_ccs_A2, args.model, args.override)

    try:
        with open(args.output, "w") as f:
            json.dump(result, f, indent=2)
        print(f"✓ Wrote updated species database to {args.output}")
    except Exception as exc:
        print(f"Failed to write output JSON: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()

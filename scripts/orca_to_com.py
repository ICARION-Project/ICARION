#!/usr/bin/env python3
"""Extract final ORCA geometries and shift them to the center-of-mass frame."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
import re
from typing import Dict, Iterable, List, Sequence, Tuple

# Basic element properties used across the existing molecule JSON files.
ATOMIC_MASSES: Dict[str, float] = {
    "H": 1.008,
    "C": 12.011,
    "N": 14.007,
    "O": 15.999,
}

LJ_PARAMS: Dict[str, Dict[str, float]] = {
    # Values mirror the existing data/molecules/*.json entries.
    "H": {"sigma": 2.261, "epsilon": 0.059579817756640736},
    "C": {"sigma": 3.0126, "epsilon": 0.21252132},
    "N": {"sigma": 3.3473, "epsilon": 0.2361348},
    "O": {"sigma": 2.4344, "epsilon": 0.1034324933362249},
}


def _parse_final_cartesians(lines: Sequence[str]) -> List[Tuple[str, List[float]]]:
    """Return the last 'CARTESIAN COORDINATES (ANGSTROEM)' block in the file."""
    header = "CARTESIAN COORDINATES (ANGSTROEM)"
    blocks: List[List[Tuple[str, List[float]]]] = []
    i = 0
    while i < len(lines):
        if header in lines[i]:
            i += 1
            while i < len(lines) and lines[i].strip().startswith("-"):
                i += 1
            atoms: List[Tuple[str, List[float]]] = []
            while i < len(lines):
                line = lines[i].strip()
                if not line:
                    break
                parts = line.split()
                if len(parts) < 4 or not parts[0][0].isalpha():
                    break
                try:
                    x, y, z = (float(parts[1]), float(parts[2]), float(parts[3]))
                except ValueError:
                    break
                atoms.append((parts[0], [x, y, z]))
                i += 1
            if atoms:
                blocks.append(atoms)
        i += 1
    if not blocks:
        raise ValueError("No 'CARTESIAN COORDINATES (ANGSTROEM)' block found.")
    return blocks[-1]


def _center_of_mass(atoms: Iterable[Tuple[str, Sequence[float]]]) -> List[float]:
    total_mass = 0.0
    weighted = [0.0, 0.0, 0.0]
    for element, (x, y, z) in atoms:
        if element not in ATOMIC_MASSES:
            raise KeyError(f"Missing atomic mass for element '{element}'.")
        mass = ATOMIC_MASSES[element]
        total_mass += mass
        weighted[0] += mass * x
        weighted[1] += mass * y
        weighted[2] += mass * z
    if total_mass == 0:
        raise ValueError("Total mass is zero; cannot compute center of mass.")
    return [value / total_mass for value in weighted]


def _to_com_frame(
    atoms: Iterable[Tuple[str, Sequence[float]]],
) -> List[Tuple[str, List[float]]]:
    atoms_list = list(atoms)
    com = _center_of_mass(atoms_list)
    shifted: List[Tuple[str, List[float]]] = []
    for element, (x, y, z) in atoms_list:
        shifted.append((element, [x - com[0], y - com[1], z - com[2]]))
    return shifted


def _parse_final_mulliken_charges(lines: Sequence[str]) -> List[float]:
    """Return the last Mulliken atomic charges block."""
    header = "MULLIKEN ATOMIC CHARGES"
    blocks: List[List[float]] = []
    i = 0
    charge_pattern = re.compile(r"^\s*\d+\s+[A-Za-z]+\s*:\s*([-+]?\d*\.\d+(?:[Ee][-+]?\d+)?)")
    while i < len(lines):
        if header in lines[i]:
            i += 1
            while i < len(lines) and lines[i].strip().startswith("-"):
                i += 1
            charges: List[float] = []
            while i < len(lines):
                line = lines[i]
                if "Sum of atomic charges" in line:
                    break
                match = charge_pattern.match(line)
                if not match:
                    if line.strip() == "":
                        break
                    # Skip non-matching detail lines until the block ends.
                    i += 1
                    continue
                charges.append(float(match.group(1)))
                i += 1
            if charges:
                blocks.append(charges)
        i += 1
    if not blocks:
        raise ValueError("No 'MULLIKEN ATOMIC CHARGES' block found.")
    return blocks[-1]


def _atom_record(element: str, pos: Sequence[float]) -> Dict[str, object]:
    if element not in ATOMIC_MASSES:
        raise KeyError(f"Unknown element '{element}' (no mass available).")
    if element not in LJ_PARAMS:
        raise KeyError(f"Unknown element '{element}' (no LJ parameters available).")
    return {
        "element": element,
        "pos": [float(coord) for coord in pos],
        "mass_u": ATOMIC_MASSES[element],
        "partial_charge_e": 0.0,
        "LJ_sigma_angstrom": LJ_PARAMS[element]["sigma"],
        "LJ_epsilon_eV": LJ_PARAMS[element]["epsilon"],
    }


def convert_orca_file(path: Path, name: str) -> Dict[str, object]:
    text = path.read_text().splitlines()
    raw_atoms = _parse_final_cartesians(text)
    mulliken_charges = _parse_final_mulliken_charges(text)
    if len(mulliken_charges) != len(raw_atoms):
        raise ValueError(
            f"Charge count ({len(mulliken_charges)}) does not match atom count ({len(raw_atoms)})."
        )
    shifted_atoms = _to_com_frame(raw_atoms)
    atom_entries = []
    for (elem, coords), charge in zip(shifted_atoms, mulliken_charges):
        atom = _atom_record(elem, coords)
        atom["partial_charge_e"] = float(charge)
        atom_entries.append(atom)
    return {
        "molecule": {
            "name": name,
            "diameter_m": None,
            "CCS_m2": None,
            "atoms": atom_entries,
        }
    }


def derive_name(path: Path, suffix: str = "") -> str:
    stem = path.stem
    if stem.endswith("_ORCA"):
        stem = stem[:-5]
    return stem + suffix


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Extract the final ORCA geometry, shift to the center-of-mass frame, "
            "and emit a molecule JSON payload."
        )
    )
    parser.add_argument("orca_out", nargs="+", type=Path, help="ORCA .out file(s)")
    parser.add_argument(
        "--name",
        help="Explicit molecule name (only used when a single ORCA file is given).",
    )
    parser.add_argument(
        "--suffix",
        default="",
        help="Suffix appended to the derived name (before writing the file).",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("data/molecules"),
        help="Where to write the JSON files.",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print JSON to stdout instead of writing files.",
    )
    args = parser.parse_args()

    if args.name and len(args.orca_out) != 1:
        raise SystemExit("--name can only be used with a single ORCA file.")

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for orca_file in args.orca_out:
        molecule_name = args.name or derive_name(orca_file, args.suffix)
        payload = convert_orca_file(orca_file, molecule_name)
        if args.dry_run:
            print(f"# {molecule_name}")
            print(json.dumps(payload, indent=2))
        else:
            out_path = args.output_dir / f"{molecule_name}.json"
            out_path.write_text(json.dumps(payload, indent=2))
            print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()

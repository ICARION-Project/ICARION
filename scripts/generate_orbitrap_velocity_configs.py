#!/usr/bin/env python3
"""
Generate three Orbitrap configs to probe trapping resilience vs velocity spread:
1) Speed spread only (Gaussian velocity, no angular spread)
2) Cone spread only (kinetic with spread_angle_deg)
3) Speed + cone spread (Gaussian velocity with transverse std set by cone)

Defaults use species from tmp/species_database_lqit.json (ReserpineH+, CaffeineH+).
"""

from __future__ import annotations

import argparse
import json
from math import tan, radians, sqrt
from pathlib import Path
from typing import Dict, List

E_CHARGE = 1.602176634e-19  # C
AMU = 1.66053906660e-27      # kg


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate Orbitrap velocity spread configs.")
    p.add_argument("--species-database", type=Path, default=Path("tmp/species_database_lqit.json"),
                   help="Species database path (used to get masses).")
    p.add_argument("--species", nargs="+", default=["ReserpineH+", "CaffeineH+"],
                   help="Species IDs to include.")
    p.add_argument("--output-dir", type=Path, default=Path("tmp/orbitrap_vel_configs"),
                   help="Directory for generated configs.")
    p.add_argument("--energy-eV", type=float, default=1600.0,
                   help="Injection kinetic energy for kinetic/cone cases [eV].")
    p.add_argument("--speed-frac", type=float, default=0.05,
                   help="Relative std for speed spread (e.g., 0.05 = 5%).")
    p.add_argument("--cone-deg", type=float, default=2.0,
                   help="Full opening angle (degrees) for cone spread.")
    p.add_argument("--radial-start-V", type=float, default=1700.0,
                   help="Orbitrap radial_V start value for ramp [V].")
    p.add_argument("--radial-end-V", type=float, default=1732.0,
                   help="Orbitrap radial_V end value for ramp [V].")
    p.add_argument("--radial-ramp-us", type=float, default=1.0,
                   help="Ramp duration in microseconds.")
    p.add_argument("--ion-count", type=int, default=50, help="Ions per species.")
    p.add_argument("--total-time", type=float, default=2e-3, help="Simulation time [s].")
    p.add_argument("--dt", type=float, default=1e-9, help="Timestep [s].")
    return p.parse_args()


def load_masses(db_path: Path, species_ids: List[str]) -> Dict[str, float]:
    with db_path.open() as fh:
        db = json.load(fh)
    masses = {}
    for sp in species_ids:
        if sp not in db.get("species", {}):
            raise KeyError(f"Species '{sp}' not found in {db_path}")
        masses[sp] = db["species"][sp]["mass_amu"]
    return masses


def speed_from_energy(energy_eV: float, mass_amu: float) -> float:
    return sqrt(2 * energy_eV * E_CHARGE / (mass_amu * AMU))


def base_config(template_species: List[Dict], sim_time: float, dt: float,
                radial_start_V: float, radial_end_V: float, ramp_us: float) -> Dict:
    return {
        "title": "Orbitrap velocity spread experiment",
        "species_database": str(Path("tmp/species_database_lqit.json").resolve()),
        "simulation": {
            "total_time_s": sim_time,
            "dt_s": dt,
            "write_interval": 100,
            "integrator": "RK4",
            "enable_gpu": False,
            "enable_openmp": True,
            "rng_seed": 42
        },
        "physics": {
            "collision_model": "NoCollisions",
            "enable_space_charge": False
        },
        "output": {
            "folder": str(Path("results/orbitrap_vel").resolve()),
            "trajectory_file": None,  # set per variant
            "print_progress": True
        },
        "ions": {
            "species": template_species
        },
        "domains": [
            {
                "name": "orbitrap_trap",
                "instrument": "Orbitrap",
                "geometry": {
                    "origin_m": [0.0, 0.0, 0.0],
                    "length_m": 0.04,
                    "radius_in_m": 0.006,
                    "radius_out_m": 0.015,
                    "radius_char_m": 0.022
                },
                "env": {
                    "temperature_K": 300.0,
                    "pressure_Pa": 1e-7,
                    "gas_species": "He",
                    "gas_velocity_m_s": [0.0, 0.0, 0.0]
                },
                "fields": {
                    "DC": {
                        "axial_V": 0.0,
                        "radial_V": {
                            "type": "linear",
                            "start": radial_start_V,
                            "end": radial_end_V,
                            "start_time_s": 0.0,
                            "end_time_s": ramp_us * 1e-6,
                            "clamp": True
                        }
                    }
                }
            }
        ]
    }


def make_species_entries(species_ids: List[str], masses: Dict[str, float],
                         count: int, v_mean: Dict[str, float],
                         speed_frac: float, cone_deg: float,
                         mode: str) -> List[Dict]:
    entries = []
    for sp in species_ids:
        if mode == "speed_spread":
            vy = v_mean[sp]
            std_y = speed_frac * vy
            cone_std = 0.0
            entries.append({
                "species_id": sp,
                "count": count,
                "position": {"type": "point", "center": [0.009, 0.0, 0.006]},
                "velocity": {
                    "type": "gaussian",
                    "mean": [0.0, vy, 0.0],
                    "std": [cone_std, std_y, cone_std]
                }
            })
        elif mode == "cone":
            entries.append({
                "species_id": sp,
                "count": count,
                "position": {"type": "point", "center": [0.009, 0.0, 0.006]},
                "velocity": {
                    "type": "kinetic",
                    "energy_eV": args.energy_eV,
                    "direction": [0.0, 1.0, 0.0],
                    "spread_angle_deg": cone_deg
                }
            })
        elif mode == "both":
            vy = v_mean[sp]
            std_y = speed_frac * vy
            cone_std = vy * tan(radians(cone_deg / 2.0))
            entries.append({
                "species_id": sp,
                "count": count,
                "position": {"type": "point", "center": [0.009, 0.0, 0.006]},
                "velocity": {
                    "type": "gaussian",
                    "mean": [0.0, vy, 0.0],
                    "std": [cone_std, std_y, cone_std]
                }
            })
        else:
            raise ValueError(f"Unknown mode {mode}")
    return entries


if __name__ == "__main__":
    args = parse_args()
    masses = load_masses(args.species_database, args.species)
    v_mean = {sp: speed_from_energy(args.energy_eV, masses[sp]) for sp in args.species}

    modes = [
        ("speed_spread", "orbitrap_speed_spread.json"),
        ("cone", "orbitrap_cone.json"),
        ("both", "orbitrap_speed_spread_cone.json"),
    ]

    args.output_dir.mkdir(parents=True, exist_ok=True)

    for mode, fname in modes:
        species_entries = make_species_entries(
            args.species, masses, args.ion_count, v_mean,
            args.speed_frac, args.cone_deg, mode
        )
        cfg = base_config(
            species_entries,
            args.total_time,
            args.dt,
            args.radial_start_V,
            args.radial_end_V,
            args.radial_ramp_us,
        )
        cfg["output"]["trajectory_file"] = fname.replace(".json", ".h5")

        out_path = args.output_dir / fname
        with out_path.open("w") as fh:
            json.dump(cfg, fh, indent=2)
        print(f"Wrote {out_path}")

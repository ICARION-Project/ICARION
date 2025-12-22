#!/usr/bin/env python3
"""Reaction kinetics validation harness.

Runs the five agreed scenarios for reaction physics validation:
1. Reversible equilibrium (H3O+ <-> PentanalH+)
2. First-order unimolecular decay
3. Competing bimolecular channels (He vs Pentanal)
4. Sequential chain (H3O+ -> PentanalH+ -> CaffeineH+)
5. Temperature-dependent Arrhenius sweep

Each scenario points to a dedicated config in validation/configs/physics/reactions.
The script only orchestrates executions and logs analytic expectations; post-run
analysis still happens via analyze_reactions.py once time-resolved outputs land
in the HDF5 files.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Dict, Iterable, List

# Physical constants
K_B = 1.380649e-23  # Boltzmann constant [J/K]
E_CHARGE = 1.602176634e-19  # Elementary charge [C]

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
CONFIG_DIR = PROJECT_ROOT / "validation" / "configs" / "physics" / "reactions"
RESULTS_ROOT = PROJECT_ROOT / "validation" / "results" / "physics" / "reactions"
LOG_DIR = PROJECT_ROOT / "validation" / "logs"
LOG_DIR.mkdir(parents=True, exist_ok=True)
RESULTS_ROOT.mkdir(parents=True, exist_ok=True)
LOG_FILE = LOG_DIR / "REACTION_KINETICS_VALIDATION.txt"
ICARION_BIN_DEFAULT = PROJECT_ROOT / "build" / "src" / "icarion_main"


@dataclass
class ScenarioRun:
    """Single configuration within a scenario."""

    label: str
    config_path: Path
    metadata: Dict[str, float] = field(default_factory=dict)


@dataclass
class Scenario:
    """Scenario definition bundling configs and metadata."""

    key: str
    description: str
    runs: List[ScenarioRun]
    metadata: Dict[str, float] = field(default_factory=dict)


SCENARIOS: Dict[str, Scenario] = {
    "equilibrium": Scenario(
        key="equilibrium",
        description="Reversible first-order pair expected to settle at a 75/25 split.",
        runs=[
            ScenarioRun(
                label="equilibrium",
                config_path=CONFIG_DIR / "reversible_equilibrium_h3o_pentanal.json",
                metadata={"k_forward": 1800.0, "k_reverse": 600.0, "n0": 1000.0},
            )
        ],
    ),
    "first_order": Scenario(
        key="first_order",
        description="Single-channel decay to benchmark first-order kinetics (k=5e3 s^-1).",
        runs=[
            ScenarioRun(
                label="first_order_decay",
                config_path=CONFIG_DIR / "first_order_decay_reference.json",
                metadata={"k": 5000.0, "n0": 1000.0},
            )
        ],
    ),
    "competing": Scenario(
        key="competing",
        description="Bimolecular competition between Pentanal and He channels.",
        runs=[
            ScenarioRun(
                label="he_vs_pentanal",
                config_path=CONFIG_DIR / "competing_channels_he_vs_pentanal.json",
                metadata={
                    "k_he": 1.25e-15,
                    "k_pentanal": 3.75e-15,
                    "n_he": 6.760105814970592e22,
                    "n_pentanal": 2.897188206415968e22,
                },
            )
        ],
    ),
    "sequential": Scenario(
        key="sequential",
        description="Two-step sequential chain ending in CaffeineH+.",
        runs=[
            ScenarioRun(
                label="chain",
                config_path=CONFIG_DIR / "sequential_chain_three_species.json",
                metadata={"k1": 1200.0, "k2": 800.0, "n0": 1000.0},
            )
        ],
    ),
    "arrhenius": Scenario(
        key="arrhenius",
        description="Arrhenius sweep for H3O+ -> PentanalH+ across 250/300/350 K.",
        metadata={"prefactor": 2.0e5, "activation_energy_eV": 0.08},
        runs=[
            ScenarioRun(
                label="250K",
                config_path=CONFIG_DIR / "arrhenius_h3o_to_pentanal_250K.json",
                metadata={"temperature_K": 250.0},
            ),
            ScenarioRun(
                label="300K",
                config_path=CONFIG_DIR / "arrhenius_h3o_to_pentanal_300K.json",
                metadata={"temperature_K": 300.0},
            ),
            ScenarioRun(
                label="350K",
                config_path=CONFIG_DIR / "arrhenius_h3o_to_pentanal_350K.json",
                metadata={"temperature_K": 350.0},
            ),
        ],
    ),
}


def log(msg: str, file_handle) -> None:
    """Print to stdout and duplicate into the log file."""

    print(msg)
    file_handle.write(msg + "\n")
    file_handle.flush()


def ensure_config_exists(run: ScenarioRun) -> None:
    if not run.config_path.exists():
        raise FileNotFoundError(f"Missing config: {run.config_path}")


def ensure_output_directory(config_path: Path) -> None:
    with open(config_path, "r", encoding="utf-8") as handle:
        cfg = json.load(handle)
    out_cfg = cfg.get("output", {})
    folder = out_cfg.get("folder")
    if not folder:
        return
    folder_path = Path(folder)
    if not folder_path.is_absolute():
        folder_path = PROJECT_ROOT / folder
    folder_path.mkdir(parents=True, exist_ok=True)


def run_simulation(config_path: Path, icarion_bin: Path, dry_run: bool, logf) -> bool:
    rel_path = config_path.relative_to(PROJECT_ROOT)
    cmd = [str(icarion_bin), str(config_path)]
    if dry_run:
        log(f"  [dry-run] {' '.join(cmd)}", logf)
        return True

    log(f"  → Running {rel_path}", logf)
    result = subprocess.run(
        cmd,
        cwd=PROJECT_ROOT,
        text=True,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        log("  ✗ Simulation failed", logf)
        err = (result.stderr or "").strip()
        log(f"    STDERR: {err}", logf)
        return False

    log("  ✓ Simulation complete", logf)
    return True


def equilibrium_notes(run: ScenarioRun) -> List[str]:
    kf = run.metadata["k_forward"]
    kr = run.metadata["k_reverse"]
    total = run.metadata.get("n0", 1000.0)
    tau = 1.0 / (kf + kr)
    pent_eq = total * kf / (kf + kr)
    h3o_eq = total - pent_eq
    return [
        f"kf={kf:.1f} s^-1, kr={kr:.1f} s^-1",
        f"Time constant tau={(tau * 1e3):.3f} ms",
        f"Equilibrium populations → H3O+={h3o_eq:.1f}, PentanalH+={pent_eq:.1f}",
    ]


def first_order_notes(run: ScenarioRun) -> List[str]:
    k = run.metadata["k"]
    tau = 1.0 / k
    half_life = math.log(2.0) / k
    return [
        f"Decay constant k={k:.1f} s^-1",
        f"Tau={(tau * 1e6):.2f} µs, t1/2={(half_life * 1e6):.2f} µs",
    ]


def competing_notes(run: ScenarioRun) -> List[str]:
    k_he_eff = run.metadata["k_he"] * run.metadata["n_he"]
    k_pent_eff = run.metadata["k_pentanal"] * run.metadata["n_pentanal"]
    total = k_he_eff + k_pent_eff
    branch_pent = k_pent_eff / total
    branch_he = k_he_eff / total
    return [
        f"Effective rates: Pentanal={k_pent_eff:.3e} s^-1, He={k_he_eff:.3e} s^-1",
        f"Expected branching → PentanalH+={branch_pent*100:.1f}%, 2,6-DTBP={branch_he*100:.1f}%",
    ]


def sequential_notes(run: ScenarioRun) -> List[str]:
    k1 = run.metadata["k1"]
    k2 = run.metadata["k2"]
    tau1 = 1.0 / k1
    tau2 = 1.0 / k2
    if k1 != k2:
        t_peak = math.log(k1 / k2) / (k1 - k2)
        peak_str = f"PentanalH+ peak near {t_peak * 1e3:.3f} ms"
    else:
        peak_str = "k1 == k2 → peak undefined"
    return [
        f"k1={k1:.1f} s^-1, tau1={(tau1 * 1e3):.3f} ms",
        f"k2={k2:.1f} s^-1, tau2={(tau2 * 1e3):.3f} ms",
        peak_str,
    ]


def arrhenius_notes(scenario: Scenario) -> List[str]:
    prefactor = scenario.metadata["prefactor"]
    ea_j = scenario.metadata["activation_energy_eV"] * E_CHARGE
    notes = [
        f"Prefactor A={prefactor:.2e} s^-1, Ea={scenario.metadata['activation_energy_eV']:.3f} eV",
        "Predicted rates per temperature:",
    ]
    for run in scenario.runs:
        temp = run.metadata["temperature_K"]
        k_t = prefactor * math.exp(-ea_j / (K_B * temp))
        notes.append(f"  {run.label}: k(T)={k_t:.3e} s^-1, tau={(1.0 / k_t * 1e6):.2f} µs")
    return notes


def scenario_notes(scenario: Scenario) -> Iterable[str]:
    if scenario.key == "equilibrium":
        return equilibrium_notes(scenario.runs[0])
    if scenario.key == "first_order":
        return first_order_notes(scenario.runs[0])
    if scenario.key == "competing":
        return competing_notes(scenario.runs[0])
    if scenario.key == "sequential":
        return sequential_notes(scenario.runs[0])
    if scenario.key == "arrhenius":
        return arrhenius_notes(scenario)
    return []


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run reaction kinetics validation scenarios.")
    parser.add_argument(
        "--scenarios",
        nargs="+",
        default=["all"],
        choices=["all"] + sorted(SCENARIOS.keys()),
        help="Subset of scenarios to execute",
    )
    parser.add_argument(
        "--icarion-bin",
        type=Path,
        default=ICARION_BIN_DEFAULT,
        help="Path to the icarion_main executable",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without running simulations",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.icarion_bin.exists() and not args.dry_run:
        raise SystemExit(f"icarion_main not found at {args.icarion_bin}")

    selected = list(SCENARIOS.keys()) if "all" in args.scenarios else args.scenarios

    start = datetime.now()
    with open(LOG_FILE, "w", encoding="utf-8") as logf:
        log("=" * 80, logf)
        log("REACTION KINETICS VALIDATION", logf)
        log("=" * 80, logf)
        log(f"Project root: {PROJECT_ROOT}", logf)
        log(f"Icarion binary: {args.icarion_bin}", logf)
        log(f"Dry run: {args.dry_run}", logf)
        log(f"Start: {start.isoformat()}", logf)
        log("", logf)

        for key in selected:
            scenario = SCENARIOS[key]
            log(f"--- Scenario: {scenario.key} ---", logf)
            log(scenario.description, logf)
            for note in scenario_notes(scenario):
                log(f"  - {note}", logf)
            log("", logf)

            for run in scenario.runs:
                ensure_config_exists(run)
                ensure_output_directory(run.config_path)
                success = run_simulation(run.config_path, args.icarion_bin, args.dry_run, logf)
                if not success:
                    log("Aborting remaining runs due to failure.", logf)
                    return
            log("", logf)

        log("All requested scenarios finished.", logf)
        log("Next step: use analyze_reactions.py for quantitative fits once time-resolved species data is available.", logf)
        log(f"End: {datetime.now().isoformat()}", logf)


if __name__ == "__main__":
    main()

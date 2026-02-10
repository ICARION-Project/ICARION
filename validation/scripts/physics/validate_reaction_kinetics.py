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
import os

# Physical constants
K_B = 1.380649e-23  # Boltzmann constant [J/K]
E_CHARGE = 1.602176634e-19  # Elementary charge [C]

# Paths
PROJECT_ROOT = Path(__file__).resolve().parents[3]
CONFIG_DIR = PROJECT_ROOT / "validation" / "configs" / "physics" / "reactions"
ICARION_BIN_DEFAULT = PROJECT_ROOT / "build" / "src" / "icarion_main"

RUN_DIR = os.environ.get("ICARION_VALIDATION_RUN_DIR")
DEFAULT_RESULTS_ROOT = PROJECT_ROOT / "validation" / "results" / "physics" / "reactions"
DEFAULT_LOG_DIR = PROJECT_ROOT / "validation" / "logs"
if RUN_DIR:
    RUN_DIR = Path(RUN_DIR)
    DEFAULT_RESULTS_ROOT = RUN_DIR / "results" / "physics" / "reactions"
    DEFAULT_LOG_DIR = RUN_DIR / "logs"



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


def prepare_config_for_run(config_path: Path, output_folder: Path, snapshot_path: Path) -> Path:
    """Write a per-run snapshot of a config with an overridden output folder."""

    with open(config_path, "r", encoding="utf-8") as handle:
        cfg = json.load(handle)

    # Many configs use paths relative to the *config file* location.
    # Once we snapshot the config into the run directory, those relative paths
    # would break unless we rewrite them.
    base_dir = config_path.parent
    for key in ("species_database", "reaction_database"):
        value = cfg.get(key)
        if isinstance(value, str) and value and not Path(value).is_absolute():
            cfg[key] = str((base_dir / value).resolve())

    cfg.setdefault("output", {})["folder"] = str(output_folder)
    snapshot_path.parent.mkdir(parents=True, exist_ok=True)
    with open(snapshot_path, "w", encoding="utf-8") as handle:
        json.dump(cfg, handle, indent=2)
    output_folder.mkdir(parents=True, exist_ok=True)
    return snapshot_path


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
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    if result.returncode != 0:
        log("  ✗ Simulation failed", logf)

        stderr = (result.stderr or "").strip()
        stdout = (result.stdout or "").strip()

        def tail_lines(text: str, limit: int = 40) -> str:
            lines = text.splitlines()
            if len(lines) <= limit:
                return text
            return "\n".join(lines[-limit:])

        log(f"    STDERR: {tail_lines(stderr) if stderr else ''}", logf)
        if stdout:
            log("    STDOUT (tail):", logf)
            log(tail_lines(stdout), logf)
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
    parser = argparse.ArgumentParser(description="Run ICARION reaction kinetics validation scenarios")
    parser.add_argument(
        "--icarion-bin",
        type=Path,
        default=ICARION_BIN_DEFAULT,
        help="Path to icarion_main (default: build/src/icarion_main)",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print what would run without executing simulations",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_RESULTS_ROOT,
        help="Root folder for reaction outputs (per-scenario subfolders will be created)",
    )
    parser.add_argument(
        "--log-dir",
        type=Path,
        default=DEFAULT_LOG_DIR,
        help="Folder for validation logs",
    )
    parser.add_argument(
        "scenarios",
        nargs="*",
        help=f"Optional subset of scenarios to run (choices: {', '.join(SCENARIOS.keys())}); default: all",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()

    icarion_bin = args.icarion_bin
    if not icarion_bin.is_absolute():
        icarion_bin = (PROJECT_ROOT / icarion_bin).resolve()

    output_root = args.output_root
    if not output_root.is_absolute():
        output_root = (PROJECT_ROOT / output_root).resolve()

    log_dir = args.log_dir
    if not log_dir.is_absolute():
        log_dir = (PROJECT_ROOT / log_dir).resolve()
    log_dir.mkdir(parents=True, exist_ok=True)

    log_file = log_dir / "REACTION_KINETICS_VALIDATION.txt"

    requested = [s for s in args.scenarios if s]
    if requested:
        unknown = [s for s in requested if s not in SCENARIOS]
        if unknown:
            raise SystemExit(f"Unknown scenario(s): {', '.join(unknown)}")
        scenarios: Iterable[Scenario] = [SCENARIOS[s] for s in requested]
    else:
        scenarios = SCENARIOS.values()

    start = datetime.now().isoformat(timespec="microseconds")

    with open(log_file, "w", encoding="utf-8") as logf:
        log("=" * 80, logf)
        log("REACTION KINETICS VALIDATION", logf)
        log("=" * 80, logf)
        log(f"Project root: {PROJECT_ROOT}", logf)
        log(f"Icarion binary: {icarion_bin}", logf)
        log(f"Dry run: {args.dry_run}", logf)
        log(f"Output root: {output_root}", logf)
        log(f"Log dir: {log_dir}", logf)
        log(f"Start: {start}", logf)
        log("", logf)

        all_ok = True
        for scenario in scenarios:
            log(f"--- Scenario: {scenario.key} ---", logf)
            log(scenario.description, logf)

            for line in scenario_notes(scenario):
                log(f"  - {line}", logf)

            for run in scenario.runs:
                ensure_config_exists(run)

                # Write per-run config snapshot into the output tree so configs and outputs stay together.
                scenario_out = output_root / scenario.key
                cfg_snapshot = scenario_out / f"{run.config_path.stem}.run_config.json"
                cfg_to_run = prepare_config_for_run(run.config_path, scenario_out, cfg_snapshot)

                log(f"\n  → Running {run.config_path}", logf)
                ok = run_simulation(cfg_to_run, icarion_bin, args.dry_run, logf)
                if not ok:
                    all_ok = False
            log("", logf)

        end = datetime.now().isoformat(timespec="microseconds")
        if all_ok:
            log("All requested scenarios finished.", logf)
            log("Next step: use analyze_reactions.py for quantitative fits once time-resolved species data is available.", logf)
        else:
            log("One or more scenarios failed.", logf)
        log(f"End: {end}", logf)

    return 0 if all_ok else 1


if __name__ == "__main__":
    raise SystemExit(main())

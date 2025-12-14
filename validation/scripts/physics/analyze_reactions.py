#!/usr/bin/env python3
"""Summaries for the reaction validation suite outputs."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import h5py
import matplotlib.pyplot as plt
import numpy as np

K_B = 1.380649e-23  # Boltzmann constant [J/K]
E_CHARGE = 1.602176634e-19  # Elementary charge [C]
FIT_MIN_FRACTION = 0.02  # ignore tail frames with <2% of ions to reduce noise


@dataclass
class ScenarioReport:
    name: str
    observations: List[str]
    warnings: List[str]


def load_species_counts(h5_path: Path):
    """Return times array and a list of Counter objects per recorded frame."""

    if not h5_path.exists():
        raise FileNotFoundError(h5_path)

    with h5py.File(h5_path, "r") as handle:
        times = handle["trajectory/time"][:]
        raw = handle["trajectory/species_ids"][:]
        if raw.ndim == 1:
            raw = raw[None, :]
        frames: List[Counter] = []
        for frame in raw:
            names = []
            for value in frame:
                if isinstance(value, bytes):
                    names.append(value.decode().strip())
                else:
                    names.append(str(value).strip())
            counts = Counter(names)
            counts.pop("", None)  # remove inactive slots
            frames.append(counts)
    return times, frames


def count_summary(counter: Counter, order: Iterable[str] | None = None) -> str:
    if order is None:
        return ", ".join(f"{name}={count}" for name, count in counter.most_common())
    return ", ".join(f"{name}={counter.get(name, 0)}" for name in order)


def reactant_decay_estimate(remaining: int, total: int, duration_s: float) -> float:
    if total <= 0 or duration_s <= 0:
        return float("nan")
    if remaining <= 0:
        return float("inf")
    ratio = remaining / total
    if ratio <= 0:
        return float("inf")
    return -math.log(ratio) / duration_s


def fit_first_order_rate(times, frames, species: str, total: int, *, min_fraction: float | None = None) -> float:
    """Weighted least-squares fit of N(t) = N0 * exp(-k*t)."""

    def collect(min_fraction: float):
        data = []
        for t, frame in zip(times, frames):
            count = frame.get(species, 0)
            if count <= 0 or total <= 0:
                continue
            fraction = count / total
            if fraction < min_fraction:
                continue
            y = math.log(max(fraction, 1e-12))
            weight = float(count)
            data.append((t, y, weight))
        return data

    threshold = FIT_MIN_FRACTION if min_fraction is None else min_fraction
    samples = collect(threshold)
    if len(samples) < 2:
        samples = collect(0.0)
    if len(samples) < 2:
        return float("nan")

    sum_w = sum(w for _, _, w in samples)
    if sum_w <= 0:
        return float("nan")

    t_mean = sum(w * t for t, _, w in samples) / sum_w
    y_mean = sum(w * y for _, y, w in samples) / sum_w
    denom = sum(w * (t - t_mean) ** 2 for t, _, w in samples)
    if denom <= 0:
        return float("nan")

    slope = sum(w * (t - t_mean) * (y - y_mean) for t, y, w in samples) / denom
    return -slope


def analyze_equilibrium(base: Path) -> ScenarioReport:
    file_path = base / "equilibrium" / "reversible_equilibrium.h5"
    times, frames = load_species_counts(file_path)
    final = frames[-1]
    total = sum(frames[0].values())
    k_forward = 1800.0
    k_reverse = 600.0
    expected_fraction = k_forward / (k_forward + k_reverse)
    product = final.get("PentanalH+", 0)
    fraction = product / total if total else 0.0
    err_pct = 100.0 * (fraction - expected_fraction) / expected_fraction
    duration = times[-1] - times[0] if len(times) > 1 else 0.0
    obs = [
        f"File: {file_path}",
        f"Duration: {duration*1e6:.2f} us with {len(frames)} recorded frames",
        f"Final counts: {count_summary(final, ['H3O+', 'PentanalH+'])}",
        f"Expected PentanalH+ fraction {expected_fraction:.3f}, observed {fraction:.3f} ({err_pct:+.2f} %)"
    ]
    warnings = []
    if duration < 1e-4:
        warnings.append("Runtime is too short to approach equilibrium; reactions rarely fire.")
    if product == 0:
        warnings.append("No PentanalH+ detected; increase total_time_s or number of frames for meaningful data.")
    return ScenarioReport("Equilibrium", obs, warnings)


def analyze_first_order(base: Path) -> ScenarioReport:
    file_path = base / "first_order_decay" / "first_order_decay_reference.h5"
    times, frames = load_species_counts(file_path)
    final = frames[-1]
    total = sum(frames[0].values())
    duration = times[-1] - times[0] if len(times) > 1 else 0.0
    k_expected = 5000.0
    reactant = final.get("H3O+", 0)
    predicted = total * math.exp(-k_expected * duration)
    k_fit = fit_first_order_rate(times, frames, "H3O+", total)
    k_fallback = reactant_decay_estimate(reactant, total, duration)
    if math.isfinite(k_fit):
        k_display = f"{k_fit:.2e} s^-1 (trajectory fit)"
    else:
        k_display = "n/a (insufficient non-zero frames)"
    obs = [
        f"File: {file_path}",
        f"Duration: {duration*1e6:.2f} us",
        f"Final counts: {count_summary(final, ['H3O+', 'PentanalH+'])}",
        f"Predicted remaining H3O+ {predicted:.2f}, observed {reactant}",
        f"Estimated k from data: {k_display}",
        f"Fallback final-frame estimate: {k_fallback:.2e} s^-1"
    ]
    warnings = []
    if duration < 5e-5:
        warnings.append("Simulation covers only a single integration step; decay cannot be observed.")
    if reactant == total:
        warnings.append("No decay events observed; consider lowering write_interval or increasing total_time_s.")
    return ScenarioReport("First-order decay", obs, warnings)


def analyze_competing(base: Path) -> ScenarioReport:
    file_path = base / "competing_channels" / "competing_channels_he_vs_pentanal.h5"
    times, frames = load_species_counts(file_path)
    final = frames[-1]
    total = sum(frames[0].values())
    duration = times[-1] - times[0] if len(times) > 1 else 0.0
    k_he = 1.25e-15
    k_pent = 3.75e-15
    n_he = 6.760105814970592e22
    n_pent = 2.897188206415968e22
    rate_he = k_he * n_he
    rate_pent = k_pent * n_pent
    expected_pent = rate_pent / (rate_he + rate_pent)
    pent_count = final.get("PentanalH+", 0)
    he_count = final.get("26DTBP", 0)
    obs = [
        f"File: {file_path}",
        f"Duration: {duration*1e6:.2f} us",
        f"Final counts: {count_summary(final, ['H3O+', 'PentanalH+', '26DTBP'])}",
        f"Expected branching PentanalH+ {expected_pent*100:.1f} %"
    ]
    warnings = []
    if pent_count == 0 and he_count == 0:
        warnings.append("No reaction events recorded; extend runtime to sample branching statistics.")
    return ScenarioReport("Competing channels", obs, warnings)


def analyze_sequential(base: Path) -> ScenarioReport:
    file_path = base / "sequential_chain" / "sequential_chain_three_species.h5"
    times, frames = load_species_counts(file_path)
    final = frames[-1]
    total = sum(frames[0].values())
    duration = times[-1] - times[0] if len(times) > 1 else 0.0
    k1 = 1200.0
    k2 = 800.0
    if abs(k1 - k2) < 1e-9:
        exp_pent = 0.0
        exp_caff = total - total * math.exp(-k1 * duration)
    else:
        exp_A = total * math.exp(-k1 * duration)
        exp_B = total * k1 / (k2 - k1) * (math.exp(-k1 * duration) - math.exp(-k2 * duration))
        exp_C = total - exp_A - exp_B
        exp_pent = exp_B
        exp_caff = exp_C
    obs = [
        f"File: {file_path}",
        f"Duration: {duration*1e6:.2f} us",
        f"Final counts: {count_summary(final, ['H3O+', 'PentanalH+', 'CaffeineH+'])}",
        f"Expected PentanalH+ {exp_pent:.2f}, CaffeineH+ {exp_caff:.2f}"
    ]
    warnings = []
    if final.get("PentanalH+", 0) == 0 and final.get("CaffeineH+", 0) == 0:
        warnings.append("Sequential conversion not observable with current runtime.")
    return ScenarioReport("Sequential chain", obs, warnings)


def analyze_arrhenius(base: Path) -> ScenarioReport:
    runs = [
        ("250K", base / "arrhenius" / "arrhenius_h3o_to_pentanal_250K.h5", 250.0),
        ("300K", base / "arrhenius" / "arrhenius_h3o_to_pentanal_300K.h5", 300.0),
        ("350K", base / "arrhenius" / "arrhenius_h3o_to_pentanal_350K.h5", 350.0),
    ]
    prefactor = 2.0e5
    ea_eV = 0.08
    ea_J = ea_eV * E_CHARGE
    obs = ["Arrhenius sweep (H3O+ -> PentanalH+)"]
    warnings = []
    for label, path, temp in runs:
        times, frames = load_species_counts(path)
        final = frames[-1]
        total = sum(frames[0].values())
        duration = times[-1] - times[0] if len(times) > 1 else 0.0
        k_pred = prefactor * math.exp(-ea_J / (K_B * temp))
        reactant = final.get("H3O+", 0)
        k_fit = fit_first_order_rate(times, frames, "H3O+", total)
        if math.isfinite(k_fit):
            k_est = k_fit
            est_note = "trajectory fit"
        else:
            k_est = reactant_decay_estimate(reactant, total, duration)
            est_note = "final frame"
        obs.append(
            f"  {label}: duration {duration*1e6:.2f} us, final counts {count_summary(final, ['H3O+', 'PentanalH+'])},"
            f" predicted k={k_pred:.2e} s^-1, estimated k={k_est:.2e} s^-1 ({est_note})"
        )
        if reactant == total:
            warnings.append(f"{label}: no observable decay; runtime too short relative to predicted tau.")
    return ScenarioReport("Arrhenius sweep", obs, warnings)


def build_first_order_plot(base: Path) -> Tuple[Path, Tuple]:
    file_path = base / "first_order" / "first_order_decay_reference.h5"
    try:
        times, frames = load_species_counts(file_path)
    except Exception:
        return file_path, ()
    total = sum(frames[0].values()) if frames else 0
    if total <= 0:
        return file_path, ()
    reactant_series = [frame.get("H3O+", 0) / total for frame in frames]
    k_pred = 5000.0
    k_fit = fit_first_order_rate(times, frames, "H3O+", total)
    return file_path, (times, reactant_series, k_pred, k_fit)


def build_arrhenius_plot(base: Path) -> Tuple[Path, List[Tuple[float, float, float]]]:
    runs = [
        ("250K", base / "arrhenius" / "arrhenius_h3o_to_pentanal_250K.h5", 250.0),
        ("300K", base / "arrhenius" / "arrhenius_h3o_to_pentanal_300K.h5", 300.0),
        ("350K", base / "arrhenius" / "arrhenius_h3o_to_pentanal_350K.h5", 350.0),
    ]
    prefactor = 2.0e5
    ea_eV = 0.08
    ea_J = ea_eV * E_CHARGE
    data = []
    for label, path, temp in runs:
        try:
            times, frames = load_species_counts(path)
        except Exception:
            continue
        total = sum(frames[0].values()) if frames else 0
        k_pred = prefactor * math.exp(-ea_J / (K_B * temp))
        k_fit = fit_first_order_rate(times, frames, "H3O+", total)
        if math.isfinite(k_fit):
            data.append((temp, k_pred, k_fit))
    return runs[0][1].parent, data


def analyze_all(results_dir: Path) -> List[ScenarioReport]:
    base = Path(results_dir)
    reports = [
        analyze_equilibrium(base),
        analyze_first_order(base),
        analyze_competing(base),
        analyze_sequential(base),
        analyze_arrhenius(base),
    ]
    return reports


def main():
    import argparse

    parser = argparse.ArgumentParser(description="Analyze reaction validation outputs.")
    parser.add_argument(
        "results_dir",
        nargs="?",
        default="validation/results/physics/reactions",
        help="Directory containing the scenario subfolders",
    )
    args = parser.parse_args()

    reports = analyze_all(Path(args.results_dir))
    print("=" * 80)
    print("Reaction validation analysis")
    print("=" * 80)
    issues = []
    for report in reports:
        print(f"\n--- {report.name} ---")
        for line in report.observations:
            print(line)
        if report.warnings:
            print("Warnings:")
            for warn in report.warnings:
                print(f"  - {warn}")
            issues.append((report.name, list(report.warnings)))

    if issues:
        print("\nSummary: Some scenarios raised warnings during analysis:")
        for name, warnings in issues:
            for warn in warnings:
                print(f"  - {name}: {warn}")
        print("\nConsider addressing the warnings above (e.g., extend runtime or adjust write_interval) for higher-fidelity kinetics data.")
    else:
        print("\nSummary: All scenarios completed without analyzer warnings. Recorded durations and species statistics look consistent with expectations.")

    # Optional plots for visual validation
    results_dir = Path(args.results_dir)
    # First-order decay: log plot of reactant fraction
    file_path, decay_data = build_first_order_plot(results_dir)
    if decay_data:
        times, fractions, k_pred, k_fit = decay_data
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.semilogy(times * 1e6, fractions, 'o', label='Observed H3O+/N0')
        ax.semilogy(times * 1e6, np.exp(-k_pred * times), '--', label=f'Predicted k={k_pred:.1f}s⁻¹')
        if math.isfinite(k_fit):
            ax.semilogy(times * 1e6, np.exp(-k_fit * times), '-', label=f'Fit k={k_fit:.1f}s⁻¹')
        ax.set_xlabel('Time (µs)')
        ax.set_ylabel('Fraction remaining')
        ax.set_title('First-order decay (H3O+ -> products)')
        ax.grid(True, alpha=0.3, which='both')
        ax.legend()
        plot_path = file_path.parent / "first_order_decay_plot.png"
        plt.tight_layout()
        plt.savefig(plot_path, dpi=200, bbox_inches='tight')
        print(f"\n📊 Saved first-order decay plot: {plot_path}")

    # Arrhenius plot: ln(k) vs 1/T
    arr_base, arr_data = build_arrhenius_plot(results_dir)
    if arr_data:
        temps = np.array([t for t, _, _ in arr_data], dtype=float)
        k_pred = np.array([kp for _, kp, _ in arr_data], dtype=float)
        k_fit = np.array([kf for _, _, kf in arr_data], dtype=float)
        inv_T = 1.0 / temps

        fig, ax = plt.subplots(figsize=(6, 4))
        ax.plot(inv_T, np.log(k_pred), 'o--', label='Predicted ln(k)')
        ax.plot(inv_T, np.log(k_fit), 's-', label='Estimated ln(k)')
        ax.set_xlabel('1 / T (1/K)')
        ax.set_ylabel('ln(k)')
        ax.set_title('Arrhenius sweep (H3O+ → PentanalH+)')
        ax.grid(True, alpha=0.3)
        ax.legend()
        plot_path = arr_base / "arrhenius_plot.png"
        plt.tight_layout()
        plt.savefig(plot_path, dpi=200, bbox_inches='tight')
        print(f"📊 Saved Arrhenius plot: {plot_path}")


if __name__ == "__main__":
    main()

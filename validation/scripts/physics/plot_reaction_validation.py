#!/usr/bin/env python3
"""Generate summary figures for the reaction validation suite."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from analyze_reactions import K_B, E_CHARGE, fit_first_order_rate, load_species_counts

PROJECT_ROOT = Path(__file__).resolve().parents[3]
RESULTS_ROOT = PROJECT_ROOT / "validation" / "results" / "physics" / "reactions"
FIG_ROOT = PROJECT_ROOT / "validation" / "figures"
FIG_ROOT.mkdir(parents=True, exist_ok=True)
TOTAL_IONS = 10_000


def frames_to_series(frames, names):
    data = {}
    for name in names:
        data[name] = np.array([frame.get(name, 0) for frame in frames], dtype=float)
    return data


def plot_overview():
    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle("Reaction Validation Overview", fontsize=16)

    # Equilibrium
    eq_times, eq_frames = load_species_counts(RESULTS_ROOT / "equilibrium" / "reversible_equilibrium.h5")
    eq_series = frames_to_series(eq_frames, ["H3O+", "PentanalH+"])
    eq_ax = axes[0, 0]
    t_ms = eq_times * 1e3
    eq_ax.plot(t_ms, eq_series["H3O+"], label="H3O+")
    eq_ax.plot(t_ms, eq_series["PentanalH+"], label="PentanalH+")
    eq_ax.axhline(TOTAL_IONS * 0.75, color="tab:orange", linestyle="--", linewidth=1, label="75% analytic")
    eq_ax.set_title("Equilibrium: k_f=1800 s$^{-1}$, k_r=600 s$^{-1}$")
    eq_ax.set_xlabel("Time (ms)")
    eq_ax.set_ylabel("Ion count")
    eq_ax.set_xlim(0, t_ms[-1])
    eq_ax.legend(frameon=False)

    # First-order decay
    fo_times, fo_frames = load_species_counts(RESULTS_ROOT / "first_order_decay" / "first_order_decay_reference.h5")
    fo_series = frames_to_series(fo_frames, ["H3O+"])
    fo_ax = axes[0, 1]
    t_ms = fo_times * 1e3
    fo_ax.plot(t_ms, fo_series["H3O+"], label="H3O+")
    total = sum(fo_frames[0].values())
    k_fit = fit_first_order_rate(fo_times, fo_frames, "H3O+", total)
    if math.isfinite(k_fit):
        fit_curve = total * np.exp(-k_fit * fo_times)
        fo_ax.plot(t_ms, fit_curve, linestyle="--", color="tab:red", label=rf"Fit k = {k_fit/1e3:.2f}·10$^3$ s$^{{-1}}$")
    fo_ax.set_title("First-order decay: H$_3$O$^+$ → PentanalH$^+$")
    fo_ax.set_xlabel("Time (ms)")
    fo_ax.set_ylabel("Ion count")
    fo_ax.set_xlim(0, t_ms[-1])
    fo_ax.legend(frameon=False)

    # Sequential chain
    seq_times, seq_frames = load_species_counts(RESULTS_ROOT / "sequential_chain" / "sequential_chain_three_species.h5")
    seq_series = frames_to_series(seq_frames, ["H3O+", "PentanalH+", "CaffeineH+"])
    seq_ax = axes[1, 0]
    t_ms = seq_times * 1e3
    for name, style in zip(["H3O+", "PentanalH+", "CaffeineH+"], ["-", "--", "-"]):
        seq_ax.plot(t_ms, seq_series[name], label=name, linestyle=style)
    seq_ax.set_title("Sequential chain: H$_3$O$^+$ → PentanalH$^+$ → CaffeineH$^+$")
    seq_ax.set_xlabel("Time (ms)")
    seq_ax.set_ylabel("Ion count")
    seq_ax.set_xlim(0, t_ms[-1])
    seq_ax.legend(frameon=False)

    # Arrhenius sweep
    arr_ax = axes[1, 1]
    temps = np.array([250.0, 300.0, 350.0])
    prefactor = 2.0e5
    ea_j = 0.08 * E_CHARGE
    k_pred = prefactor * np.exp(-ea_j / (K_B * temps))
    labels = ["250K", "300K", "350K"]
    k_fit = []
    for lbl in ["250K", "300K", "350K"]:
        file_path = RESULTS_ROOT / "arrhenius" / f"arrhenius_h3o_to_pentanal_{lbl}.h5"
        times, frames = load_species_counts(file_path)
        total = sum(frames[0].values())
        k_fit.append(fit_first_order_rate(times, frames, "H3O+", total))
    k_fit = np.array(k_fit)
    arr_ax.plot(temps, k_pred, label="Theory", color="tab:blue")
    arr_ax.scatter(temps, k_fit, label="Simulation fit", color="tab:orange", zorder=5)
    arr_ax.set_title("Arrhenius sweep: H$_3$O$^+$ → PentanalH$^+$")
    arr_ax.set_xlabel("Temperature (K)")
    arr_ax.set_ylabel("Rate constant k (s$^{-1}$)")
    arr_ax.grid(True, linestyle=":", alpha=0.5)
    arr_ax.legend(frameon=False)

    fig.tight_layout(rect=[0, 0.01, 1, 0.97])
    out_path = FIG_ROOT / "reaction_validation_overview.png"
    fig.savefig(out_path, dpi=200)
    print(f"Saved {out_path}")


def plot_competing():
    times, frames = load_species_counts(RESULTS_ROOT / "competing_channels" / "competing_channels_he_vs_pentanal.h5")
    final = frames[-1]
    total = sum(frames[0].values())
    actual = np.array([final.get("PentanalH+", 0), final.get("26DTBP", 0)], dtype=float)
    expected_fraction = np.array([0.562, 0.438])
    expected = total * expected_fraction

    fig, ax = plt.subplots(figsize=(6, 4))
    x = np.arange(2)
    width = 0.35
    ax.bar(x - width / 2, expected, width, color="lightgray", label="Expected")
    ax.bar(x + width / 2, actual, width, color="tab:green", label="Simulation")
    ax.set_xticks(x, ["PentanalH+", "2,6-DTBP"])
    ax.set_ylabel("Ion count")
    ax.set_title("Competing bimolecular channels")
    ax.legend(frameon=False)
    ax.set_ylim(0, total)
    fig.tight_layout()
    out_path = FIG_ROOT / "reaction_competing_branching.png"
    fig.savefig(out_path, dpi=200)
    print(f"Saved {out_path}")


def main():
    plot_overview()
    plot_competing()


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Analyze TOF arrival spectra or trap frequency spectra from trajectory output."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

if __package__ is None or __package__ == "":
    import sys

    sys.path.append(str(Path(__file__).resolve().parent.parent))

from analysis.physics import tof_mass_to_charge_da_per_e
from analysis.spectra import fft_frequency_spectrum, histogram_spectrum, resample_to_uniform_grid, time_grid_stats, top_peaks
from analysis.trajectory import (
    ensure_nonempty_selection,
    find_arrival_times,
    normalize_species_filter,
    open_trajectory,
    read_domain_geometry,
    read_positions_subset,
    select_ions_from_trajectory,
    species_for_indices,
)


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Create TOF m/z spectra or FFT frequency spectra from ICARION trajectories.")
    p.add_argument("--traj", "--input", dest="traj", type=Path, required=True, help="Trajectory HDF5 file.")
    p.add_argument("--mode", choices=("tof", "frequency"), required=True, help="Analysis mode.")
    p.add_argument("--out", "--output", dest="out", type=Path, default=Path("analysis/output/spectrum.png"), help="Output plot path.")
    p.add_argument("--out-csv", type=Path, default=Path("analysis/output/spectrum.csv"), help="Output CSV path.")
    p.add_argument("--species", nargs="+", default=None, help="Optional species filter.")
    p.add_argument("--domain-index", type=int, default=0, help="Domain index for geometry defaults.")
    p.add_argument("--max-ions", type=int, default=1000, help="Maximum ions to analyze.")
    p.add_argument("--rng-seed", type=int, default=0, help="Subsampling seed.")
    p.add_argument("--bins", type=int, default=100, help="Histogram bins for TOF mode.")
    p.add_argument("--arrival-z-m", type=float, default=None, help="Detector/arrival plane z [m].")
    p.add_argument("--tol-frac", type=float, default=0.02, help="Arrival tolerance fraction.")
    p.add_argument("--acceleration-voltage", type=float, default=None, help="TOF acceleration voltage [V].")
    p.add_argument("--flight-distance", type=float, default=None, help="TOF flight distance [m].")
    p.add_argument("--coordinate", choices=("x", "y", "z", "r"), default="x", help="Signal coordinate for frequency mode.")
    p.add_argument("--min-frequency", type=float, default=0.0, help="Minimum plotted FFT frequency [Hz].")
    p.add_argument("--top-peaks", type=int, default=5, help="Number of peak rows to write.")
    p.add_argument(
        "--resample",
        choices=("none", "uniform"),
        default="none",
        help="Frequency mode only: resample nonuniform trajectories onto a uniform grid before FFT.",
    )
    p.add_argument("--resample-dt", type=float, default=None, help="Uniform resampling step [s]; default uses median input dt.")
    p.add_argument("--max-dt-jitter-frac", type=float, default=1e-6, help="Allowed relative dt jitter before nonuniform input is rejected.")
    p.add_argument("--out-meta", type=Path, default=None, help="Optional JSON metadata output path.")
    return p.parse_args()


def _selected_positions(args: argparse.Namespace):
    with open_trajectory(args.traj) as h5:
        traj = h5["trajectory"]
        species_ids, ion_indices = select_ions_from_trajectory(
            traj,
            species_filter=normalize_species_filter(args.species),
            max_ions=args.max_ions,
            max_per_species=None,
            rng_seed=args.rng_seed,
        )
        ensure_nonempty_selection(ion_indices)
        time, positions = read_positions_subset(traj, ion_indices, time_stride=1, max_frames=None)
        length_m, radius_m, z_origin_m = read_domain_geometry(h5, args.domain_index, traj, ion_indices)
    return np.asarray(time), np.asarray(positions), species_for_indices(species_ids, ion_indices), length_m, radius_m, z_origin_m


def _write_rows(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def _tof(args: argparse.Namespace) -> None:
    time, positions, species, length_m, radius_m, z_origin_m = _selected_positions(args)
    arrivals, valid_mask, _ = find_arrival_times(
        times=time,
        positions=positions,
        z_origin_m=z_origin_m,
        length_m=length_m,
        radius_m=radius_m,
        tol_frac=args.tol_frac,
        arrival_plane_z_m=args.arrival_z_m,
    )
    if arrivals.size == 0:
        raise RuntimeError("No ions reached the detector plane.")

    flight_distance = args.flight_distance if args.flight_distance is not None else length_m
    x_label = "flight_time_s"
    values = arrivals
    if args.acceleration_voltage is not None:
        values = tof_mass_to_charge_da_per_e(arrivals, args.acceleration_voltage, flight_distance)
        x_label = "mz_Da_per_e"

    centers, counts = histogram_spectrum(values, bins=args.bins)
    species_valid = species[valid_mask]
    rows = [
        {"kind": "ion", "species": sp, "flight_time_s": t, "mz_Da_per_e": mz if x_label == "mz_Da_per_e" else ""}
        for sp, t, mz in zip(species_valid, arrivals, values)
    ]
    for peak in top_peaks(centers, counts, args.top_peaks):
        rows.append(
            {
                "kind": "peak",
                "species": "",
                "flight_time_s": peak.x if x_label == "flight_time_s" else "",
                "mz_Da_per_e": peak.x if x_label == "mz_Da_per_e" else "",
            }
        )

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(centers, counts, lw=1.8)
    ax.set_xlabel("m/z [Da/e]" if x_label == "mz_Da_per_e" else "flight time [s]")
    ax.set_ylabel("counts")
    ax.grid(alpha=0.25)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=180)
    plt.close(fig)
    _write_rows(args.out_csv, rows, ["kind", "species", "flight_time_s", "mz_Da_per_e"])


def _frequency(args: argparse.Namespace) -> None:
    time, positions, _species, _length_m, _radius_m, _z_origin_m = _selected_positions(args)
    if args.coordinate == "r":
        signal_by_ion = np.sqrt(positions[:, :, 0] ** 2 + positions[:, :, 1] ** 2)
    else:
        axis = {"x": 0, "y": 1, "z": 2}[args.coordinate]
        signal_by_ion = positions[:, :, axis]
    signal = np.nanmean(signal_by_ion, axis=1)
    stats = time_grid_stats(time, max_jitter_frac=args.max_dt_jitter_frac)
    fft_time = time
    fft_signal = signal
    resampled = False
    resample_dt = None
    if not stats.is_uniform:
        if args.resample != "uniform":
            raise ValueError(
                "Nonuniform time grid for FFT: "
                f"dt_min={stats.dt_min_s:.6g}s, dt_median={stats.dt_median_s:.6g}s, "
                f"dt_max={stats.dt_max_s:.6g}s, jitter={stats.dt_jitter_frac:.6g}; "
                "use --resample uniform"
            )
        resampled_signal = resample_to_uniform_grid(time, signal, dt_s=args.resample_dt)
        fft_time = resampled_signal.times_s
        fft_signal = resampled_signal.signal
        resampled = True
        resample_dt = resampled_signal.dt_s

    freq, amp = fft_frequency_spectrum(
        fft_time,
        fft_signal,
        min_frequency_hz=args.min_frequency,
        max_jitter_frac=args.max_dt_jitter_frac,
    )
    peaks = top_peaks(freq, amp, args.top_peaks, min_separation_bins=2)

    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(freq, amp, lw=1.5)
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("amplitude")
    ax.grid(alpha=0.25)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(args.out, dpi=180)
    plt.close(fig)
    rows = [{"kind": "peak", "frequency_Hz": p.x, "amplitude": p.intensity} for p in peaks]
    _write_rows(args.out_csv, rows, ["kind", "frequency_Hz", "amplitude"])
    _write_frequency_metadata(
        args,
        stats=stats,
        resampled=resampled,
        resample_dt=resample_dt,
        fft_time=fft_time,
        frequency=freq,
    )


def _write_frequency_metadata(
    args: argparse.Namespace,
    stats,
    resampled: bool,
    resample_dt: float | None,
    fft_time: np.ndarray,
    frequency: np.ndarray,
) -> None:
    meta_path = args.out_meta or args.out_csv.with_suffix(args.out_csv.suffix + ".meta.json")
    dt = float(np.median(np.diff(fft_time))) if fft_time.size > 1 else float("nan")
    payload = {
        "source_file": str(args.traj),
        "mode": "frequency",
        "coordinate": args.coordinate,
        "resampled": resampled,
        "resample_dt_s": resample_dt,
        "max_dt_jitter_frac": args.max_dt_jitter_frac,
        "input_time_grid": stats.__dict__,
        "fft": {
            "n_samples": int(fft_time.size),
            "dt_s": dt,
            "frequency_resolution_Hz": float(frequency[1] - frequency[0]) if frequency.size > 1 else None,
            "nyquist_Hz": float(0.5 / dt) if dt > 0.0 else None,
        },
        "warnings": [] if stats.is_uniform or resampled else ["nonuniform_time_grid"],
    }
    meta_path.parent.mkdir(parents=True, exist_ok=True)
    meta_path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def main() -> int:
    args = parse_args()
    try:
        if args.mode == "tof":
            _tof(args)
        else:
            _frequency(args)
    except ValueError as exc:
        raise SystemExit(f"error: {exc}") from None
    print(f"Wrote {args.out}")
    print(f"Wrote {args.out_csv}")
    if args.mode == "frequency":
        print(f"Wrote {args.out_meta or args.out_csv.with_suffix(args.out_csv.suffix + '.meta.json')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

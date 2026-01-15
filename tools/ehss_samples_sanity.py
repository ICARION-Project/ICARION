#!/usr/bin/env python3
import argparse
import json
import math
import statistics
from pathlib import Path


GAS_RADII_M = {
    "He": 1.3e-10,
    "Ar": 1.70e-10,
    "CO2": 1.65e-10,
    "Ne": 1.38e-10,
    "N2": 1.82e-10,
    "O2": 1.73e-10,
    "H2O": 1.58e-10,
}
ANGSTROM_TO_M = 1.0e-10


def summarize(values, bins):
    if not values:
        return None
    mean = statistics.mean(values)
    stdev = statistics.pstdev(values) if len(values) > 1 else 0.0
    vmin = min(values)
    vmax = max(values)

    if bins <= 0:
        bins = 10
    if vmax <= vmin:
        edges = [vmin, vmax]
        counts = [len(values)]
    else:
        step = (vmax - vmin) / bins
        edges = [vmin + i * step for i in range(bins + 1)]
        counts = [0] * bins
        for v in values:
            idx = int((v - vmin) / step)
            if idx >= bins:
                idx = bins - 1
            counts[idx] += 1

    return mean, stdev, vmin, vmax, edges, counts


def print_summary(name, values, bins, oapa_mean=None):
    summary = summarize(values, bins)
    if summary is None:
        print(f"{name}: no values")
        return
    mean, stdev, vmin, vmax, edges, counts = summary
    print(f"{name}:")
    print(f"  count: {len(values)}")
    print(f"  mean: {mean:.6e} m^2")
    print(f"  std:  {stdev:.6e} m^2")
    print(f"  min:  {vmin:.6e} m^2")
    print(f"  max:  {vmax:.6e} m^2")
    if oapa_mean is not None:
        diff = mean - oapa_mean
        rel = diff / oapa_mean if oapa_mean != 0 else float("nan")
        print(f"  oapa_mean: {oapa_mean:.6e} m^2")
        print(f"  delta:     {diff:.6e} m^2 ({rel:+.2%})")
    print("  histogram:")
    for i, c in enumerate(counts):
        lo = edges[i]
        hi = edges[i + 1]
        print(f"    [{lo:.6e}, {hi:.6e}): {c}")


def load_geometry(geom_path):
    data = json.loads(Path(geom_path).read_text())
    if "molecule" not in data:
        raise ValueError("geometry file missing 'molecule' object")
    mol = data["molecule"]
    atoms = mol.get("atoms", [])
    if not atoms:
        raise ValueError("geometry file has no atoms")

    centers = []
    radii = []
    for atom in atoms:
        pos = atom.get("pos")
        if not isinstance(pos, list) or len(pos) != 3:
            raise ValueError("atom position must be a 3-element array")
        centers.append([pos[0] * ANGSTROM_TO_M, pos[1] * ANGSTROM_TO_M, pos[2] * ANGSTROM_TO_M])

        sigma_a = atom.get("LJ_sigma_angstrom", 0.0)
        radii.append(0.5 * sigma_a * ANGSTROM_TO_M)

    return centers, radii


def quat_to_rotation(q):
    w, x, y, z = q
    norm = math.sqrt(w * w + x * x + y * y + z * z)
    if norm <= 0.0:
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    w /= norm
    x /= norm
    y /= norm
    z /= norm

    return [
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - z * w), 2.0 * (x * z + y * w)],
        [2.0 * (x * y + z * w), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - x * w)],
        [2.0 * (x * z - y * w), 2.0 * (y * z + x * w), 1.0 - 2.0 * (x * x + y * y)],
    ]


def rotate_vec(v, R):
    return [
        R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
        R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
        R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2],
    ]


def compute_oapa_mean(orientations, centers, radii, gas_radius):
    if gas_radius is None:
        return None
    areas = []
    for q in orientations:
        R = quat_to_rotation(q)
        rmax = 0.0
        for c, r in zip(centers, radii):
            rc = rotate_vec(c, R)
            r_xy = math.sqrt(rc[0] * rc[0] + rc[1] * rc[1]) + r + gas_radius
            if r_xy > rmax:
                rmax = r_xy
        areas.append(math.pi * rmax * rmax)
    if not areas:
        return None
    return statistics.mean(areas)


def plot_histogram(values, bins, title, out_path):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception as exc:
        print(f"Plot skipped (matplotlib unavailable): {exc}")
        return

    plt.figure(figsize=(6, 4))
    plt.hist(values, bins=bins, edgecolor="black")
    plt.title(title)
    plt.xlabel("Area (m^2)")
    plt.ylabel("Count")
    plt.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(out_path, dpi=150)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Sanity checks for EHSS samples JSON")
    parser.add_argument("--input", required=True, help="Path to EHSS samples JSON")
    parser.add_argument("--gas", default="", help="Gas key to inspect (default: all)")
    parser.add_argument("--bins", type=int, default=20, help="Histogram bins (default: 20)")
    parser.add_argument("--geometry-file", default="", help="Geometry JSON for OAPA comparison")
    parser.add_argument("--plot", default="", help="Write histogram plot to this PNG (requires --gas)")
    parser.add_argument("--plot-dir", default="", help="Write one histogram per gas into this directory")
    args = parser.parse_args()

    path = Path(args.input)
    with path.open("r", encoding="utf-8") as f:
        data = json.load(f)

    areas = data.get("areas_by_gas_m2", {})
    if not isinstance(areas, dict) or not areas:
        raise SystemExit("areas_by_gas_m2 is missing or empty")

    orientations = data.get("orientations_quat", [])
    centers = radii = None
    if args.geometry_file:
        centers, radii = load_geometry(args.geometry_file)

    if args.gas:
        if args.gas not in areas:
            raise SystemExit(f"gas '{args.gas}' not found in areas_by_gas_m2")
        oapa_mean = None
        if centers is not None and orientations:
            gas_radius = GAS_RADII_M.get(args.gas)
            oapa_mean = compute_oapa_mean(orientations, centers, radii, gas_radius)
        print_summary(args.gas, areas[args.gas], args.bins, oapa_mean)
        if args.plot:
            plot_histogram(areas[args.gas], args.bins, args.gas, Path(args.plot))
    else:
        for gas, values in areas.items():
            oapa_mean = None
            if centers is not None and orientations:
                gas_radius = GAS_RADII_M.get(gas)
                oapa_mean = compute_oapa_mean(orientations, centers, radii, gas_radius)
            print_summary(gas, values, args.bins, oapa_mean)
            if args.plot_dir:
                out_path = Path(args.plot_dir) / f"ehss_hist_{gas}.png"
                plot_histogram(values, args.bins, gas, out_path)


if __name__ == "__main__":
    main()

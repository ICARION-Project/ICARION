#!/usr/bin/env python3
# Plot ICARION thread scaling from summary CSV

from __future__ import annotations

import argparse
import csv
import pathlib
from collections import defaultdict


def load_summary(path: pathlib.Path):
    rows = []
    with path.open("r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            rows.append(row)
    return rows


def to_float(value: str) -> float | None:
    if value is None or value == "":
        return None
    return float(value)


def to_int(value: str) -> int:
    return int(value)


def main() -> int:
    parser = argparse.ArgumentParser(description="Plot ICARION thread scaling")
    parser.add_argument("--summary", required=True, help="Summary CSV from benchmark script")
    parser.add_argument("--metric", choices=["speedup", "efficiency"], default="speedup",
                        help="Plot metric (speedup or efficiency)")
    parser.add_argument("--output", default=None, help="Output image path (png)")
    parser.add_argument("--pdf", action="store_true", help="Also write PDF")
    parser.add_argument("--title", default=None, help="Plot title override")
    args = parser.parse_args()

    summary_path = pathlib.Path(args.summary)
    rows = load_summary(summary_path)
    if not rows:
        raise SystemExit("Summary CSV is empty")

    grouped = defaultdict(list)
    for row in rows:
        ions = to_int(row["ions"])
        threads = to_int(row["threads"])
        speedup = to_float(row["speedup"])
        efficiency = to_float(row["efficiency_percent"])
        grouped[ions].append({
            "threads": threads,
            "speedup": speedup,
            "efficiency": efficiency,
        })

    import matplotlib.pyplot as plt

    plt.figure(figsize=(6.8, 4.2))
    max_threads = 0
    for ions, series in sorted(grouped.items()):
        series = sorted(series, key=lambda r: r["threads"])
        threads = [r["threads"] for r in series]
        max_threads = max(max_threads, max(threads))
        if args.metric == "speedup":
            values = [r["speedup"] for r in series]
            label = f"N={ions:,}"
        else:
            values = [r["efficiency"] for r in series]
            label = f"N={ions:,}"
        plt.plot(threads, values, marker="o", label=label)

    if args.metric == "speedup":
        ideal = list(range(1, max_threads + 1))
        plt.plot(ideal, ideal, linestyle="--", color="gray", linewidth=1, label="Ideal")

    plt.xlabel("Threads")
    plt.ylabel("Speedup" if args.metric == "speedup" else "Efficiency (%)")
    if args.title:
        plt.title(args.title)
    else:
        plt.title("ICARION Thread Scaling")
    plt.grid(True, linestyle=":", linewidth=0.7, alpha=0.7)
    plt.legend()

    if args.output:
        out_path = pathlib.Path(args.output)
    else:
        out_path = summary_path.parent / f"thread_scaling_{args.metric}.png"
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    if args.pdf:
        pdf_path = out_path.with_suffix(".pdf")
        plt.savefig(pdf_path)

    print(f"Wrote: {out_path}")
    if args.pdf:
        print(f"Wrote: {pdf_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

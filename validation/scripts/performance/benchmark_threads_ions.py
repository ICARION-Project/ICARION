#!/usr/bin/env python3
# ICARION thread/ion scaling benchmark (CPU-only)

from __future__ import annotations

import argparse
import csv
import json
import os
import pathlib
import platform
import re
import statistics
import subprocess
import sys
import time
from typing import List, Tuple


CPU_TIME_RE = re.compile(r"CPU time:\s+([0-9.]+)\s*s")


def parse_int_list(value: str) -> List[int]:
    items = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        items.append(int(part))
    return items


def load_config(path: pathlib.Path) -> dict:
    with path.open("r", encoding="utf-8") as f:
        return json.load(f)


def scaled_counts(counts: List[int], target_total: int) -> List[int]:
    total = sum(counts)
    if total <= 0:
        raise ValueError("Sum of ion counts must be > 0")
    if target_total <= 0:
        raise ValueError("Target ion count must be > 0")

    scale = target_total / total
    raw = [c * scale for c in counts]
    base = [int(x) for x in raw]
    remainder = target_total - sum(base)
    if remainder > 0:
        frac = sorted(
            range(len(raw)),
            key=lambda i: raw[i] - base[i],
            reverse=True,
        )
        for idx in frac[:remainder]:
            base[idx] += 1
    return base


def resolve_path(value: str, base_dir: pathlib.Path) -> str:
    if not value:
        return value
    path = pathlib.Path(value)
    if path.is_absolute():
        return value
    return str((base_dir / path).resolve())


def normalize_database_paths(cfg: dict, base_cfg_path: pathlib.Path) -> None:
    base_dir = base_cfg_path.parent
    if "species_database" in cfg:
        cfg["species_database"] = resolve_path(cfg["species_database"], base_dir)
    if "reaction_database" in cfg:
        cfg["reaction_database"] = resolve_path(cfg["reaction_database"], base_dir)

    database = cfg.get("database")
    if isinstance(database, dict):
        if "species" in database:
            database["species"] = resolve_path(database["species"], base_dir)
        if "reactions" in database:
            database["reactions"] = resolve_path(database["reactions"], base_dir)


def write_scaled_config(base_cfg: dict, target_total: int, out_path: pathlib.Path,
                        base_cfg_path: pathlib.Path) -> None:
    cfg = json.loads(json.dumps(base_cfg))
    normalize_database_paths(cfg, base_cfg_path)
    species = cfg.get("ions", {}).get("species", [])
    if not species:
        raise ValueError("Config has no ions.species entries")

    counts = [int(s.get("count", 0)) for s in species]
    new_counts = scaled_counts(counts, target_total)
    for s, new_count in zip(species, new_counts):
        s["count"] = int(new_count)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8") as f:
        json.dump(cfg, f, indent=2)
        f.write("\n")


def run_cmd(cmd: List[str]) -> Tuple[float | None, float, int, str]:
    start = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    end = time.perf_counter()
    output = proc.stdout
    cpu_time = None
    match = CPU_TIME_RE.search(output)
    if match:
        cpu_time = float(match.group(1))
    return cpu_time, end - start, proc.returncode, output


def run_cmd_env(cmd: List[str], env: dict) -> Tuple[float | None, float, int, str]:
    start = time.perf_counter()
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    end = time.perf_counter()
    output = proc.stdout
    cpu_time = None
    match = CPU_TIME_RE.search(output)
    if match:
        cpu_time = float(match.group(1))
    return cpu_time, end - start, proc.returncode, output


def find_repo_root(start: pathlib.Path) -> pathlib.Path | None:
    cur = start.resolve()
    if cur.is_file():
        cur = cur.parent
    for parent in [cur] + list(cur.parents):
        if (parent / ".git").exists():
            return parent
    return None


def run_git(args: List[str], cwd: pathlib.Path) -> str | None:
    try:
        proc = subprocess.run(
            ["git"] + args, cwd=str(cwd),
            stdout=subprocess.PIPE, stderr=subprocess.DEVNULL, text=True, check=False
        )
    except FileNotFoundError:
        return None
    if proc.returncode != 0:
        return None
    return proc.stdout.strip()


def read_cpu_model() -> str | None:
    cpuinfo = pathlib.Path("/proc/cpuinfo")
    if not cpuinfo.exists():
        return None
    with cpuinfo.open("r", encoding="utf-8") as f:
        for line in f:
            if line.lower().startswith("model name"):
                return line.split(":", 1)[1].strip()
    return None


def read_build_type(repo_root: pathlib.Path | None) -> str | None:
    if not repo_root:
        return None
    cache_path = repo_root / "build" / "CMakeCache.txt"
    if not cache_path.exists():
        return None
    with cache_path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith("CMAKE_BUILD_TYPE:STRING="):
                return line.split("=", 1)[1].strip()
    return None


def collect_metadata(config_path: pathlib.Path, args: argparse.Namespace) -> dict:
    repo_root = find_repo_root(config_path)
    metadata = {
        "timestamp": time.strftime("%Y-%m-%d %H:%M:%S"),
        "hostname": platform.node(),
        "platform": platform.platform(),
        "python": platform.python_version(),
        "cpu_count": os.cpu_count(),
        "cpu_model": read_cpu_model(),
        "repo_root": str(repo_root) if repo_root else None,
        "git_commit": run_git(["rev-parse", "HEAD"], repo_root) if repo_root else None,
        "git_branch": run_git(["rev-parse", "--abbrev-ref", "HEAD"], repo_root) if repo_root else None,
        "build_type": read_build_type(repo_root),
        "bin": str(args.bin),
        "config": str(args.config),
        "threads": args.threads,
        "ions": args.ions,
        "repeats": args.repeats,
        "warmup": args.warmup,
        "write_interval": args.write_interval,
        "omp_proc_bind": args.omp_proc_bind,
        "omp_places": args.omp_places,
        "omp_dynamic": args.omp_dynamic,
        "omp_schedule": args.omp_schedule,
    }
    return metadata


def main() -> int:
    parser = argparse.ArgumentParser(description="ICARION thread/ion scaling benchmark")
    parser.add_argument("--bin", default="./build/src/icarion_main", help="Path to icarion_main")
    parser.add_argument("--config", default="examples/ims/ims_basic.json", help="Base config JSON")
    parser.add_argument("--threads", default="1,2,4,8", help="Comma-separated thread counts")
    parser.add_argument("--ions", default="5000,10000,50000,100000", help="Comma-separated ion totals")
    parser.add_argument("--repeats", type=int, default=3, help="Repeats per configuration")
    parser.add_argument("--warmup", type=int, default=1, help="Warmup runs per configuration")
    parser.add_argument("--output-dir", default="results/bench_threads_ions", help="Output directory")
    parser.add_argument("--write-interval", type=int, default=None, help="Override simulation.write_interval")
    parser.add_argument("--keep-output", action="store_true", help="Keep run output (stdout)")
    parser.add_argument("--omp-proc-bind", default=None, help="Set OMP_PROC_BIND")
    parser.add_argument("--omp-places", default=None, help="Set OMP_PLACES")
    parser.add_argument("--omp-dynamic", default="false", help="Set OMP_DYNAMIC")
    parser.add_argument("--omp-schedule", default=None, help="Set OMP_SCHEDULE")
    args = parser.parse_args()

    bin_path = pathlib.Path(args.bin)
    config_path = pathlib.Path(args.config)
    output_dir = pathlib.Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    threads_list = parse_int_list(args.threads)
    ions_list = parse_int_list(args.ions)
    if not threads_list or not ions_list:
        raise ValueError("threads and ions lists must be non-empty")

    base_cfg = load_config(config_path)
    config_dir = output_dir / "configs"
    configs = {}
    for ions in ions_list:
        cfg_path = config_dir / f"{config_path.stem}_N{ions}.json"
        write_scaled_config(base_cfg, ions, cfg_path, config_path)
        configs[ions] = cfg_path

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    csv_path = output_dir / f"benchmark_{timestamp}.csv"
    summary_path = output_dir / f"summary_{timestamp}.csv"
    metadata_path = output_dir / f"metadata_{timestamp}.json"
    metadata = collect_metadata(config_path, args)
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    results = []
    with csv_path.open("w", encoding="utf-8") as csv_file:
        csv_file.write("ions,threads,repeat,cpu_time_s,wall_time_s,returncode\n")
        for ions in ions_list:
            cfg_path = configs[ions]
            for threads in threads_list:
                base_env = os.environ.copy()
                base_env["OMP_NUM_THREADS"] = str(threads)
                if args.omp_proc_bind:
                    base_env["OMP_PROC_BIND"] = args.omp_proc_bind
                if args.omp_places:
                    base_env["OMP_PLACES"] = args.omp_places
                if args.omp_dynamic is not None:
                    base_env["OMP_DYNAMIC"] = args.omp_dynamic
                if args.omp_schedule:
                    base_env["OMP_SCHEDULE"] = args.omp_schedule

                for warm in range(1, args.warmup + 1):
                    output_file = f"{config_path.stem}_N{ions}_T{threads}_warm{warm}.h5"
                    cmd = [
                        str(bin_path),
                        "--threads", str(threads),
                        "--output-dir", str(output_dir),
                        "--output", output_file,
                        "--set", "output.print_progress=false",
                    ]
                    if args.write_interval is not None:
                        cmd.extend(["--set", f"simulation.write_interval={args.write_interval}"])
                    cmd.append(str(cfg_path))

                    cpu_time, wall_time, rc, output = run_cmd_env(cmd, base_env)
                    if rc != 0:
                        log_dir = output_dir / "logs"
                        log_dir.mkdir(parents=True, exist_ok=True)
                        log_path = log_dir / f"{config_path.stem}_N{ions}_T{threads}_warm{warm}.log"
                        log_path.write_text(output, encoding="utf-8")
                        sys.stderr.write(
                            f"Warmup failed (ions={ions}, threads={threads}): rc={rc}\n"
                        )
                        return 1

                for rep in range(1, args.repeats + 1):
                    output_file = f"{config_path.stem}_N{ions}_T{threads}_r{rep}.h5"
                    cmd = [
                        str(bin_path),
                        "--threads", str(threads),
                        "--output-dir", str(output_dir),
                        "--output", output_file,
                        "--set", "output.print_progress=false",
                    ]
                    if args.write_interval is not None:
                        cmd.extend(["--set", f"simulation.write_interval={args.write_interval}"])
                    cmd.append(str(cfg_path))

                    cpu_time, wall_time, rc, output = run_cmd_env(cmd, base_env)
                    csv_file.write(
                        f"{ions},{threads},{rep},"
                        f"{'' if cpu_time is None else cpu_time},"
                        f"{wall_time:.3f},{rc}\n"
                    )
                    csv_file.flush()
                    results.append({
                        "ions": ions,
                        "threads": threads,
                        "repeat": rep,
                        "cpu_time_s": cpu_time,
                        "wall_time_s": wall_time,
                        "returncode": rc,
                    })
                    if rc != 0:
                        log_dir = output_dir / "logs"
                        log_dir.mkdir(parents=True, exist_ok=True)
                        log_path = log_dir / f"{config_path.stem}_N{ions}_T{threads}_r{rep}.log"
                        log_path.write_text(output, encoding="utf-8")
                    if args.keep_output:
                        sys.stdout.write(output)

    summary_rows = []
    grouped = {}
    for row in results:
        key = (row["ions"], row["threads"])
        grouped.setdefault(key, []).append(row)

    baseline_map = {}
    for ions in ions_list:
        base_key = (ions, 1)
        if base_key not in grouped:
            continue
        ok = [r for r in grouped[base_key] if r["returncode"] == 0]
        if not ok:
            continue
        wall = [r["wall_time_s"] for r in ok]
        baseline_map[ions] = statistics.median(wall)

    for ions in ions_list:
        for threads in threads_list:
            rows = grouped.get((ions, threads), [])
            ok = [r for r in rows if r["returncode"] == 0]
            fail_count = len(rows) - len(ok)
            if ok:
                wall = [r["wall_time_s"] for r in ok]
                wall_median = statistics.median(wall)
                wall_min = min(wall)
                wall_max = max(wall)
                cpu_vals = [r["cpu_time_s"] for r in ok if r["cpu_time_s"] is not None]
                cpu_median = statistics.median(cpu_vals) if cpu_vals else None
                cpu_min = min(cpu_vals) if cpu_vals else None
                cpu_max = max(cpu_vals) if cpu_vals else None
                speedup = None
                efficiency = None
                baseline = baseline_map.get(ions)
                if baseline and wall_median > 0:
                    speedup = baseline / wall_median
                    efficiency = (speedup / threads) * 100.0
            else:
                wall_median = wall_min = wall_max = None
                cpu_median = cpu_min = cpu_max = None
                speedup = efficiency = None

            summary_rows.append({
                "ions": ions,
                "threads": threads,
                "repeat_count": len(ok),
                "fail_count": fail_count,
                "median_wall_time_s": wall_median,
                "min_wall_time_s": wall_min,
                "max_wall_time_s": wall_max,
                "median_cpu_time_s": cpu_median,
                "min_cpu_time_s": cpu_min,
                "max_cpu_time_s": cpu_max,
                "speedup": speedup,
                "efficiency_percent": efficiency,
            })

    with summary_path.open("w", encoding="utf-8", newline="") as summary_file:
        writer = csv.writer(summary_file)
        writer.writerow([
            "ions", "threads", "repeat_count", "fail_count",
            "median_wall_time_s", "min_wall_time_s", "max_wall_time_s",
            "median_cpu_time_s", "min_cpu_time_s", "max_cpu_time_s",
            "speedup", "efficiency_percent",
        ])
        for row in summary_rows:
            writer.writerow([
                row["ions"],
                row["threads"],
                row["repeat_count"],
                row["fail_count"],
                "" if row["median_wall_time_s"] is None else f"{row['median_wall_time_s']:.3f}",
                "" if row["min_wall_time_s"] is None else f"{row['min_wall_time_s']:.3f}",
                "" if row["max_wall_time_s"] is None else f"{row['max_wall_time_s']:.3f}",
                "" if row["median_cpu_time_s"] is None else f"{row['median_cpu_time_s']:.3f}",
                "" if row["min_cpu_time_s"] is None else f"{row['min_cpu_time_s']:.3f}",
                "" if row["max_cpu_time_s"] is None else f"{row['max_cpu_time_s']:.3f}",
                "" if row["speedup"] is None else f"{row['speedup']:.3f}",
                "" if row["efficiency_percent"] is None else f"{row['efficiency_percent']:.1f}",
            ])

    print(f"Raw results: {csv_path}")
    print(f"Summary: {summary_path}")
    print(f"Metadata: {metadata_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

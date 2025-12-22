#!/usr/bin/env python3
"""
Generate IMS configs sweeping E/N (Td) based on a template ims_basic.json.

Default grid: 5 Td .. 40 Td in 5 Td steps.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any, Dict, List


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Generate IMS E/N sweep configs.")
    p.add_argument("--template", type=Path, default=Path("../ICARION_Inputs/ims_basic.json"),
                   help="Template config path (ims_basic.json).")
    p.add_argument("--output-dir", type=Path, default=Path("tmp/ims_en_sweep/configs"),
                   help="Directory for generated configs.")
    p.add_argument("--en-start", type=float, default=5.0, help="Start E/N [Td].")
    p.add_argument("--en-stop", type=float, default=100.0, help="Stop E/N [Td].")
    p.add_argument("--en-step", type=float, default=5.0, help="Step E/N [Td].")
    p.add_argument("--output-folder", type=Path, default=Path("results/ims_en_sweep"),
                   help="Output folder for trajectories.")
    p.add_argument("--species-db", type=Path, default=Path("data/species_database_v1.json"),
                   help="Species database path (written as absolute path).")
    p.add_argument("--domain-name", type=str, default="drift_region",
                   help="Domain name to apply EN_Td to.")
    p.add_argument("--run-script", type=Path, default=None,
                   help="Optional run script to execute all configs.")
    p.add_argument("--binary", type=str, default="./build/src/icarion_main",
                   help="Binary to call in the run script.")
    return p.parse_args()


def load_template(path: Path) -> Dict[str, Any]:
    with path.open() as fh:
        return json.load(fh)


def find_domain(domains: List[Dict[str, Any]], name: str) -> Dict[str, Any]:
    for dom in domains:
        if dom.get("name") == name:
            return dom
    raise KeyError(f"Domain '{name}' not found in template.")


def generate_values(start: float, stop: float, step: float) -> List[float]:
    values = []
    v = start
    while v <= stop + 1e-9:
        values.append(round(v, 6))
        v += step
    return values


def write_run_script(paths: List[Path], binary: str, script_path: Path) -> None:
    script_path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["#!/usr/bin/env bash", "set -euo pipefail", ""]
    for cfg in paths:
        lines.append(f'{binary} "{cfg}"')
    script_path.write_text("\n".join(lines) + "\n")
    script_path.chmod(0o755)


def main() -> None:
    args = parse_args()
    cfg_template = load_template(args.template)
    cfg_template["species_database"] = str(args.species_db.resolve())
    cfg_template["output"]["folder"] = str(args.output_folder.resolve())

    values = generate_values(args.en_start, args.en_stop, args.en_step)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    written: List[Path] = []

    for en_td in values:
        cfg = json.loads(json.dumps(cfg_template))
        domain = find_domain(cfg["domains"], args.domain_name)
        domain.setdefault("fields", {}).setdefault("DC", {})["EN_Td"] = en_td

        cfg["output"]["trajectory_file"] = f"ims_en_{en_td:g}Td.h5"
        out_path = args.output_dir / f"ims_en_{en_td:g}Td.json"
        with out_path.open("w") as fh:
            json.dump(cfg, fh, indent=2)
        written.append(out_path)

    print(f"Wrote {len(written)} configs to {args.output_dir}")

    if args.run_script:
        write_run_script(written, args.binary, args.run_script)
        print(f"Wrote run script: {args.run_script}")


if __name__ == "__main__":
    main()

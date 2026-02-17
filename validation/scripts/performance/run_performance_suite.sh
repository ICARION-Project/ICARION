#!/bin/bash
# Unified ICARION performance benchmark runner (CPU-only in v1.0.0; GPU runtime-disabled)

set -euo pipefail

ORIGINAL_ARGS=("$@")

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(cd "$SCRIPT_DIR/../../.." && pwd)"
BUILD_DIR="$PROJECT_ROOT/build"
CPU_CONFIG_DIR="$PROJECT_ROOT/validation/configs/performance"
GPU_CONFIG_DIR="$PROJECT_ROOT/validation/configs/performance/gpu"
CPU_RESULT_DIR="$PROJECT_ROOT/validation/results/v1.0.0_test/performance/logs"
GPU_RESULT_DIR="$PROJECT_ROOT/validation/results/v1.0.0_test/performance/gpu_logs"
ICARION_BIN="$BUILD_DIR/src/icarion_main"
CMAKE_CACHE="$BUILD_DIR/CMakeCache.txt"

RUN_ID=${RUN_ID:-""}
RUN_DIR=${RUN_DIR:-""}
BASELINE_OUTPUT=false

RUN_ID_SET=false
RUN_DIR_SET=false

declare -A CPU_CATEGORY_GLOBS=(
  [baseline]="scaling_baseline_*.json"
  [collision]="collision_overhead_*.json"
  [spacecharge]="spacecharge_overhead_*.json"
  [longrun]="scaling_longrun_*.json"
  [mixed]="scaling_mixedphysics_*.json"
)
declare -a CPU_DEFAULT_CATEGORIES=(baseline collision spacecharge longrun mixed)

declare -A GPU_CATEGORY_GLOBS=(
  [cpu_scaling]="RK4_cpu_*.json"
  [gpu_scaling]="RK4_gpu_*.json"
  [integrators]="integrator_*.json"
  [threshold]="threshold_*.json"
  [long]="long_*.json"
)
declare -a GPU_DEFAULT_CATEGORIES=(cpu_scaling gpu_scaling integrators threshold long)

print_usage() {
  cat <<EOF
Usage: $(basename "$0") [options] [category ...]

Options:
  --cpu-only        Run only CPU benchmark categories (default)
  --gpu-only        (ignored in v1.0.0) GPU runtime is disabled; GPU categories are skipped
  --run-id ID        Run identifier (default: YYYYmmdd_HHMMSS)
  --run-dir PATH     Output directory for this run (default: validation/runs/<run-id>)
  --baseline-output  Write into validation/results/v1.0.0_test/performance (legacy)
  -h, --help        Show this help and exit

Categories:
  CPU: ${CPU_DEFAULT_CATEGORIES[*]}
  GPU: ${GPU_DEFAULT_CATEGORIES[*]} (retained for future releases; skipped in v1.0.0)
  Special aliases:
    all       Run all CPU categories (GPU section skipped in v1.0.0)
    cpu_all   Run all CPU categories
    gpu_all   (ignored in v1.0.0) Run GPU categories when backend is re-enabled
EOF
}

dedupe_array() {
  local -n _arr=$1
  local -A _seen=()
  local _unique=()
  for item in "${_arr[@]:-}"; do
    if [[ -n "$item" && -z "${_seen[$item]:-}" ]]; then
      _unique+=("$item")
      _seen[$item]=1
    fi
  done
  _arr=("${_unique[@]:-}")
}

MODE="cpu" # all|cpu|gpu (default to CPU-only; GPU runtime-disabled in v1.0.0)
USER_CATEGORIES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --cpu-only)
      MODE="cpu"
      shift
      ;;
    --gpu-only)
      MODE="gpu"
      shift
      ;;
    --run-id)
      RUN_ID="$2"
      RUN_ID_SET=true
      shift 2
      ;;
    --run-dir)
      RUN_DIR="$2"
      RUN_DIR_SET=true
      shift 2
      ;;
    --baseline-output)
      BASELINE_OUTPUT=true
      shift
      ;;
    -h|--help)
      print_usage
      exit 0
      ;;
    *)
      USER_CATEGORIES+=("$1")
      shift
      ;;
  esac
done

if [[ -z "$RUN_ID" ]]; then
  RUN_ID=$(date +%Y%m%d_%H%M%S)
fi

RUN_DIR_ABS=""
if ! $BASELINE_OUTPUT || $RUN_ID_SET || $RUN_DIR_SET; then
  if [[ -z "$RUN_DIR" ]]; then
    RUN_DIR="$PROJECT_ROOT/validation/runs/$RUN_ID"
  fi

  RUN_DIR_ABS=$(cd "$(dirname "$RUN_DIR")" && pwd)/"$(basename "$RUN_DIR")"
  mkdir -p "$RUN_DIR_ABS/logs" "$RUN_DIR_ABS/figures" "$RUN_DIR_ABS/results"
  export ICARION_VALIDATION_RUN_DIR="$RUN_DIR_ABS"

  if ! $BASELINE_OUTPUT; then
    CPU_RESULT_DIR="$RUN_DIR_ABS/results/performance/logs"
    GPU_RESULT_DIR="$RUN_DIR_ABS/results/performance/gpu_logs"
  fi

  if command -v python3 >/dev/null 2>&1; then
    cmd=$(printf '%q ' "$0" "${ORIGINAL_ARGS[@]}")
    SUITE_NAME="performance" \
    RUN_ID="$RUN_ID" \
    RUN_DIR="$RUN_DIR_ABS" \
    REPO_ROOT="$PROJECT_ROOT" \
    VALIDATION_DIR="$PROJECT_ROOT/validation" \
    ICARION_BIN="$ICARION_BIN" \
    COMMAND_LINE="$cmd" \
    python3 - "$RUN_DIR_ABS/manifest.performance.json" <<'PY'
import json
import os
import subprocess
import sys
from datetime import datetime, timezone

out_file = sys.argv[1]

def _cmd(argv):
    try:
        return subprocess.check_output(argv, stderr=subprocess.DEVNULL, text=True).strip()
    except Exception:
        return None

repo_root = os.environ.get("REPO_ROOT")

data = {
    "suite": os.environ.get("SUITE_NAME"),
    "run_id": os.environ.get("RUN_ID"),
    "run_dir": os.environ.get("RUN_DIR"),
    "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    "repo_root": repo_root,
    "validation_dir": os.environ.get("VALIDATION_DIR"),
    "command": os.environ.get("COMMAND_LINE"),
    "icarion_bin": os.environ.get("ICARION_BIN") or None,
}

data["git_commit"] = _cmd(["git", "-C", repo_root, "rev-parse", "HEAD"]) if repo_root else None
git_status = _cmd(["git", "-C", repo_root, "status", "--porcelain"]) if repo_root else None
data["git_dirty"] = bool(git_status)

os.makedirs(os.path.dirname(out_file), exist_ok=True)
with open(out_file, "w", encoding="utf-8") as handle:
    json.dump(data, handle, indent=2)
PY
  fi
fi

declare -a CPU_SELECTED=()
declare -a GPU_SELECTED=()

if [ ${#USER_CATEGORIES[@]} -eq 0 ]; then
  CPU_SELECTED=("${CPU_DEFAULT_CATEGORIES[@]}")
  GPU_SELECTED=()
else
  for category in "${USER_CATEGORIES[@]}"; do
    case "$category" in
      all)
        CPU_SELECTED+=("${CPU_DEFAULT_CATEGORIES[@]}")
        ;;
      cpu_all)
        CPU_SELECTED+=("${CPU_DEFAULT_CATEGORIES[@]}")
        ;;
      gpu_all)
        GPU_SELECTED+=("${GPU_DEFAULT_CATEGORIES[@]}")
        ;;
      *)
        if [[ -n "${CPU_CATEGORY_GLOBS[$category]:-}" ]]; then
          CPU_SELECTED+=("$category")
        elif [[ -n "${GPU_CATEGORY_GLOBS[$category]:-}" ]]; then
          GPU_SELECTED+=("$category")
        else
          echo "Unknown category: $category" >&2
          print_usage
          exit 1
        fi
        ;;
    esac
  done
fi

if [[ "$MODE" == "cpu" ]]; then
  GPU_SELECTED=()
elif [[ "$MODE" == "gpu" ]]; then
  CPU_SELECTED=()
fi

dedupe_array CPU_SELECTED
dedupe_array GPU_SELECTED

# GPU runtime is disabled in v1.0.0; honor CPU selections only and skip GPU categories.
if [[ "$MODE" == "gpu" ]] || [[ ${#GPU_SELECTED[@]} -gt 0 ]]; then
  echo "GPU benchmarks are runtime-disabled for v1.0.0 (enable_gpu falls back to CPU). Skipping GPU categories." >&2
  if [[ "$MODE" == "gpu" ]]; then
    exit 0
  fi
  GPU_SELECTED=()
fi

if [ ! -x "$ICARION_BIN" ]; then
  echo "Error: ICARION binary not found at $ICARION_BIN" >&2
  exit 1
fi

GPU_ENABLED=0
if [ -f "$CMAKE_CACHE" ] && grep -q '^USE_GPU_ACCEL:BOOL=ON' "$CMAKE_CACHE"; then
  GPU_ENABLED=1
fi

if [ ${#GPU_SELECTED[@]} -gt 0 ] && [ $GPU_ENABLED -ne 1 ]; then
  if [[ "$MODE" == "gpu" ]]; then
    echo "Error: GPU benchmarks requested but this build was configured without USE_GPU_ACCEL=ON" >&2
    exit 1
  fi
  echo "Skipping GPU benchmarks because USE_GPU_ACCEL=OFF in $CMAKE_CACHE"
  GPU_SELECTED=()
fi

if [ ${#CPU_SELECTED[@]} -gt 0 ] && [ ! -d "$CPU_CONFIG_DIR" ]; then
  echo "Error: CPU config directory not found at $CPU_CONFIG_DIR" >&2
  exit 1
fi

if [ ${#GPU_SELECTED[@]} -gt 0 ] && [ ! -d "$GPU_CONFIG_DIR" ]; then
  echo "Error: GPU config directory not found at $GPU_CONFIG_DIR" >&2
  exit 1
fi

mkdir -p "$CPU_RESULT_DIR"
mkdir -p "$GPU_RESULT_DIR"

echo ""
echo "============================================================"
echo "ICARION Performance Benchmark Suite"
echo "============================================================"
echo "Binary:  $ICARION_BIN"
if [[ -n "$RUN_DIR_ABS" ]]; then
  echo "Run dir: $RUN_DIR_ABS"
else
  echo "Run dir: (baseline-output)"
fi
echo "CPU configs: $CPU_CONFIG_DIR"
echo "GPU configs: $GPU_CONFIG_DIR"
echo "============================================================"
echo ""

TOTAL_RUNS=0
CPU_RUNS=0
GPU_RUNS=0

run_config() {
  local config_path="$1"
  local result_dir="$2"
  local name="$(basename "$config_path" .json)"
  local log_file="$result_dir/$name.log"
  local time_file="$result_dir/$name.time"

  local effective_config_path="$config_path"
  if [[ -n "${RUN_DIR_ABS:-}" ]] && command -v python3 >/dev/null 2>&1; then
    local snapshot_dir="$RUN_DIR_ABS/results/performance/config_snapshots"
    local sim_output_dir="$RUN_DIR_ABS/results/performance/sim_output/$name"
    mkdir -p "$snapshot_dir" "$sim_output_dir"
    local snapshot_path="$snapshot_dir/$name.run_config.json"

    PROJECT_ROOT="$PROJECT_ROOT" \
    OUT_DIR="$sim_output_dir" \
    python3 - "$config_path" "$snapshot_path" <<'PY'
import json
import os
import sys
from pathlib import Path

src = Path(sys.argv[1])
dst = Path(sys.argv[2])
project_root = Path(os.environ["PROJECT_ROOT"])  # repo root
out_dir = Path(os.environ["OUT_DIR"])            # per-config sim output dir

with src.open("r", encoding="utf-8") as handle:
    cfg = json.load(handle)

output = cfg.get("output") or {}
output["folder"] = str(out_dir)
cfg["output"] = output

for key in ("species_database_path", "species_database", "reaction_database", "reaction_database_path"):
    if key in cfg and isinstance(cfg[key], str):
        p = Path(cfg[key])
        if not p.is_absolute():
            cfg[key] = str((project_root / p).resolve())

dst.parent.mkdir(parents=True, exist_ok=True)
with dst.open("w", encoding="utf-8") as handle:
    json.dump(cfg, handle, indent=2)
    handle.write("\n")
PY

    effective_config_path="$snapshot_path"
  fi

  echo "→ Running $name"
  if (cd "$PROJECT_ROOT" && /usr/bin/time -f "%e" -o "$time_file" "$ICARION_BIN" "$effective_config_path" > "$log_file" 2>&1); then
    local elapsed
    elapsed="$(cat "$time_file")"
    echo "  ✓ Completed (${elapsed}s)"
  else
    local status=$?
    echo "  ✗ Failed (exit $status)"
    echo "    Log: $log_file"
  fi
}

run_category() {
  local config_dir="$1"
  local result_dir="$2"
  local category="$3"
  local pattern="$4"
  local kind="$5"

  echo ""
  echo "========================================"
  echo "Category: $category"
  echo "========================================"
  echo ""

  local matches=0
  shopt -s nullglob
  for config in "$config_dir"/$pattern; do
    matches=$((matches + 1))
    run_config "$config" "$result_dir"
  done
  shopt -u nullglob

  if [ $matches -eq 0 ]; then
    echo "(no configs matched pattern '$pattern')"
  else
    echo ""
    echo "Completed $matches runs for '$category'"
  fi

  TOTAL_RUNS=$((TOTAL_RUNS + matches))
  if [[ "$kind" == "cpu" ]]; then
    CPU_RUNS=$((CPU_RUNS + matches))
  else
    GPU_RUNS=$((GPU_RUNS + matches))
  fi
}

for category in "${CPU_SELECTED[@]}"; do
  [[ -n "$category" ]] || continue
  pattern=${CPU_CATEGORY_GLOBS[$category]:-}
  if [ -n "$pattern" ]; then
    run_category "$CPU_CONFIG_DIR" "$CPU_RESULT_DIR" "$category" "$pattern" "cpu"
  fi
done

if [ ${#GPU_SELECTED[@]} -gt 0 ]; then
  echo ""
  echo "GPU benchmarks enabled (USE_GPU_ACCEL=ON)"
fi

for category in "${GPU_SELECTED[@]}"; do
  [[ -n "$category" ]] || continue
  pattern=${GPU_CATEGORY_GLOBS[$category]:-}
  if [ -n "$pattern" ]; then
    run_category "$GPU_CONFIG_DIR" "$GPU_RESULT_DIR" "$category" "$pattern" "gpu"
  fi
done

if [ $TOTAL_RUNS -eq 0 ]; then
  echo "No benchmarks were executed."
  exit 0
fi

echo ""
echo "============================================================"
echo "Summary"
echo "============================================================"
echo "CPU runs: $CPU_RUNS (logs: $CPU_RESULT_DIR)"
if [ $GPU_RUNS -gt 0 ]; then
  echo "GPU runs: $GPU_RUNS (logs: $GPU_RESULT_DIR)"
else
  echo "GPU runs: 0"
fi
echo "Total runs: $TOTAL_RUNS"
echo "============================================================"
echo ""
